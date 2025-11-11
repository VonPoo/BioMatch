#!/usr/bin/env python3
"""
BioMatch Analysis Script: filterRepetiveSNP

Purpose
- Filter and stratify unique k-mers from a SAM file according to uniqueness
  tags (`X0`, `X1`) to reproduce the original Perl script behavior.
- Output stratified FASTA files named `<prefix>_n<i>.fa` where `i` is the
  uniqueness bucket index.

Inputs
- `sam_file`: SAM file from the previous alignment step.
- `prefix`: Output file prefix.
- `kmer_size`: K-mer window size.
- `subkmer_size`: Sub k-mer size; `max_count = kmer_size - subkmer_size + 1`.

Outputs
- Multiple FASTA files: `<prefix>_n0.fa`, `<prefix>_n1.fa`, ..., one per
  uniqueness bucket, each containing paired sequences (`ref` and `var`).

Usage Example
    python filterRepetiveSNP.py \
        input.sam \
        output_prefix \
        31 \
        19

Notes
- A k-mer is considered unique if `X0` tag is missing (treated unique), or
  `X0 + X1 == 1` when present.
- The script writes informative progress to STDERR and handles malformed tags.
"""
import sys
import argparse
import re
import os
from collections import defaultdict

def parse_sam_tags(tags_list):
    """Parse SAM tags to determine k-mer uniqueness.

    A k-mer is unique if `X0` is missing, or `X0 + X1 == 1` when present.
    Malformed tag values are ignored.
    """
    x0 = None
    x1 = 0
    
    for tag in tags_list:
        if tag.startswith("X0:i:"):
            try:
                x0 = int(tag[5:])
            except ValueError:
                pass # Malformed tag, ignore
        elif tag.startswith("X1:i:"):
            try:
                x1 = int(tag[5:])
            except ValueError:
                pass # Malformed tag, ignore
    
    # 1. Corresponds to the Perl script's 'else' branch (X0 tag does not exist)
    if x0 is None:
        return True
    
    # 2. Corresponds to the Perl script's 'if' branch (X0 tag exists)
    total_count = x0 + x1
    if total_count == 1:
        return True
    
    return False

def main():
    parser = argparse.ArgumentParser(
        description="Filter and stratify unique k-mers from a SAM file.",
        epilog="This script is a Python refactor of the original Perl script."
    )
    parser.add_argument("sam_file", help="Input SAM file (from the previous alignment step)")
    parser.add_argument("prefix", help="Output file prefix (e.g., 'output')")
    parser.add_argument("kmer_size", type=int, help="K-mer window size (e.g., 31)")
    parser.add_argument("subkmer_size", type=int, help="Sub k-mer size (e.g., 19)")
    
    args = parser.parse_args()

    try:
        max_count = args.kmer_size - args.subkmer_size + 1
    except ValueError:
        print("Error: kmer_size and subkmer_size must be integers", file=sys.stderr)
        sys.exit(1)
        
    if max_count <= 0:
        print(f"Error: kmer_size ({args.kmer_size}) must be greater than or equal to subkmer_size ({args.subkmer_size})", file=sys.stderr)
        sys.exit(1)

    # Compile QNAME regex for reuse
    # Format: id|pos|type (e.g.: rs4970383|0|AT)
    qname_regex = re.compile(r'([^\|]+)\|(\d+)\|(AT|CG)')

    # Equivalent to Perl's %idToUniqCount and %idToStr
    # Using defaultdict(dict) would be cleaner, but we manually create nested dicts
    # to maintain a 1:1 logical mapping to the Perl script.
    id_to_uniq_count = {}
    id_to_str = {}

    print(f"Starting to parse SAM file: {args.sam_file}", file=sys.stderr)
    print(f"Maximum k-mer count (maxCount) = {max_count}", file=sys.stderr)

    # --- Phase 1: Parse the SAM file ---
    with open(args.sam_file, 'r') as fh:
        for line in fh:
            # Skip SAM header
            if line.startswith('@'):
                continue
            
            try:
                tmpAr = line.strip().split('\t')
                
                # Ensure there are enough columns to parse
                if len(tmpAr) < 11:
                    continue 

                qname = tmpAr[0]
                seq = tmpAr[9]
                tags = tmpAr[11:]

                match = qname_regex.match(qname)
                if not match:
                    print(f"Warning: Could not parse QNAME: {qname}", file=sys.stderr)
                    continue
                
                snp_id, pos, allele_type = match.groups()

                # Initialize data structures for this SNP
                if snp_id not in id_to_uniq_count:
                    id_to_uniq_count[snp_id] = {}
                if snp_id not in id_to_str:
                    id_to_str[snp_id] = {}
                
                # Initialize the counter (equivalent to Perl's 'if (!exists($idToUniqCount{$id}{$type}))')
                if allele_type not in id_to_uniq_count[snp_id]:
                    id_to_uniq_count[snp_id][allele_type] = max_count
                
                # Check if the k-mer is unique
                is_unique = parse_sam_tags(tags)
                
                if is_unique:
                    # Decrement the counter
                    id_to_uniq_count[snp_id][allele_type] -= 1
                    
                    # Concatenate the sequence
                    # Note: We are reproducing the *active* logic from the Perl script,
                    # not the commented-out 'currPos' logic.
                    if allele_type in id_to_str[snp_id]:
                        id_to_str[snp_id][allele_type] += 'N' + seq
                    else:
                        id_to_str[snp_id][allele_type] = seq

            except Exception as e:
                print(f"Error processing line: {e}\nLine content: {line.strip()}", file=sys.stderr)

    print("SAM file parsing complete. Writing stratified files...", file=sys.stderr)

    # --- Phase 2: Stratified Output ---
    fhs = [] # List of file handles
    try:
        # Pre-create all files
        for i in range(max_count):
            filename = f"{args.prefix}_n{i}.fa"
            fhs.append(open(filename, 'w'))
            
        # Iterate over all SNPs
        for snp_id in sorted(id_to_uniq_count.keys()):
            
            # Check for each stratification level i
            for i in range(max_count):
                
                count_at = id_to_uniq_count[snp_id].get("AT")
                count_cg = id_to_uniq_count[snp_id].get("CG")
                
                # Check if both AT and CG alleles are present
                if count_at is not None and count_cg is not None:
                    
                    # Check if the threshold for level i is met
                    if count_at <= i and count_cg <= i:
                        
                        seq_at = id_to_str[snp_id].get("AT")
                        seq_cg = id_to_str[snp_id].get("CG")
                        
                        # Check if the sequences also exist
                        if seq_at is not None and seq_cg is not None:
                            fhs[i].write(f">{snp_id} ref\n{seq_at}\n")
                            fhs[i].write(f">{snp_id} var\n{seq_cg}\n")
                        else:
                            # Replicate the error message from the Perl script
                            exists_at = 1 if seq_at is not None else 0
                            exists_cg = 1 if seq_cg is not None else 0
                            print(f"Possible file truncation. Missing: "
                                  f"{snp_id} {i} {count_at} {count_cg} {exists_at} {exists_cg}", 
                                  file=sys.stderr)
                            
    except IOError as e:
        print(f"File writing error: {e}", file=sys.stderr)
    finally:
        # Ensure all files are closed
        for fh in fhs:
            fh.close()
            
    print("Processing complete.", file=sys.stderr)

if __name__ == "__main__":
    main()
