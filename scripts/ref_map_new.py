#!/usr/bin/env python3
"""
BioMatch Analysis Script: ref_map_new

Purpose
- Rename FASTA sequence headers to match chromosome naming expected by VCF /
  downstream tools. Supports multiple extraction modes and duplicate handling.

Key Features
- Modes: `auto`, `chromosome`, `template`, `field` to extract chromosome names
  from FASTA headers.
- Duplicate handling strategies: `primary` (keep longest per chromosome) and
  `unique` (append unique suffix for duplicates).
- Optional VCF validation to ensure renamed chromosomes align with VCF records.

Inputs
- `-ref/--reference`: Input reference genome (FASTA).
- `-out/--output`: Output path for renamed reference FASTA.
- Optional: `-vcf/--vcf-file` for validation.
- Mode-related options: `-mode`, `-template`, `-field`, `-prefix`.

Outputs
- Renamed FASTA file with standardized chromosome names. Logs are printed to
  STDOUT/STDERR depending on `--quiet/--verbose`.

Usage Examples
    # Auto mode (recommended)
    python ref_map_new.py -ref genome.fasta -out genome_renamed.fasta -mode auto -prefix chr

    # Validate with VCF
    python ref_map_new.py -ref genome.fasta -out genome_renamed.fasta -vcf variants.vcf

    # Template mode
    python ref_map_new.py -ref genome.fasta -out genome_renamed.fasta \
        -mode template -template "NC_\\d+\\.\\d+=>\\d+"

Notes
- `auto` mode samples headers and selects the best heuristic with confidence.
- `chromosome` mode detects patterns like `chromosome N`, `chr N`, or `chrN`.
- `template` mode extracts via `ID_PATTERN=>CHROM_PATTERN` pairing.
- `field` mode uses a specified header token.
"""

import argparse
import re
import os
import sys
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description='Rename FASTA headers to match VCF chromosome naming',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            """
Example:
  # Auto mode (recommended)
  %(prog)s -ref genome.fasta -out genome_renamed.fasta -mode auto -prefix chr

  # Validate with VCF
  %(prog)s -ref genome.fasta -out genome_renamed.fasta -vcf variants.vcf

  # Template mode
  %(prog)s -ref genome.fasta -out genome_renamed.fasta -mode template -template "NC_\\d+\\.\\d+=>\\d+"
            """
        ),
    )
    parser.add_argument('-ref', '--reference', required=True, help='Input reference genome (FASTA)')
    parser.add_argument('-vcf', '--vcf-file', help='Optional VCF file for validation')
    parser.add_argument('-out', '--output', required=True, help='Output path for renamed reference genome')
    parser.add_argument('-strategy', '--strategy', choices=['primary', 'unique'], default='primary',
                      help='Duplicate handling: primary=keep longest per chromosome; unique=add unique suffix for duplicates')
    parser.add_argument('-mode', '--mode', choices=['auto', 'chromosome', 'template', 'field'], default='auto',
                      help='Extraction mode: auto (recommended), chromosome/chr, template, or field index')
    parser.add_argument('-template', '--template',
                      help='Extraction template "ID_PATTERN=>CHROM_PATTERN", e.g., "NC_\\d+\\.\\d+=>\\d+"')
    parser.add_argument('-field', '--field', type=int, default=1,
                      help='Header field index (0-based; negative counts from right)')
    parser.add_argument('-clean', '--clean', action='store_true', default=True,
                      help='Clean extracted chromosome names (remove trailing punctuation)')
    parser.add_argument('-no-clean', '--no-clean', action='store_true',
                      help='Disable cleaning of extracted chromosome names')
    parser.add_argument('-prefix', '--prefix', default='', help='Prefix for output chromosome names (e.g., chr or auto)')
    parser.add_argument('-quiet', '--quiet', action='store_true', help='Suppress non-essential log output')
    parser.add_argument('-verbose', '--verbose', action='store_true', help='Show detailed processing log')
    return parser.parse_args()

def clean_chrom_name(chrom):
    """Sanitize chromosome name by removing punctuation"""
    chrom = re.sub(r'[,;:\s]+$', '', chrom)
    chrom = re.sub(r'^[,;:\s]+', '', chrom)
    chrom = chrom.strip('\"\'')
    return chrom

def _log(msg, quiet=False):
    if not quiet:
        print(msg)

def _extract_chrom_by_chromosome_mode(header):
    patterns = [
        r'chromosome[ _\s]+(\d+|X|Y|MT|M)(?:[,\s\)]|$)',
        r'chr[ _\s]+(\d+|X|Y|MT|M)(?:[,\s\)]|$)',
        r'\bchr(\d+|X|Y|MT|M)(?:[,\s\)]|$)'
    ]
    for p in patterns:
        m = re.search(p, header, re.IGNORECASE)
        if m:
            return m.group(1)
    return None

def _extract_chrom_by_template(header, seq_id, id_regex, chrom_regex):
    try:
        id_match = id_regex.search(header)
        chrom_match = chrom_regex.search(header)
        if id_match and chrom_match:
            extracted_id = id_match.group(0)
            if extracted_id in seq_id:
                return chrom_match.group(0)
    except Exception:
        pass
    return None

def _extract_chrom_by_field(header, field_index):
    fields = header[1:].split()
    if not fields:
        return None
    try:
        if field_index >= 0 and field_index < len(fields):
            return fields[field_index]
        elif field_index < 0 and abs(field_index) <= len(fields):
            return fields[field_index]
    except IndexError:
        return None
    return None

def _try_common_patterns(header, seq_id):
    chrom = _extract_chrom_by_chromosome_mode(header)
    if chrom:
        return chrom, "chromosome/chr"
    for idx in [0, 1, 2, -1]:
        chrom = _extract_chrom_by_field(header, idx)
        if chrom and re.match(r'^(\d+|X|Y|MT|M|chr\d+|chr[XYM])', chrom, re.IGNORECASE):
            return chrom, f"field {idx}"
    acc = re.search(r'NC_\d+\.\d+', header)
    num = re.search(r'\b(\d+|X|Y|MT|M)\b', header)
    if acc and num and acc.group(0) in seq_id:
        return num.group(1), "NC_=>num"
    return None, None

def _auto_detect_mode(reference_file, sample_size=100, clean=True, quiet=False):
    samples = []
    results = defaultdict(list)
    try:
        with open(reference_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()
                    seq_id = header[1:].split()[0] if len(header) > 1 else ''
                    if seq_id:
                        samples.append((seq_id, header))
                        if len(samples) >= sample_size:
                            break
    except Exception as e:
        print(f"Error sampling reference: {e}", file=sys.stderr)
        return None, {}, 0.0
    if not samples:
        return None, {}, 0.0
    for seq_id, header in samples:
        chrom, method = _try_common_patterns(header, seq_id)
        if chrom:
            if clean:
                chrom = clean_chrom_name(chrom)
            if chrom:
                results[method].append((seq_id, chrom))
    best_method, best = None, []
    for m, vals in results.items():
        if len(vals) > len(best):
            best_method, best = m, vals
    if not best_method:
        return None, {}, 0.0
    success = len(best) / len(samples)
    if best_method.startswith('field'):
        field_idx = int(best_method.split()[-1])
        return 'field', {'field': field_idx}, success
    if best_method == 'chromosome/chr':
        return 'chromosome', {}, success
    if best_method == 'NC_=>num':
        return 'template', {'template': r'NC_\d+\.\d+=>\d+|X|Y|MT|M'}, success
    return 'chromosome', {}, success

def extract_chromosome_mapping(reference_file, mode, template=None, field=None, clean=True, quiet=False):
    """Extract mapping from sequence IDs to chromosome names from FASTA headers according to mode"""
    mapping = {}
    seq_lengths = {}
    id_regex = None
    chrom_regex = None
    if mode == 'auto':
        _log("Auto-detecting extraction mode...", quiet)
        detected, params, rate = _auto_detect_mode(reference_file, clean=clean, quiet=quiet)
        if not detected or rate < 0.5:
            _log("Auto-detection low confidence; fallback to chromosome mode", quiet)
            mode = 'chromosome'
        else:
            mode = detected
            field = params.get('field', field)
            template = params.get('template', template)
            _log(f"Selected mode: {mode} (confidence: {rate:.1%})", quiet)
    if mode == 'template' and template:
        try:
            id_pattern, chrom_pattern = template.split('=>')
            id_regex = re.compile(id_pattern)
            chrom_regex = re.compile(chrom_pattern)
        except Exception:
            print(f"Error: Invalid template format: {template}. Expected 'ID_PATTERN=>CHROM_PATTERN'", file=sys.stderr)
            sys.exit(1)
    with open(reference_file, 'r') as f:
        current_seq_id = None
        current_length = 0
        for line in f:
            if line.startswith('>'):
                if current_seq_id is not None:
                    seq_lengths[current_seq_id] = current_length
                header = line.strip()
                current_seq_id = header[1:].split()[0]
                current_length = 0
                seq_id = header[1:].split()[0]
                chrom_num = None
                if mode == 'chromosome':
                    chrom_num = _extract_chrom_by_chromosome_mode(header)
                elif mode == 'template':
                    chrom_num = _extract_chrom_by_template(header, seq_id, id_regex, chrom_regex)
                elif mode == 'field':
                    chrom_num = _extract_chrom_by_field(header, field if field is not None else 0)
                if chrom_num and clean:
                    chrom_num = clean_chrom_name(chrom_num)
                if chrom_num:
                    mapping[seq_id] = chrom_num
            else:
                current_length += len(line.strip())
        if current_seq_id is not None:
            seq_lengths[current_seq_id] = current_length
    return mapping, seq_lengths

def find_primary_sequences(mapping, seq_lengths):
    chrom_to_seqs = defaultdict(list)
    for seq_id, chrom in mapping.items():
        if seq_id in seq_lengths:
            chrom_to_seqs[chrom].append((seq_id, seq_lengths[seq_id]))
    primary_seqs = {}
    for chrom, seqs in chrom_to_seqs.items():
        seqs.sort(key=lambda x: x[1], reverse=True)
        primary_seqs[chrom] = seqs[0][0]
    return primary_seqs

def validate_vcf_chromosomes(vcf_file, mapping):
    if not vcf_file or not os.path.exists(vcf_file):
        return
    vcf_chroms = set()
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                chrom = line.split('\t')[0]
                vcf_chroms.add(chrom)
    reversed_mapping = set(mapping.values())
    missing_chroms = []
    for chrom in vcf_chroms:
        if chrom not in reversed_mapping:
            missing_chroms.append(chrom)
    if missing_chroms:
        print(f"Note: {len(missing_chroms)} chromosomes in VCF not found in reference mapping", file=sys.stderr)
        print(f"Missing chromosomes: {', '.join(missing_chroms[:5])}{'...' if len(missing_chroms) > 5 else ''}", file=sys.stderr)
    return vcf_chroms

def detect_vcf_prefix(vcf_file, quiet=False):
    """Detect whether VCF uses 'chr' prefix. Returns 'chr' or ''.

    Heuristic: if >=60% of chromosome names start with 'chr', use 'chr'; else ''.
    """
    if not vcf_file or not os.path.exists(vcf_file):
        return ''
    total = 0
    chr_prefixed = 0
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if not parts:
                    continue
                total += 1
                chrom = parts[0].strip()
                if chrom.lower().startswith('chr'):
                    chr_prefixed += 1
                # stop early if enough lines
                if total >= 2000:
                    break
    except Exception:
        return ''
    if total == 0:
        return ''
    ratio = chr_prefixed / total
    return 'chr' if ratio >= 0.6 else ''

def _apply_prefix_strategy(chrom: str, prefix: str) -> str:
    """Normalize chromosome naming based on desired prefix.

    - If prefix == 'chr' and chrom doesn't start with 'chr', add it.
    - If prefix == '' and chrom starts with 'chr', strip it.
    - Otherwise, keep chrom as-is, but if prefix is a non-empty custom string, prepend it.
    """
    if not chrom:
        return chrom
    if prefix == 'chr':
        return chrom if chrom.lower().startswith('chr') else f"chr{chrom}"
    if prefix == '':
        return chrom[3:] if chrom.lower().startswith('chr') else chrom
    # custom prefix: do not double-add if already present
    return chrom if chrom.lower().startswith(prefix.lower()) else f"{prefix}{chrom}"

def create_renamed_reference(reference_file, output_file, mapping, primary_seqs, strategy, prefix=""):
    processed_chroms = defaultdict(int)
    with open(reference_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_seq_id = None
        include_current_seq = False
        for line in infile:
            if line.startswith('>'):
                header = line.strip()
                seq_id = header[1:].split()[0]
                current_seq_id = seq_id
                if seq_id in mapping:
                    chrom = mapping[seq_id]
                    out_chrom = _apply_prefix_strategy(chrom, prefix)
                    if strategy == 'primary':
                        if seq_id == primary_seqs.get(chrom):
                            new_header = f">{out_chrom}"
                            outfile.write(new_header + '\n')
                            include_current_seq = True
                        else:
                            include_current_seq = False
                    else:
                        processed_chroms[chrom] += 1
                        if processed_chroms[chrom] == 1:
                            new_header = f">{out_chrom}"
                        else:
                            new_header = f">{out_chrom}_{processed_chroms[chrom]}"
                        outfile.write(new_header + '\n')
                        include_current_seq = True
                else:
                    include_current_seq = False
            elif include_current_seq:
                outfile.write(line)

def main():
    args = parse_args()
    _log(f"Extracting chromosome mapping from {args.reference}...", args.quiet)
    _log(f"Mode: {args.mode}", args.quiet)
    if args.mode == 'template' and not args.template:
        print("Error: '-template' is required when using template mode", file=sys.stderr)
        sys.exit(1)
    use_clean = args.clean and not getattr(args, 'no_clean', False)
    mapping, seq_lengths = extract_chromosome_mapping(
        args.reference, args.mode, args.template, args.field, use_clean, quiet=args.quiet
    )
    if not mapping:
        print("Error: Failed to extract any chromosome mapping from FASTA headers.", file=sys.stderr)
        sys.exit(1)
    primary_seqs = find_primary_sequences(mapping, seq_lengths)
    chrom_counts = defaultdict(int)
    for chrom in mapping.values():
        chrom_counts[chrom] += 1
    duplicate_chroms = [chrom for chrom, count in chrom_counts.items() if count > 1]
    if duplicate_chroms and not args.quiet:
        _log(f"Duplicate chromosomes detected: {', '.join(duplicate_chroms[:5])}{'...' if len(duplicate_chroms) > 5 else ''}", args.quiet)
        _log(f"Using strategy '{args.strategy}' for duplicates", args.quiet)
    if not args.quiet:
        _log(f"Extracted {len(mapping)} seq IDs mapped to {len(set(mapping.values()))} unique chromosomes", args.quiet)
    if args.vcf_file and not args.quiet:
        _log(f"Validating chromosome names in VCF {args.vcf_file}...", args.quiet)
        vcf_chroms = validate_vcf_chromosomes(args.vcf_file, mapping)
        if vcf_chroms:
            display_count = min(10, len(vcf_chroms))
            sorted_vcf_chroms = sorted(vcf_chroms)
            _log(f"VCF chromosome examples: {', '.join(sorted_vcf_chroms[:display_count])}", args.quiet)
    # Determine prefix: if 'auto', infer from VCF; else use provided string
    prefix_to_apply = args.prefix
    if prefix_to_apply.lower() == 'auto':
        inferred = detect_vcf_prefix(args.vcf_file, quiet=args.quiet)
        if not args.quiet:
            _log(f"Auto-detected VCF prefix: '{inferred or ''}'", args.quiet)
        prefix_to_apply = inferred
    _log(f"Creating renamed reference file {args.output}...", args.quiet)
    create_renamed_reference(args.reference, args.output, mapping, primary_seqs, args.strategy, prefix=prefix_to_apply)
    _log(f"Done! Renamed reference saved to {args.output}", args.quiet)
    _log(f"Note: Only sequences with extracted chromosome names are included", args.quiet)

if __name__ == "__main__":
    main()
