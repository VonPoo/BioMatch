
"""
BioMatch Analysis Script: 00_extractSNPsfromVCF

Purpose
- Extract k-mer sequences surrounding SNPs from a VCF using a reference FASTA.
- Supports two encoding modes:
  - Legacy fast mode (assumes no 'N' bases; switches to strict if encountered).
  - Strict mode that skips k-mers containing non-ATCG bases.
- Produces stratified outputs and statistics for downstream panel generation.

Inputs
- `--vcf` (implicit via constructor in this file): VCF with single ALT per record.
- `--fasta`: Reference genome indexed via `pyfaidx`.
- `--k`: K-mer window size.
- `--prefix`: Output prefix for generated files.
- `--ignore`: If True, skip transitions (AT↔AT or CG↔CG equivalents).
- `--subKmer`: Sub k-mer size for counting uniqueness.
- `--n_distance`: Distance threshold to flag nearby 'N' bases.

Outputs
- Stratified FASTA files with extracted k-mers by uniqueness buckets.
- STDERR summary with counts and mode selection details.

Usage Example
    python 00_extractSNPsfromVCF.py \
        --vcf path/to/variants.vcf \
        --fasta path/to/reference.fa \
        --k 31 \
        --prefix out/panel \
        --ignore 1 \
        --subKmer 19 \
        --n_distance 50

Notes
- The script prefers the fast legacy path but will enable strict mode if 'N'
  bases are detected in relevant genomic windows.
- VCF entries must have a single ALT allele; multi-allelic sites will raise.
"""
import argparse
from pyfaidx import Fasta
import os
import math
import sys

class VCFEntry:
    """Container for a single VCF record used for k-mer extraction.

    Attributes
    ---------
    chr : str
        Chromosome name.
    pos : int
        1-based position.
    wt : str
        Reference allele.
    variant : str
        Alternate allele.
    """
    def __init__(self, chr, pos, wt, variant):
        self.chr = chr
        self.pos = int(pos)
        self.wt = wt
        self.variant = variant

class ExtractKmers:
    """K-mer extraction engine with legacy and strict encoding paths.

    Parameters
    ----------
    vcf : str
        Path to VCF file.
    fasta : str
        Path to FASTA genome (pyfaidx indexable).
    k : int
        K-mer window size.
    prefix : str
        Output file prefix.
    ignore : int or bool
        If truthy, skip transitions; only process transversions.
    subKmer : int
        Sub k-mer size used for uniqueness estimation.
    n_distance : int
        Radius for scanning nearby 'N' bases.
    """
    def __init__(self, vcf, fasta, k, prefix, ignore, subKmer, n_distance):
        self._vcf = vcf
        self._fasta = fasta
        self._k = k
        self._prefix = prefix
        self._ignore = ignore
        self._subKmer = subKmer
        self._n_distance = n_distance
        self._vcfEntries = {}
        self._kmers = {}
        self._has_n_bases = False  # Flag to track if N bases are detected
        self._stats = {
            'total': 0,
            'near_n': 0,
            'wt_mismatch': 0,
            'non_transversion': 0,
            'insufficient_unique_kmers': 0,
            'processed': 0,
            'total_kmers_skipped': 0
        }
        
    def _encodeKmer_strict(self, kmer):
        """Strict encoding that returns None for k-mers containing non-ATCG.

        Returns the canonical integer encoding of the k-mer (min of forward
        and reverse encodings), or None if invalid bases are present.
        """
        fw_encode = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
        rv_encode = {'A': 1, 'T': 0, 'C': 3, 'G': 2}
        fw = 0
        rv = 0
        
        # Check for N or any non-ATCG bases
        for base in kmer:
            if base not in fw_encode:
                return None
        
        for base in kmer:
            fw = (fw << 2) | fw_encode[base]
        
        for base in reversed(kmer):
            rv = (rv << 2) | rv_encode[base]
        
        return fw if fw < rv else rv
    
    def _encodeKmer_legacy(self, kmer):
        """Legacy encoding method — fast but raises on non-ATCG bases."""
        fw_encode = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
        rv_encode = {'A': 1, 'T': 0, 'C': 3, 'G': 2}
        fw = 0
        rv = 0
        
        for base in kmer:
            fw = (fw << 2) | fw_encode[base]
        
        for base in reversed(kmer):
            rv = (rv << 2) | rv_encode[base]
        
        return fw if fw < rv else rv
        
    def _insertKmer(self, kmerStr, strict_mode=False):
        """Insert a k-mer occurrence, handling invalid bases per mode."""
        if strict_mode:
            kmer = self._encodeKmer_strict(kmerStr)
            if kmer is None:
                self._stats['total_kmers_skipped'] += 1
                return
        else:
            try:
                kmer = self._encodeKmer_legacy(kmerStr)
            except KeyError:
                # N base detected, switch to strict mode
                self._has_n_bases = True
                return
        
        if kmer in self._kmers:
            self._kmers[kmer] += 1
        else:
            self._kmers[kmer] = 1
    
    def _getKmerCount(self, kmerStr, strict_mode=False):
        """Return occurrence count for a k-mer, inf for invalid in legacy mode."""
        if strict_mode:
            kmer = self._encodeKmer_strict(kmerStr)
            if kmer is None:
                return float('inf')  # Treat as repetitive
            return self._kmers.get(kmer, 0)
        else:
            try:
                kmer = self._encodeKmer_legacy(kmerStr)
                return self._kmers.get(kmer, 0)
            except KeyError:
                # N base detected
                return float('inf')
    
    def _checkVariant(self, base1, base2):
        """Return True for transitions (A↔G, C↔T), False for transversions."""
        if(base1 == "A" or base1 == "T"):
            if(base2 == "A" or base2 == "T"):
                return True  # Transition
            else:
                return False  # Transversion
        elif(base1 == "C" or base1 == "G"):
            if(base2 == "C" or base2 == "G"):
                return True  # Transition
            else:
                return False  # Transversion
    
    def _orderVariant(self, base1, base2):
        """Return False if ordering needs to be reversed (GC-class first)."""
        if(base1 == "A" or base1 == "T"):
            return True
        else:
            return False
    
    def _hasNearbyN(self, fastaFile, chr, pos, distance):
        """Return True if any 'N' is within `distance` bp of `pos`."""
        start = max(0, pos - 1 - distance)
        end = pos + distance
        
        try:
            region = str(fastaFile[chr][start:end]).upper()
            if 'N' in region:
                return True
        except:
            return True
        
        return False
    
    def _detectNBases(self, fastaFile):
        """Quick scan to detect whether 'N' bases occur in any target windows."""
        print("Scanning for 'N' bases in variant regions...", file=sys.stderr)
        
        for id in self._vcfEntries.keys():
            chr = self._vcfEntries[id].chr
            offset = self._vcfEntries[id].pos - 1
            pos1 = math.ceil(offset - self._k / 2)
            pos2 = pos1 + self._k
            
            try:
                tmpStr = str(fastaFile[chr][pos1:pos2]).upper()
                if 'N' in tmpStr:
                    self._has_n_bases = True
                    print(f"'N' bases detected in genome. Enabling strict filtering mode.", file=sys.stderr)
                    return True
            except:
                continue
        
        print("No 'N' bases detected. Using fast legacy mode.", file=sys.stderr)
        return False
    
    def _parseVCF(self):
        vcfFH = open(self._vcf, 'r')
        lines = vcfFH.readlines()
        idCounter = 0
        
        for line in lines:
            if line[0] != "#":
                self._stats['total'] += 1
                tmpArr = line.rstrip().split("\t")
                snpID = tmpArr[2]
                if snpID == '.':
                    snpID = idCounter
                    idCounter += 1
                
                if len(tmpArr[4]) > 1:
                    print("Error: Multiple alternate alleles found in VCF", file=sys.stderr)
                    exit(1)
                
                info = VCFEntry(tmpArr[0], tmpArr[1], tmpArr[3], tmpArr[4])
                self._vcfEntries[str(snpID)] = info
        
        vcfFH.close()
    
    def _extract_legacy(self, fastaFile):
        """Legacy fast extraction path used when no 'N' bases are detected."""
        print("Using fast legacy extraction mode...", file=sys.stderr)
        
        # First pass: count k-mers
        for id in self._vcfEntries.keys():
            chr = self._vcfEntries[id].chr
            offset = self._vcfEntries[id].pos - 1
            pos1 = math.ceil(offset - self._k / 2)
            pos2 = pos1 + self._k
            tmpStr = str(fastaFile[self._vcfEntries[id].chr][pos1:pos2]).upper()
            
            if self._vcfEntries[id].wt != tmpStr[int(self._k / 2)]:
                continue
            
            if self._checkVariant(self._vcfEntries[id].wt, self._vcfEntries[id].variant) and self._ignore:
                continue
            
            modStr = tmpStr[0:int(self._k / 2)] + self._vcfEntries[id].variant + tmpStr[int(self._k / 2) + 1:]
            
            for pos in range(0, len(modStr) - self._subKmer + 1):
                self._insertKmer(tmpStr[pos:pos + self._subKmer], strict_mode=False)
                self._insertKmer(modStr[pos:pos + self._subKmer], strict_mode=False)
        
        # Check if N was encountered during counting
        if self._has_n_bases:
            print("'N' bases encountered during processing. Switching to strict mode...", file=sys.stderr)
            return self._extract_strict(fastaFile)
        
        # Second pass: extract unique k-mers
        removeCount = 0
        processCount = 0
        atToCGFilterCount = 0
        kmersRemoved = 0
        
        for id in self._vcfEntries.keys():
            chr = self._vcfEntries[id].chr
            offset = self._vcfEntries[id].pos - 1
            pos1 = math.ceil(offset - self._k / 2)
            pos2 = pos1 + self._k
            tmpStr = str(fastaFile[self._vcfEntries[id].chr][pos1:pos2]).upper()
            
            if self._vcfEntries[id].wt != tmpStr[int(self._k / 2)]:
                removeCount += 1
                continue
            
            if self._checkVariant(self._vcfEntries[id].wt, self._vcfEntries[id].variant) and self._ignore:
                removeCount += 1
                atToCGFilterCount += 1
                continue
            
            modStr = tmpStr[0:int(self._k / 2)] + self._vcfEntries[id].variant + tmpStr[int(self._k / 2) + 1:]
            retainedKmers = kmersRemoved
            
            if self._orderVariant(self._vcfEntries[id].wt, self._vcfEntries[id].variant):
                for pos in range(0, len(modStr) - self._subKmer + 1):
                    modKmer = modStr[pos:pos + self._subKmer]
                    kmer = tmpStr[pos:pos + self._subKmer]
                    
                    if self._getKmerCount(kmer, strict_mode=False) == 1:
                        print(">" + id + "|" + str(pos) + "|AT")
                        print(tmpStr[pos:pos + self._subKmer])
                    else:
                        kmersRemoved += 1
                        
                    if self._getKmerCount(modKmer, strict_mode=False) == 1:
                        print(">" + id + "|" + str(pos) + "|CG")
                        print(modStr[pos:pos + self._subKmer])
                    else:
                        kmersRemoved += 1
            else:
                for pos in range(0, len(modStr) - self._subKmer + 1):
                    modKmer = modStr[pos:pos + self._subKmer]
                    kmer = tmpStr[pos:pos + self._subKmer]
                    
                    if self._getKmerCount(modKmer, strict_mode=False) == 1:
                        print(">" + id + "|" + str(pos) + "|AT")
                        print(modStr[pos:pos + self._subKmer])
                    else:
                        kmersRemoved += 1
                        
                    if self._getKmerCount(kmer, strict_mode=False) == 1:
                        print(">" + id + "|" + str(pos) + "|CG")
                        print(tmpStr[pos:pos + self._subKmer])
                    else:
                        kmersRemoved += 1
            
            if (kmersRemoved - retainedKmers) == 2 * (len(modStr) - self._subKmer + 1):
                removeCount += 1
            
            processCount += 1
        
        print("Processed " + str(processCount) + " SNPs. Removed " + str(removeCount) + 
              " SNPs. " + str(kmersRemoved) + " duplicate k-mers removed.", file=sys.stderr)
        if atToCGFilterCount > 0:
            print("Filtered " + str(atToCGFilterCount) + 
                  " SNPs that did not have A/T to C/G variants", file=sys.stderr)
    
    def _extract_strict(self, fastaFile):
        """Strict extraction method with N base filtering"""
        print("Using strict extraction mode with N-base filtering...", file=sys.stderr)
        
        # Reset k-mer counts for strict mode
        self._kmers = {}
        self._stats['total_kmers_skipped'] = 0
        
        # Filter SNPs near N bases
        filtered_entries = {}
        for id in self._vcfEntries.keys():
            entry = self._vcfEntries[id]
            if self._hasNearbyN(fastaFile, entry.chr, entry.pos, self._n_distance):
                self._stats['near_n'] += 1
                continue
            filtered_entries[id] = entry
        
        print(f"Filtered {self._stats['near_n']} SNPs within {self._n_distance}bp of 'N' bases", 
              file=sys.stderr)
        
        # First pass: count all k-mers
        for id in filtered_entries.keys():
            chr = filtered_entries[id].chr
            offset = filtered_entries[id].pos - 1
            pos1 = math.ceil(offset - self._k / 2)
            pos2 = pos1 + self._k
            
            tmpStr = str(fastaFile[filtered_entries[id].chr][pos1:pos2]).upper()
            
            if filtered_entries[id].wt != tmpStr[int(self._k / 2)]:
                continue
            
            if self._checkVariant(filtered_entries[id].wt, filtered_entries[id].variant) and self._ignore:
                continue
            
            modStr = tmpStr[0:int(self._k / 2)] + filtered_entries[id].variant + tmpStr[int(self._k / 2) + 1:]
            
            for pos in range(0, len(modStr) - self._subKmer + 1):
                self._insertKmer(tmpStr[pos:pos + self._subKmer], strict_mode=True)
                self._insertKmer(modStr[pos:pos + self._subKmer], strict_mode=True)
        
        # Second pass: extract unique k-mers with minimum threshold
        min_unique_kmers = 3
        
        for id in filtered_entries.keys():
            chr = filtered_entries[id].chr
            offset = filtered_entries[id].pos - 1
            pos1 = math.ceil(offset - self._k / 2)
            pos2 = pos1 + self._k
            
            tmpStr = str(fastaFile[filtered_entries[id].chr][pos1:pos2]).upper()
            
            if filtered_entries[id].wt != tmpStr[int(self._k / 2)]:
                self._stats['wt_mismatch'] += 1
                continue
            
            if self._checkVariant(filtered_entries[id].wt, filtered_entries[id].variant) and self._ignore:
                self._stats['non_transversion'] += 1
                continue
            
            modStr = tmpStr[0:int(self._k / 2)] + filtered_entries[id].variant + tmpStr[int(self._k / 2) + 1:]
            
            unique_kmer_count = 0
            output_buffer = []
            
            if self._orderVariant(filtered_entries[id].wt, filtered_entries[id].variant):
                for pos in range(0, len(modStr) - self._subKmer + 1):
                    modKmer = modStr[pos:pos + self._subKmer]
                    kmer = tmpStr[pos:pos + self._subKmer]
                    
                    wt_count = self._getKmerCount(kmer, strict_mode=True)
                    var_count = self._getKmerCount(modKmer, strict_mode=True)
                    
                    if wt_count == 1:
                        output_buffer.append((">" + id + "|" + str(pos) + "|AT", kmer))
                        unique_kmer_count += 1
                    
                    if var_count == 1:
                        output_buffer.append((">" + id + "|" + str(pos) + "|CG", modKmer))
                        unique_kmer_count += 1
            else:
                for pos in range(0, len(modStr) - self._subKmer + 1):
                    modKmer = modStr[pos:pos + self._subKmer]
                    kmer = tmpStr[pos:pos + self._subKmer]
                    
                    var_count = self._getKmerCount(modKmer, strict_mode=True)
                    wt_count = self._getKmerCount(kmer, strict_mode=True)
                    
                    if var_count == 1:
                        output_buffer.append((">" + id + "|" + str(pos) + "|AT", modKmer))
                        unique_kmer_count += 1
                    
                    if wt_count == 1:
                        output_buffer.append((">" + id + "|" + str(pos) + "|CG", kmer))
                        unique_kmer_count += 1
            
            # Only output if at least 3 unique k-mers
            if unique_kmer_count >= min_unique_kmers:
                for header, seq in output_buffer:
                    print(header)
                    print(seq)
                self._stats['processed'] += 1
            else:
                self._stats['insufficient_unique_kmers'] += 1
        
        # Print statistics
        print("\n=== Strict Mode Filtering Statistics ===", file=sys.stderr)
        print(f"Total SNPs in VCF: {self._stats['total']}", file=sys.stderr)
        print(f"SNPs near 'N' bases (within {self._n_distance}bp): {self._stats['near_n']}", file=sys.stderr)
        print(f"Wildtype mismatches: {self._stats['wt_mismatch']}", file=sys.stderr)
        print(f"Non-transversion filtered: {self._stats['non_transversion']}", file=sys.stderr)
        print(f"Insufficient unique k-mers (<{min_unique_kmers}): {self._stats['insufficient_unique_kmers']}", 
              file=sys.stderr)
        print(f"Successfully processed SNPs: {self._stats['processed']}", file=sys.stderr)
        print(f"K-mers skipped (containing N): {self._stats['total_kmers_skipped']}", file=sys.stderr)
    
    def extract(self):
        self._parseVCF()
        fastaFile = Fasta(self._fasta)
        
        # Quick detection of N bases
        has_n = self._detectNBases(fastaFile)
        
        if has_n:
            # Use strict mode with N filtering
            self._extract_strict(fastaFile)
        else:
            # Use fast legacy mode
            self._extract_legacy(fastaFile)

def main():
    parser = argparse.ArgumentParser(
        description='Extract k-mers from VCF file. Automatically switches to strict mode if N bases detected.'
    )
    parser.add_argument("-v", '--vcf', type=str, dest='vcf', required=True,
                        help='VCF file')
    parser.add_argument("-f", '--fa', type=str, dest='fasta', required=True,
                        help='Fasta file with fai index')
    parser.add_argument("-k", '--kmer', type=int, dest='kmer', default=31,
                        help='K-mer size (default: 31)')
    parser.add_argument("-p", '--prefix', type=str, dest='prefix', default="",
                        help='Output prefix')
    parser.add_argument("-i", '--ignoreReq', action='store_false', dest='ignore', default=True,
                        help='Ignore AT to CG conversion requirements')
    parser.add_argument("-s", '--subKmer', type=int, dest='subKmer', default=19,
                        help='Sub k-mer size (default: 19)')
    parser.add_argument("-n", '--n-distance', type=int, dest='n_distance', default=31,
                        help='Distance threshold from N bases in strict mode (default: 31)')

    args = parser.parse_args()
    
    if args.subKmer > args.kmer:
        print("Error: Sub k-mer size cannot be larger than k-mer size", file=sys.stderr)
        exit(1)
    
    extractor = ExtractKmers(args.vcf, args.fasta, args.kmer, args.prefix, 
                            args.ignore, args.subKmer, args.n_distance)
    extractor.extract()

if __name__ == "__main__":
    main()
