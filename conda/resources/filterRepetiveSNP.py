#!/usr/bin/env python3
import sys

def load_snps(path):
    snps = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            chrom, pos, ref, alt = parts[:4]
            snps.append((chrom, pos, ref, alt))
    return snps

def filter_repetitive(snps):
    seen = set()
    unique = []
    for s in snps:
        key = (s[0], s[1])
        if key in seen:
            continue
        seen.add(key)
        unique.append(s)
    return unique

def write_snps(snps, path):
    with open(path, 'w') as f:
        for s in snps:
            f.write('\t'.join(map(str, s)) + '\n')

def main():
    if len(sys.argv) < 3:
        print("Usage: filterRepetiveSNP.py <input.tsv> <output.tsv>")
        sys.exit(1)
    inp = sys.argv[1]
    outp = sys.argv[2]
    snps = load_snps(inp)
    uniq = filter_repetitive(snps)
    write_snps(uniq, outp)

if __name__ == '__main__':
    main()