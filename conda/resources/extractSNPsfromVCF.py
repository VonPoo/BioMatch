#!/usr/bin/env python3
import argparse
from pyfaidx import Fasta
import os
import math
import sys

# Optimized extractor aligned with workspace 00_extractSNPsfromVCF.py
# Slimmed to work as a resource script for ntsm-scripts

class VCFEntry:
    def __init__(self, chr, pos, wt, variant):
        self.chr = chr
        self.pos = int(pos)
        self.wt = wt
        self.variant = variant

def is_transversion(base1, base2):
    if (base1 in "AT" and base2 in "CG") or (base1 in "CG" and base2 in "AT"):
        return True
    return False

def encode_kmer(kmer):
    mapping = {'A':0,'T':1,'C':2,'G':3}
    try:
        fw = 0
        rv = 0
        for b in kmer:
            fw = (fw<<2)|mapping[b]
        for b in reversed(kmer):
            rv = (rv<<2)|({ 'A':1,'T':0,'C':3,'G':2 }[b])
        return fw if fw < rv else rv
    except KeyError:
        return None

def main():
    p = argparse.ArgumentParser(description="Optimized SNP/k-mer extractor")
    p.add_argument("vcf")
    p.add_argument("fasta")
    p.add_argument("k", type=int)
    p.add_argument("prefix")
    p.add_argument("ignore_transversion", type=int, choices=[0,1])
    p.add_argument("subKmer", type=int)
    p.add_argument("n_distance", type=int)
    args = p.parse_args()

    fasta = Fasta(args.fasta)
    entries = {}
    total = 0

    with open(args.vcf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            total += 1
            tmp = line.rstrip().split('\t')
            snpID = tmp[2]
            if snpID == '.':
                snpID = total
            entries[str(snpID)] = VCFEntry(tmp[0], tmp[1], tmp[3], tmp[4])

    kmers = {}
    def insert_kmer(s):
        enc = encode_kmer(s)
        if enc is None:
            return
        kmers[enc] = kmers.get(enc,0)+1

    # count sub-kmers
    for id, e in entries.items():
        offset = e.pos-1
        pos1 = math.ceil(offset - args.k/2)
        pos2 = pos1 + args.k
        seq = str(fasta[e.chr][pos1:pos2]).upper()
        if e.wt != seq[int(args.k/2)]:
            continue
        if args.ignore_transversion and not is_transversion(e.wt, e.variant):
            continue
        mod = seq[:int(args.k/2)] + e.variant + seq[int(args.k/2)+1:]
        for pos in range(0, len(mod)-args.subKmer+1):
            insert_kmer(seq[pos:pos+args.subKmer])
            insert_kmer(mod[pos:pos+args.subKmer])

    # emit unique sub-kmers
    for id, e in entries.items():
        offset = e.pos-1
        pos1 = math.ceil(offset - args.k/2)
        pos2 = pos1 + args.k
        seq = str(fasta[e.chr][pos1:pos2]).upper()
        if e.wt != seq[int(args.k/2)]:
            continue
        if args.ignore_transversion and not is_transversion(e.wt, e.variant):
            continue
        mod = seq[:int(args.k/2)] + e.variant + seq[int(args.k/2)+1:]
        for pos in range(0, len(mod)-args.subKmer+1):
            ref_s = seq[pos:pos+args.subKmer]
            var_s = mod[pos:pos+args.subKmer]
            enc_r = encode_kmer(ref_s)
            enc_v = encode_kmer(var_s)
            if enc_r is None or enc_v is None:
                continue
            if kmers.get(enc_r,0) == 1:
                print(f">{id}|{pos}|AT")
                print(ref_s)
            if kmers.get(enc_v,0) == 1:
                print(f">{id}|{pos}|CG")
                print(var_s)

if __name__ == "__main__":
    main()