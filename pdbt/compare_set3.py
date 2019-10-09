#!/usr/bin/env python3
# build 05082019_XXX

import argparse
import csv
import random
import time

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from collections import Counter
from matplotlib import pyplot as plt
from matplotlib_venn import venn3

def parse_args(args):
    parser = argparse.ArgumentParser(description="Takes three FASTA formatted file and compares the three sets.")

    parser.add_argument("input1", help="FASTA input file", metavar="INPUT1.fasta")
    parser.add_argument("A_label", help="Label A")
    parser.add_argument("input2", help="FASTA input file", metavar="INPUT2.fasta")
    parser.add_argument("B_label", help="Label B")
    parser.add_argument("input3", help="FASTA input file", metavar="INPUT3.fasta")
    parser.add_argument("C_label", help="Label C")
    parser.add_argument("output", help="csv output file header, formatted as 'OUTPUT-#.csv'", metavar="OUTPUT")

    return parser.parse_args(args)

def file_type_checker(file):
    file = file.lower()
    if file[-4:] == '.txt':
        return ("txt")
    if file[-6:] == '.fasta':
        return ("fasta")
    if file[-6:] == '.fastq':
        return ('fastq')

def parse_fasta(file):
    start = time.time()
    print("Parsing FASTA file using SeqIO.parse...")
    sequences = Counter()
    type = file_type_checker(file)
    entries = 0
    for record in SeqIO.parse(file, type):
        entries += 1
        sequences[str(record.seq)] += 1
    interval = (time.time() - start)
    print("Processed " + str(entries) + " entries in " + str(format(interval, '.2f')) + " sec")
    print("Processing rate: " + str(format(entries/interval, '.2f')) + " entries/sec\n")
    return (sequences, entries)

def nuc_to_pep(nuc_seq):
    start = time.time()
    print("Calculating unique peptide sequences...")
    entries = 0
    bad_reads = 0
    pep_seq = Counter()
    for i in nuc_seq:
        entries += 1
        if len(str(i)) < 27:
            bad_reads += 1
        else:
            s = str(Seq(str(i), IUPAC.unambiguous_dna).translate())
            if "*" in s:
                bad_reads += 1
            else:
                pep_seq[s] += 1
    interval = (time.time() - start)
    print("Processed " + str(entries) + " entries in " + str(format(interval, '.2f')) + " sec")
    print("Processing rate: " + str(format(entries/interval, '.2f')) + " entries/sec")
    print("Number of peptide sequences with bad translations: " + str(bad_reads) + "\n")
    return (pep_seq, entries)

def main(args):
    args = parse_args(args)
    sort_type = ""
    #print(file_type_checker(args.input))

    seq1 = nuc_to_pep(parse_fasta(args.input1)[0])
    seq2 = nuc_to_pep(parse_fasta(args.input2)[0])
    seq3 = nuc_to_pep(parse_fasta(args.input3)[0])
#    print(seq[0])
    list_pep1 = list(seq1[0])
    list_pep2 = list(seq2[0])
    list_pep3 = list(seq3[0])
    A = set(list_pep1)
    B = set(list_pep2)
    C = set(list_pep3)
    venn3([A, B, C], set_labels = (args.A_label, args.B_label, args.C_label))
    plt.savefig(str(args.output + '.png'))

    #ABC = A.intersection(B,C)
    #print(ABC)
    print(args.B_label + ":" + args.C_label + " Intersection...")
    BC = B.intersection(C)
    for i in sorted(list(BC)):
        print(i + ": " + args.B_label + " (" + str(seq2[0][i]) + ") | " + args.C_label + " (" + str(seq3[0][i]) + ")")
#        print(seq2[0][i])
#        print(seq3[0][i])
#    print(sorted(list(BC)))
    print("\n")

    print(args.A_label + ":" + args.C_label + " Intersection...")
    AC = A.intersection(C)
    for i in sorted(list(AC)):
        print(i + ": " + args.A_label + " (" + str(seq1[0][i]) + ") | " + args.C_label + " (" + str(seq3[0][i]) + ")")
#    print(sorted(list(AC)))
    print("\n")

    print(args.A_label + ":" + args.B_label + " Intersection...")
    AB = A.intersection(B)
    for i in sorted(list(AB)):
        print(i + ": " + args.A_label + " (" + str(seq1[0][i]) + ") | " + args.B_label + " (" + str(seq2[0][i]) + ")")
#    print(sorted(list(AB)))
    print("\n")
#    NotAB = A.symmetric_difference(B)
#    NotABlist = list(NotAB)

    #print(list_pep1)

    sys.exit(0)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
