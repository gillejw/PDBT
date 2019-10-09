#!/usr/bin/env python3

import argparse
import csv
import time

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from collections import Counter, OrderedDict

def parse_args(args):
    parser = argparse.ArgumentParser(description="Takes a FASTA formatted file with peptide sequences and a CSV file containing motifs and outputs a list of unique peptide sequences with each k-mer in csv format.")

    parser.add_argument("kmer_input", help="CSV input file containing kmer sequences", metavar="INPUT.csv")
    parser.add_argument("pep_input", help="FASTA input file containing peptide sequences", metavar="INPUT.fasta")
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

def parse_csv(file):
    start = time.time()
    entries = 0
    sequences = []
    skip = ['DMP','MPG','PGT','GTV','TVL','VLP','LPD']
    skipped_rows = 0
    print("Parsing CSV file " + str(file))
    with open(file) as f:
        csv_reader = csv.reader(f, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            if str(row[0]) in skip:
                skipped_rows += 1
            else:
                entries += 1
#               print(str(row[0]))
                sequences.append(str(row[0]))
    f.close()
    interval = (time.time() - start)
    print("Processed " + str(entries) + " entries in " + str(format(interval, '.2f')) + " sec")
    print("Procesing rate: " + str(format(entries/interval, '.2f')) + " entries/sec")
    print("Skipped: " + str(skipped_rows) + " rows\n")
    return (sequences, entries)

def parse_fasta(file):
    start = time.time()
    print("Parsing FASTA file using SeqIO.parse...")
    sequences = Counter()
    type = file_type_checker(file)
    entries = 0
    bad_reads = 0
    for record in SeqIO.parse(file, type):
        entries += 1
        if "*" in str(record.seq):
            bad_reads += 1
        else:
            sequences[str(record.seq)] += 1
    interval = (time.time() - start)
    print("Processed " + str(entries) + " entries in " + str(format(interval, '.2f')) + " sec")
    print("Processing rate: " + str(format(entries/interval, '.2f')) + " entries/sec\n")
    return (sequences, entries)

def pep_with_kmer(kmers, peps):
    start = time.time()
    sequences = {}
    entries = 0
    print("Collecting peptide sequences containing k-mers...")
#    print(peps.keys())
    for k in kmers:
#        print(k)
        key_list = []
        for key in peps.keys():
            if str(k) in str(key):
                key_list.append(key)
                entries += 1
#                print(key)
            else:
                pass
        sequences[str(k)] = (key_list)
    interval = (time.time() - start)
    print("Processed " + str(entries) + " entries in " + str(format(interval, '.2f')) + " sec")
    return sequences

def main(args):
    args = parse_args(args)

    kmer_list = parse_csv(args.kmer_input)
    peptide_list = parse_fasta(args.pep_input)
    kmer_peps = pep_with_kmer(kmer_list[0], peptide_list[0])

    with open(args.output + "-sequence-summary.txt", 'w') as f:
        for k in kmer_peps:
            f.write(k + ": " + str(len(kmer_peps[k])))
            f.write("\n")
            for s in kmer_peps[k]:
                f.write(s)
                f.write("\n")
    f.close()

#    print(kmer_list[0])
#    print(type(peptide_list[0]))

    sys.exit(0)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
