#!/usr/bin/env python3

import argparse
import csv
import random
import time

from collections import Counter
from matplotlib import pyplot as plt
from matplotlib_venn import venn2

def parse_args(args):
    parser = argparse.ArgumentParser(description="Takes two csv file and compares the three sets.")

    parser.add_argument("input1", help="FASTA input file", metavar="INPUT1.fasta")
    parser.add_argument("A_label", help="Label A")
    parser.add_argument("input2", help="FASTA input file", metavar="INPUT2.fasta")
    parser.add_argument("B_label", help="Label B")
    parser.add_argument("output", help="csv output file header, formatted as 'OUTPUT-#.csv'", metavar="OUTPUT")

    return parser.parse_args(args)

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

def extract_rows(file, kmer_list):
    start = time.time()
    entries = 0
    header_row = []
    rows = []
    print("Extracting rows from " + str(file))
    with open(file) as f:
        csv_reader = csv.reader(f, delimiter=',')
        header_row = next(csv_reader)
        for row in csv_reader:
            if str(row[0]) in kmer_list:
                entries += 1
                rows.append(row)
            else:
                pass
    f.close()
    interval = (time.time() - start)
    print("Processed " + str(entries) + " entries in " + str(format(interval, '.2f')) + " sec")
    print("Procesing rate: " + str(format(entries/interval, '.2f')) + " entries/sec\n")
    return (header_row, rows, entries)

def write_extracts(input, header_row, output, name_append):
    start = time.time()
    with open((str(output) + "-" + str(name_append) + ".csv"), 'w', newline='') as f:
        fwriter = csv.writer(f)
        fwriter.writerow(header_row)
        for item in input:
            fwriter.writerow(item)
    f.close()

def main(args):
    args = parse_args(args)

    kmers1 = parse_csv(args.input1)
    kmers2 = parse_csv(args.input2)
#    print(kmers1[0])

    A = set(kmers1[0])
    B = set(kmers2[0])

    venn2([A, B], set_labels = (args.A_label, args.B_label))
    plt.savefig(str(args.output + '.png'))

    AB = A.intersection(B)

    kmer1_extract = extract_rows(args.input1, AB)
    kmer2_extract = extract_rows(args.input2, AB)

    kmer1_extract_out = write_extracts(kmer1_extract[1], kmer1_extract[0], args.output, str(args.A_label + "fromABint"))
    kmer2_extract_out = write_extracts(kmer2_extract[1], kmer2_extract[0], args.output, str(args.B_label + "fromABint"))

    sys.exit(0)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
