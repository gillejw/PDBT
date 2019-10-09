#!/usr/bin/env python3
# build 05082019_1414

import argparse
import csv
import random
import time

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from collections import Counter

def parse_args(args):
    parser = argparse.ArgumentParser(description="Takes a FASTA formatted file and outputs a list of unique peptide sequences in csv format.")

    parser.add_argument("input", help="FASTA input file", metavar="INPUT.fasta")
    parser.add_argument("output", help="csv output file header, formatted as 'OUTPUT-#.csv'", metavar="OUTPUT")
    parser.add_argument("-a", "--abc", help="Sort peptide list alphabetically", action="store_true")
    parser.add_argument("-n", "--number", metavar="num", help="Number of peptides per output file", default=int(0))
    parser.add_argument("-r", "--randomize", help="Shuffle peptide list", action="store_true")
    parser.add_argument("--rank", help="Sort peptide list by rank", action="store_true")

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

def chunk_file(input, output, chunk_size, s_type):
    start = time.time()
    print("Chunking peptides " + str(s_type) + " into files containing " + str(chunk_size) + " sequences each...")
    files = 0
    if int(chunk_size) > 0:
        for i in range(0, int(len(input)), int(chunk_size)):
            files += 1
    #        print(i)
    #        print(input[i:i+int(chunk_size)])
            with open((output + "-" + str(files) + ".csv"), 'w', newline='') as f:
                fwriter = csv.writer(f)
                for k in input[i:i+int(chunk_size)]:
                    fwriter.writerow([k])
            f.close()
    else:
        files = 1
        with open((output + ".csv"), 'w', newline='') as f:
            fwriter = csv.writer(f)
            for k in input:
                fwriter.writerow([k])
        f.close()

    interval = (time.time() - start)
    print("Processed " + str(len(input)) + " sequences into " + str(files) + " different files in " + str(format(interval, '.2f')) + " sec")
    print("Procesing rate: " + str(format(files/interval, '.2f')) + " files/sec")
    return

def main(args):
    args = parse_args(args)
    sort_type = ""
    #print(file_type_checker(args.input))

#    nuc_seq = parse_fasta(args.input)
#    pep_seq = nuc_to_pep(seq[0])
    seq = nuc_to_pep(parse_fasta(args.input)[0])
#    print(seq[0])
    list_pep = list(seq[0])
#    print("Original List:")
#    print(list_pep)
    if args.abc == True and args.rank == False and args.randomize == False:
        list_pep.sort()
        sort_type = "alphabetically"
    elif args.randomize == True and args.abc == False and args.rank == False:
        random.shuffle(list_pep)
        sort_type = "randomly"
    elif args.rank == True and args.abc == False and args.randomize == False:
        list_pep = [pair[0] for pair in sorted(seq[0].items(), key=lambda item: item[1])]
        sort_type = "ranked"
    elif args.rank == False and args.abc == False and args.randomize == False:
        pass
    else:
        print("Please select a single sorting option.")
        sys.exit(0)

#    print("Sorted List:")
#    print(list_pep)

    list_pep_chunk = chunk_file(list_pep, args.output, args.number, sort_type)

    sys.exit(0)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
