#!/usr/bin/env python3

import argparse
import csv
import time
import scipy.stats as stats
import numpy as np

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from collections import Counter, OrderedDict
from itertools import product
from math import log2, log10
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from statsmodels.stats.multitest import fdrcorrection

def parse_args(args):
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("input", help="FASTA input file", metavar="INPUT.fasta")
    parser.add_argument("output", help="output file header", metavar="OUTPUT")
    parser.add_argument("unselected", help="FASTA file with non-target peptides", metavar="UNSELECTED.fasta")
    parser.add_argument("-k", "--kmer", help="Length of k-mers to analyze", metavar="size", default=int(0))

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

def nuc_to_pep(nuc_seq):
    start = time.time()
    print("Calculating unique peptide sequences...")
    entries = 0
    bad_reads = 0
    pep_seq = Counter()
#    skip = ['DMP','MPG','PGT','GTV','TVL','VLP','LPD']
    for i in nuc_seq:
        entries += 1
        if len(str(i)) < 27:
            bad_reads += 1
        else:
            s = str(Seq(str(i), IUPAC.unambiguous_dna).translate())
#            if "*" in s or any(ele in s for ele in skip):
            if "*" in s or "DMPGTVLPD" in s:
                bad_reads += 1
#                print(s)
            else:
                pep_seq[s] += 1
    interval = (time.time() - start)
    print("Processed " + str(entries) + " entries in " + str(format(interval, '.2f')) + " sec")
    print("Processing rate: " + str(format(entries/interval, '.2f')) + " entries/sec")
    print("Number of peptide sequences with bad translations: " + str(bad_reads) + "\n")
    return (pep_seq, entries)

def k_mer_count(sample, k_mer_size):
    sequences = 0
    k_mers = 0
    k_mer_dict = Counter()
    for line in sample:
        sequences += 1
        for i in range(len(line) - k_mer_size + 1):
            k_mer = line[i:i+k_mer_size]
            k_mer_dict[k_mer] += 1
            k_mers += 1
#    print(k_mer_dict.keys())
    return (k_mer_dict, k_mers)

def main(args):
    args = parse_args(args)

    seq = nuc_to_pep(parse_fasta(args.input)[0])
    seq_list = list(seq[0])
    unselected = parse_fasta(args.unselected)
    unselected_list = list(unselected[0])
    total_list = seq_list + unselected_list
#    print(seq_list)
#    print(unselected_list)

    k_mers = k_mer_count(seq_list, int(args.kmer))
    unselected_k_mers = k_mer_count(unselected_list, int(args.kmer))
    total_kmers = k_mer_count(total_list, int(args.kmer))
#    print(k_mers[0])
#    print(unselected_k_mers)
#    print(k_mers[1])

    k_mer_pval = {}
    significant_k_mers = {}
    significant_k_mers2 = {}
    occurance_k_mers = {}
    occurance_k_mers2 = {}
    fold_k_mers = {}
    key_count_k_mers = {}
    not_key_count_k_mers = {}
    key_unselected_count_k_mers = {}
    not_unselected_key_count_k_mers = {}
    skip = ['DMP','MPG','PGT','GTV','TVL','VLP','LPD']
    skipped_rows = 0
    for key in k_mers[0]:
        if str(key) in skip:
            skipped_rows += 1
        else:
            key_count = k_mers[0][key]
            not_key_count = k_mers[1]-key_count
            key_unselected_count = unselected_k_mers[0][key]
            not_unselected_key_count = unselected_k_mers[1]-key_unselected_count
#           print(key)
#           print(key_count)
#           print(not_key_count)
#           print(key_unselected_count)
#           print(not_unselected_key_count)
            pvalue = stats.fisher_exact([[key_count, not_key_count], [key_unselected_count, not_unselected_key_count]], alternative='greater')
            pvalue2 = stats.fisher_exact([[key_count, not_key_count], [key_unselected_count, not_unselected_key_count]], alternative='two-sided')
            k_mer_pval[key] = pvalue[1]
#           print(pvalue[1])
            if pvalue[1] <= 0.05:
                significant_k_mers[key] = pvalue[1]
                occurance_k_mers[key] = (((key_count)/(not_key_count + key_count))*100)
                occurance_k_mers2[key] = (((key_count)/(key_count + key_unselected_count))*100)
                key_count_k_mers[key] = key_count
                not_key_count_k_mers[key] = not_key_count
                key_unselected_count_k_mers[key] = key_unselected_count
                not_unselected_key_count_k_mers[key] = not_unselected_key_count
                try:
                    fold_k_mers[key] = log2((key_count/(not_key_count+key_count))/(key_unselected_count/not_unselected_key_count))
                except:
                    fold_k_mers[key] = str("NA")

            if pvalue2[1] <= 0.05:
                significant_k_mers2[key] = pvalue2[1]
                occurance_k_mers[key] = (((key_count)/(not_key_count + key_count))*100)
                occurance_k_mers2[key] = (((key_count)/(key_count + key_unselected_count))*100)
                key_count_k_mers[key] = key_count
                not_key_count_k_mers[key] = not_key_count
                key_unselected_count_k_mers[key] = key_unselected_count
                not_unselected_key_count_k_mers[key] = not_unselected_key_count
                try:
                    fold_k_mers[key] = log2((key_count/(not_key_count+key_count))/(key_unselected_count/not_unselected_key_count))
                except:
                    fold_k_mers[key] = str("0")

#            if str(key) == str("HDQ") or str(key) == str("RLF"):
#                print(key)
#                print(key_count)
#                print(not_key_count)
#                print(key_unselected_count)
#                print(not_unselected_key_count)
#                print(pvalue[1])

#    print(k_mer_pval.values())
    k_mer_key = list(k_mer_pval.keys())
#    print(k_mer_key)
    adj_pval = fdrcorrection(np.array(list(k_mer_pval.values())), alpha=0.05, method='indep')
    k_mer_adj_pval = dict(zip(k_mer_key, adj_pval[1]))
#    print(k_mer_adj_pval)
    sorted_significant_k_mers = OrderedDict(sorted(significant_k_mers.items(), key=lambda x: x[1]))
    sorted_significant_k_mers2 = OrderedDict(sorted(significant_k_mers2.items(), key=lambda x: x[1]))

    AA_alpha = 'ACDEFGHIKLMNPQRSTVWY'
    trimer_AA = []
    for sequence in product(AA_alpha, repeat=int(args.kmer)):
        seq_str = ''.join(sequence)
        trimer_AA.append(seq_str)

    with open(args.output + "-" + args.kmer + "-mer-summary.csv", 'w', newline='') as f:
        fwriter = csv.writer(f)
        for item in trimer_AA:
            try:
                fwriter.writerow((item, k_mers[0][item], (k_mers[0][item]/k_mers[1])))
            except:
                fwriter.writerow((item, 0, 0))
    f.close()

#    print(k_mers[0])
    with open(args.output + "-" + args.kmer + "-mer-summary.txt", 'w') as f:
        f.write("Summary of K-mers")
        f.write("\n")
        for item in sorted_significant_k_mers:
            f.write(item + ": " + str(sorted_significant_k_mers[item]))
            f.write("\n")
    f.close()

    with open(args.output + "-" + args.kmer + "-mer-summary-single-tailed.csv", 'w', newline='') as f:
        fwriter = csv.writer(f)
        fwriter.writerow(("Motif", "Motif Count in Tissue", "Non-Motif Count in Tissue","Motif Count in Other", "Non-Motif Count in Other", "Occurance in Tissue (%)", "Occurance in Collection (%)", "Fold Change (log2)", "P-value", "P-value (-log10(p))", "Adjusted P-value", "Adj. P-value (-log10(p))"))
        for item in sorted_significant_k_mers:
            try:
                fwriter.writerow((item, str(key_count_k_mers[item]), str(not_key_count_k_mers[item]), str(key_unselected_count_k_mers[item]), str(not_unselected_key_count_k_mers[item]), str(occurance_k_mers[item]), str(occurance_k_mers2[item]), str(fold_k_mers[item]), str(sorted_significant_k_mers[item]), str(-log10(sorted_significant_k_mers[item])), str(k_mer_adj_pval[item]), str(-log10(k_mer_adj_pval[item]))))
            except:
                fwriter.writerow((item, str(key_count_k_mers[item]), str(not_key_count_k_mers[item]), str(key_unselected_count_k_mers[item]), str(not_unselected_key_count_k_mers[item]), str(occurance_k_mers[item]), str(occurance_k_mers2[item]), str(fold_k_mers[item]), str(sorted_significant_k_mers[item]), str(0)))
    f.close()

    with open(args.output + "-" + args.kmer + "-mer-adj-summary-single-tailed.csv", 'w', newline='') as f:
        fwriter = csv.writer(f)
        fwriter.writerow(("Motif", "Motif Count in Tissue", "Non-Motif Count in Tissue","Motif Count in Other", "Non-Motif Count in Other", "Occurance in Tissue (%)", "Occurance in Collection (%)", "Fold Change (log2)", "P-value", "P-value (-log10(p))", "Adjusted P-value", "Adj. P-value (-log10(p))"))
        for item in sorted_significant_k_mers:
            if k_mer_adj_pval[item] < 0.05:
                try:
                    fwriter.writerow((item, str(key_count_k_mers[item]), str(not_key_count_k_mers[item]), str(key_unselected_count_k_mers[item]), str(not_unselected_key_count_k_mers[item]), str(occurance_k_mers[item]), str(occurance_k_mers2[item]), str(fold_k_mers[item]), str(sorted_significant_k_mers[item]), str(-log10(sorted_significant_k_mers[item])), str(k_mer_adj_pval[item]), str(-log10(k_mer_adj_pval[item]))))
                except:
                    fwriter.writerow((item, str(key_count_k_mers[item]), str(not_key_count_k_mers[item]), str(key_unselected_count_k_mers[item]), str(not_unselected_key_count_k_mers[item]), str(occurance_k_mers[item]), str(occurance_k_mers2[item]), str(fold_k_mers[item]), str(sorted_significant_k_mers[item]), str(0)))
            else:
                pass
    f.close()

    with open(args.output + "-" + args.kmer + "-mer-summary-two-tailed.csv", 'w', newline='') as f:
        fwriter = csv.writer(f)
        fwriter.writerow(("Motif", "Motif Count in Tissue", "Non-Motif Count in Tissue","Motif Count in Other", "Non-Motif Count in Other", "Occurance in Tissue (%)", "Occurance in Collection (%)", "Fold Change (log2)", "P-value", "P-value (-log10(p))", "Adjusted P-value"))
        for item in sorted_significant_k_mers2:
            try:
                fwriter.writerow((item, str(key_count_k_mers[item]), str(not_key_count_k_mers[item]), str(key_unselected_count_k_mers[item]), str(not_unselected_key_count_k_mers[item]), str(occurance_k_mers[item]), str(occurance_k_mers2[item]), str(fold_k_mers[item]), str(sorted_significant_k_mers2[item]), str(-log10(sorted_significant_k_mers2[item])), str(k_mer_adj_pval[item])))
            except:
                fwriter.writerow((item, str(key_count_k_mers[item]), str(not_key_count_k_mers[item]), str(key_unselected_count_k_mers[item]), str(not_unselected_key_count_k_mers[item]), str(occurance_k_mers[item]), str(occurance_k_mers2[item]), str(fold_k_mers[item]), str(sorted_significant_k_mers2[item]), str(0)))

#    print(sorted_significant_k_mers.items())
#    print(len(sorted_significant_k_mers))
    print("Processed " + str(len(k_mers[0])) + " unique k-mers in the sample from " + str(len(total_kmers[0])) + " total k-mers")
    print("Skipped " + str(skipped_rows) + " k-mers for invalid keys")
#    for key in k_mers:
#        print(key)
#        print(k_mers[key])

    A = set(list(k_mers[0]))
    B = set(list(unselected_k_mers[0]))
    venn2([A, B], set_labels = ("A", "B"))
    plt.savefig(str(args.output + "-" + args.kmer + '-mer.png'))

    sys.exit(0)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
