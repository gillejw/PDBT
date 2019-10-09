#!/usr/bin/env python3

import argparse
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import time

from itertools import product
from matplotlib import pyplot as plt

def parse_args(args):
    parser = argparse.ArgumentParser(description="Takes a collection of csv files containing motifs and returns a heatmap.")

    parser.add_argument("input", help="csv input file", metavar="INPUT.csv")
#    parser.add_argument("output", help="csv output file header, formatted as 'OUTPUT-#.csv'", metavar="OUTPUT")

    return parser.parse_args(args)

def parse_csv(csv_file, label):
    df = pd.read_csv(csv_file, header=None)
    extracted_col = df[df.columns[2]]
    return (list(extracted_col), label)

def parse_input(list):
    files = []
    with open(list, 'r') as f:
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            tuple = (str(row[0]), str(row[1]))
            files.append(tuple)
    f.close()
    return files

def main(args):
    args = parse_args(args)

    infiles = parse_input(args.input)

#    brain_df = parse_csv(args.input, "Brain")
#    brain_df2 = parse_csv(args.input, "Brain2")

    AA_alpha = 'ACDEFGHIKLMNPQRSTVWY'
    c = []
    for seq in product(AA_alpha, repeat=3):
        seq_str = ''.join(seq)
        c.append(seq_str)

    df = pd.DataFrame(index=c)
    for file in infiles:
        s = parse_csv(file[0], file[1])
        df[s[1]] = s[0]
        df[s[1]] *= 100
#    df[brain_df[1]] = brain_df[0]
#    df[brain_df2[1]] = brain_df2[0]
    print(df)
    df_no0 = df.loc[(df!=0).any(1)]
    print(df_no0)
    df_less = df_no0.loc[(df>0.0125).any(1)]
    print(df_less)
#    df_transposed = df.transpose()
    heat_map = sns.clustermap(df, yticklabels=False, robust=True, cmap="YlOrRd")
#    heat_map = sns.heatmap(df, yticklabels=False, robust=True,cmap="cubehelix",)
    plt.savefig(str('heatmap.png'))

    heat_map2 = sns.clustermap(df_no0, yticklabels=False, robust=True, cmap="YlOrRd", method='ward')
    plt.savefig(str('heatmap2.png'))

    heat_map3 = sns.clustermap(df_less, yticklabels=False, robust=True, cmap="YlOrRd", metric="euclidean", method="ward")
    heat_map3.ax_row_dendrogram.set_visible(False)
    plt.savefig(str('heatmap3.png'))

    sys.exit(0)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
