#!/usr/bin/env python3

import argparse
import csv
import networkx as nx
import time

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from collections import Counter, defaultdict

def parse_args(args):
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("motif_input", help="csv input file of motifs", metavar="INPUT.csv")
    parser.add_argument("pep_input", help="FASTA input file of peptides", metavar="INPUT2.fasta")

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

def csv_node_parser(file):
    start = time.time()
    nodes = []
    print("Parsing CSV file " + str(file) + " for nodes...")
    with open(file) as f:
        csv_reader = csv.reader(f, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            nodes.append(str(row[0]))
    interval = (time.time() - start)
    print("Processed " + str(len(nodes)) + " nodes in " + str(format(interval, '.2f')) + " sec")
    return nodes

def edge_parser(peptides, nodes):
    edges = []
    d = defaultdict(Counter)
    for key in peptides.keys():
        for node1 in nodes:
            for node2 in nodes:
                if node1 in key and node2 in key and (key.index(node1) is not key.index(node2) and key.index(node1) < key.index(node2) or node2 in key[key.index(node1)+3:]):
                    if node1 < node2:
                        d[node2][node1] += 1
                    if node2 <= node1:
                        d[node1][node2] += 1
#                    print(node1 + " " + node2)
#                    print(key)
#                    print(key.index(node1))
#                    print(key.index(node2))
#                    d[node1][node2] += 1
#    print(d)
    for u in d:
        for v in d[u]:
#            print(u)
#            print(v)
#            print(d[u][v])
            edge = (u,v,d[u][v])
#            print(edge)
            edges.append(edge)
    return edges

def main(args):
    args = parse_args(args)

    peptide_list = parse_fasta(args.pep_input)
#    peptide_list = {'HDQETYCCC':1,'HDQETYMYW':1,'ETYHDQETY':1,'HDQHDQEHY':1}
    node_list = csv_node_parser(args.motif_input)
#    node_list = ['DMP','MPG','PGT','GTV','TVL','VLP','LPD']
    edge_list = edge_parser(peptide_list[0], node_list)
#    edge_list = [("HDQ","VPY"), ("HKA","ETY")]

#    print(edge_list)

    G = nx.Graph()
    G.add_nodes_from(node_list)
    G.add_weighted_edges_from(edge_list)

#    for (u,v,d) in G.edges(data='weight'):
#        print(d)

    with open("./test/test_graph.graphml", "wb") as f:
        nx.write_graphml(G, f)
        f.close()

    sys.exit(0)

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
