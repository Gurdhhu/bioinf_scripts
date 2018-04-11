#!/usr/bin/env python3
from Bio import SeqIO
import re
from argparse import ArgumentParser, RawDescriptionHelpFormatter

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description='''
This script takes TSV table produced by VSEARCH after search of OTU representative
sequences across SILVA database (--usearch_global command) and for each database match 
adds 7-level taxonomic annotation, taking taxonomy lines from sequence names in
FASTA-formatted SILVA database. Writes output TSV table adding ".tax.tsv" instead of ".tsv"
to the original input filename.
Example usage: ./add_silva_taxonomy.py all.otus.991.silva.tsv ./DBs/SILVA_132_upper_T.fas
                                    ''')
parser.add_argument("tsv",
                    type=str,
                    help="Pathway to tsv-formatted hit table file with --usearch_global results")
parser.add_argument("fasta",
                    type=str,
                    help="FASTA file with taxonomy in sequence names")

args = parser.parse_args()
tsvfile = args.tsv
fastafile = args.fasta

seqdict = {}
with open(tsvfile, "r") as f:
    lines = f.readlines()

for line in lines:
    OTU, refID, sim = line.strip().split()
    if refID not in seqdict:
        seqdict[refID] = [[OTU, sim]]
    else:
        seqdict[refID].append([OTU, sim])

record_dict = SeqIO.index(fastafile, "fasta")
for key in record_dict:
    if key in seqdict:
        seqdict[key].append(record_dict[key].description)

tsvlist = []
for ref in seqdict:
    for match in seqdict[ref][:-1]:
        tsvlist.append("\t".join([match[0], ref, match[1], seqdict[ref][-1]]))

with open(tsvfile[:-3] + "tax.tsv", "w") as out:
    out.write("\n".join(tsvlist))
