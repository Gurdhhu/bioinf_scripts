#!/usr/bin/env python3
from argparse import ArgumentParser, RawDescriptionHelpFormatter

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description='''
This script filters an input TSV table produced by VSEARCH (--usearch_global command) leaving only the OTUs
that are present in another input TSV table of the same format.
Writes output TSV table adding ".myxonly.tsv" instead of ".tsv" to the original input filename.
Example usage: ./filter_tsv_by_other_tsv.py all.otus.991.silva.tax.myxonly.tsv .all.otus.991.myxodb.tsv
                                    ''')
parser.add_argument("tsv_retain",
                    type=str,
                    help="Pathway to tsv-formatted hit table file with OTUs to retain")
parser.add_argument("tsv_filter",
                    type=str,
                    help="Pathway to tsv-formatted hit table file with OTUs to filter")

args = parser.parse_args()
tsv_retain = args.tsv_retain
tsv_filter = args.tsv_filter


with open(tsv_retain, "r") as f:
    lines = f.readlines()
idlist = [line.split()[0] for line in lines]

with open(tsv_filter, "r") as f:
    lines = f.readlines()
tsvlist = [line for line in lines if line.split()[0] in idlist]

with open(tsv_filter[:-3] + "myxonly.tsv", "w") as out:
    out.write("".join(tsvlist))
