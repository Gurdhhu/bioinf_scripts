#!/usr/bin/env python3
from Bio import SeqIO
from argparse import ArgumentParser, RawDescriptionHelpFormatter

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description='''
This script filters an input FASTA file with OTU centroids retaining only OTUs
that are present in filtered TSV OTU abundance table.
Output FASTA ends with ".60.fasta" instead of ".fasta".
Example usage:
./filter_fasta.py all.otutab.98.myxonly.10.60.txt all.otus.98.myxonly.10.fasta
                                    ''')
parser.add_argument("tsv_retain",
                    type=str,
                    help="Pathway to tsv-formatted OTU table with OTUs to retain")
parser.add_argument("fasta_filter",
                    type=str,
                    help="Pathway to FASTA file with OTU centroids to filter")

args = parser.parse_args()
tsv_retain = args.tsv_retain
fasta_filter = args.fasta_filter

with open(tsv_retain, "r") as f:
    lines = f.readlines()

OTU_list = [i.split("\t")[0] for i in lines[1:]]

fasta_list = []
with open(fasta_filter, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        if record.id.split(";")[0] in OTU_list:
            fasta_list.append(">{}\n{}".format(record.id, str(record.seq)))

with open(fasta_filter[:-5] + "60.fasta", "w") as out:
    out.write("\n".join(fasta_list))
