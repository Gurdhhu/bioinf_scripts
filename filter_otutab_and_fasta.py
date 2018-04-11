#!/usr/bin/env python3
from Bio import SeqIO
import re
from argparse import ArgumentParser, RawDescriptionHelpFormatter

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description='''
This script filters an input TSV OTU abundance table produced by VSEARCH
(--otutabout option in --cluster_size command) leaving only the OTUs
that are present in the first column of another input TSV table.
Writes output TSV OTU abundance table adding information for every OTU about similarity to its best
hits in reference database and best hit consensus taxonomy. Output filename ends with ".myxonly.10.txt"
instead of ".txt" of the original input OTU table filename.
FASTA file with OTU centroids is also filtered accordingly, output FASTA ends with ".myxonly.10.fasta"
instead of ".fasta".
Example usage:
./filter_otutab_and_fasta.py all_acronyms_tax.csv \
    all.otus.991.myxodb.myxonly.10.tsv \
    all.otus.991.fasta \
    all.otutab.991.txt
                                    ''')
parser.add_argument("tsv_acr",
                    type=str,
                    help="Pathway to tsv-formatted table with acronyms and 7-level taxonomy for myxomycetes")
parser.add_argument("tsv_retain",
                    type=str,
                    help="Pathway to tsv-formatted hit table file with OTUs to retain")
parser.add_argument("fasta_filter",
                    type=str,
                    help="Pathway to FASTA file with OTU centroids to filter")
parser.add_argument("tsv_filter",
                    type=str,
                    help="Pathway to tsv-formatted hit table file with OTUs to filter")

args = parser.parse_args()
tsv_acr = args.tsv_acr
tsv_retain = args.tsv_retain
fasta_filter = args.fasta_filter
tsv_filter = args.tsv_filter

def common_start(*strings):
    """ Returns the longest common substring
        from the beginning of the `strings`
    """
    def _iter():
        for z in zip(*strings):
            if z.count(z[0]) == len(z):  # check all elements in `z` are the same
                yield z[0]
            else:
                return

    return ''.join(_iter())

with open(tsv_acr, "r") as f:
    lines = f.readlines()

taxdict = {}
for line in lines[1:]:
    acr, tax = line.split("\t")[0], "\t".join(line.split("\t")[1:8])
    taxdict[acr] = tax

with open(tsv_retain, "r") as f:
    lines = f.readlines()

centroiddict = {}
for line in lines:
    otuID, refID, sim = line.strip().split()
    otuID = otuID.split(";")[0]
    acr = refID.split("_")[0]
    if re.match("[A-Z]cf$", acr) or re.match("[A-Z]{3}[a-z]{2,3}cf$", acr) or re.match("[A-Z]{3}[a-z]{6}cf$", acr):
        acr = acr[:-2]
    if otuID not in centroiddict:
        centroiddict[otuID] = [[refID, sim, taxdict[acr]]]
    else:
        centroiddict[otuID].append([refID, sim, taxdict[acr]])

fasta_dict = SeqIO.index(fasta_filter, "fasta")
fasta_list = []
for key in fasta_dict:
    if key.split(";")[0] in centroiddict:
        fasta_list.append(">{}\n{}".format(key, str(fasta_dict[key].seq)))
with open(fasta_filter[:-5] + "myxonly.10.fasta", "w") as out:
    out.write("\n".join(fasta_list))

with open(tsv_filter, "r") as f:
    lines = f.readlines()

tsvlist = ["\t".join([lines[0].strip(),
                      "Best_match",
                      "Similarity",
                      "Kingdom",
                      "Phylum",
                      "Class",
                      "Order",
                      "Family",
                      "Genus",
                      "Species"])]
for line in lines[1:]:
    OTU = line.split("\t")[0]
    if OTU in centroiddict:
        if len(centroiddict[OTU]) == 1:
            tsvlist.append("{}\t{}".format(line.strip(), "\t".join(centroiddict[OTU][0])))
        else:
            taxlist = [i[2] for i in centroiddict[OTU]]
            consensus_tax = common_start(*taxlist)
            refidlist = [i[0] for i in centroiddict[OTU]]
            sim = centroiddict[OTU][0][1]
            tsvlist.append("{}\t{}".format(line.strip(), "\t".join([",".join(refidlist), sim, consensus_tax])))
with open(tsv_filter[:-3] + "myxonly.10.txt", "w") as out:
    out.write("\n".join(tsvlist))
