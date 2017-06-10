#!/usr/bin/env python3
from Bio import Entrez
from os import remove
from argparse import ArgumentParser, RawDescriptionHelpFormatter


parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description='''
This script is intended for annotation of blastn results saved in Hit Table CSV format where GenBank accession numbers are in the 4th column.
The output is an annotated CSV file "*_annotated.csv" with the following columns added:

Record name, Species, Date of update, Reference, Full taxonomy

Example usage: ./annotate_blast_hits.py Documents/AlignmentHitTable.csv yourname@mail.com

https://github.com/Gurdhhu/public_scripts
                                    ''')
parser.add_argument("csv",
                    type=str,
                    help="Pathway to csv-formatted Hit Table file with blastn results")
parser.add_argument("email",
                    type=str,
                    help="Your e-mail address. It is important when you make many requests to NCBI databases "
                         "because it helps NCBI to contact you if something goes wrong. Otherwise they can "
                         "just silently block you")

args = parser.parse_args()
file = args.csv
Entrez.email = args.email

with open(file, "r") as f:
    lines = f.readlines()

idlist = []  # list of sequence identifiers that will be sent to GenBank Nucleotide
filelist = []  # list of temporary xml files with descriptions of sequences
gblist = []  # list of annotations
csvout = []  # csv-formatted list of annotations

# Reading GenBank accessions from CSV file
for line in lines:
    line = line.strip().split(",")
    if len(line) > 1:  # contition to deal with blank lines
        idlist.append(line[3])


# Defining a function to get xml-formatted descriptions of accessions from GenBank Nucleotide
def get(ids, step):
    handle = Entrez.efetch(db="nucleotide", id=ids, retmode='xml')
    tmpfile = "".join([file, "_tmp", str(step), ".xml"])
    filelist.append(tmpfile)
    with open(tmpfile, "w") as out:  # Writing a temporary file that will be removed after parsing
        out.write(handle.read())


print("Fetching annotation for", len(idlist), "accessions from GenBank...")

'''
Requests to GenBank databases with efetch function cannot be longer than 200 accessions.
Do teal with more than 200 accessions, the list of accessions is split into shorter lists
that are queried sequentially.
'''
a = 0
b = 200
step = 1
try:
    for i in range(len(idlist)//200):
        ids_to_fetch = ",".join(idlist[a:b])
        a += 200
        b += 200
        get(ids_to_fetch, step)
        step += 1
    ids_to_fetch = ",".join(idlist[a:])
    get(ids_to_fetch, step)

    print("Parsing the result...")
    for tmp in filelist:
        with open(tmp, "r") as f:
            records = Entrez.parse(f)
            for i in records:
                gblist.append([i['GBSeq_definition'], i['GBSeq_organism'], i['GBSeq_update-date'],
                              ", ".join(i["GBSeq_references"][0]["GBReference_authors"]) + " " +
                               i["GBSeq_references"][0]["GBReference_title"],
                               "\",\"".join(i['GBSeq_taxonomy'].split("; "))])
        remove(tmp)  # removing temporary file
        filelist.remove(tmp)

except Exception as err:  # removing all temporary files in case of any error
    for i in filelist:
        remove(i)
    raise err

print("Writing annotations into file", "\"" + file[:-4] + "_annotated.csv\"")

for i in range(len(idlist)):  # formatting annotation into csv
    csvout.append(lines[i].strip() + ",\"" + "\",\"".join(gblist[i]) + "\"")

with open(file[:-4] + "_annotated.csv", "w") as out:
    out.write("\n".join(csvout))


print("Done!")
