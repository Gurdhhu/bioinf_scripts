#!/usr/bin/env python3
from Bio import Entrez
from os import remove
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from urllib.error import HTTPError


parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description='''
This script is intended for annotation of blast results saved in Hit Table CSV format.
The output is an annotated CSV file "*_annotated.csv" with the following columns added:
Record name, Species, Date of update, Reference, Full taxonomy
https://github.com/Gurdhhu/bioinf_scripts
                                    ''',
                        epilog='Example usage: ./annotate_blast_hits.py example_input.csv n yourname@mail.com')
parser.add_argument("csv",
                    type=str,
                    help="Pathway to csv-formatted Hit Table file with blast results")
parser.add_argument("db",
                    choices=['n', 'p'],
                    help="For the output of nucleotide blast or tblastn, use \'n\'. "
                         "For the output of protein blast or blastx, use \'p\'.")
parser.add_argument("email",
                    type=str,
                    help="Your e-mail address. It is important when you make many requests to NCBI databases "
                         "because it helps NCBI to contact you if something goes wrong. Otherwise they can "
                         "just silently block you")
parser.add_argument("-c", "--column",
                    type=int,
                    help="Number of the column with GenBank accession numbers, counting from zero. Default value is 1")

args = parser.parse_args()
file = args.csv
db = args.db
Entrez.email = args.email
col = 1
if args.column != None:
    col = args.column

with open(file, "r") as f:
    lines = f.readlines()

idlist = []  # list of sequence identifiers that will be sent to GenBank Nucleotide
filelist = []  # list of temporary xml files with descriptions of sequences
gblist = []  # list of annotations
csvout = []  # csv-formatted list of annotations
sep = "," # column separator
ext = "csv"

# Determining column separator
if len(lines[0].strip().split("\t")) > 1:
    sep = "\t"
    ext = "tsv"
    print("Column separator: tabulation")
elif len(lines[0].strip().split(",")) > 1:
    sep = ","
    print("Column separator: comma")
elif len(lines[0].strip().split(";")) > 1:
    sep = ";"
    print("Column separator: semicolon")
else:
    print("Column separator: comma")

# Reading GenBank accessions from CSV file
for line in lines:
    line = line.strip().split(sep)
    if len(line) > 1:  # condition to deal with blank lines
        idlist.append(line[col])

dbdict = {"n": "nuccore","p": "protein"}

# Defining a function to get xml-formatted descriptions of accessions from GenBank Nucleotide
def get(ids, step):
    handle = Entrez.efetch(db=dbdict[db], id=ids, retmode='xml')
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
        get(ids_to_fetch, step)
        a += 200
        b += 200
        step += 1
    ids_to_fetch = ",".join(idlist[a:])
    if ids_to_fetch != "":
        get(ids_to_fetch, step)


except HTTPError as err:
    for i in filelist:
        remove(i)    
    if err.code == 400:
        print("HTTP Error 400: Bad Request.\n"
              f"First query element was: {idlist[0]}\n"
              "Wrong database? Use 'n' for Nuccore and 'p' for Protein.\n"
              "Wrong column with GenBank accession numbers? "
              f"Currently, column # {col} is used. You can change it using --column or -c")
        exit()
    else:
        raise err

except RuntimeError:
    print(f"RuntimeError occurred during step {step-1}. Trying one more time...") #list the ID numbers of current step - e.g. 201-400

    try:
        for i in range(len(idlist[a:])//200):
            ids_to_fetch = ",".join(idlist[a:b])
            a += 200
            b += 200
            get(ids_to_fetch, step)
            step += 1
        ids_to_fetch = ",".join(idlist[a:])
        if ids_to_fetch != "":
            get(ids_to_fetch, step)

    except Exception as err:  # removing all temporary files in case of any error
        for i in filelist:
            remove(i)
        raise err

except Exception as err:  # removing all temporary files in case of any error
    for i in filelist:
        remove(i)
    raise err


print("Parsing the result...")

try:
    for tmp in filelist:
        with open(tmp, "r") as f:
            records = Entrez.parse(f)
            for i in records:
                if "GBSeq_references" in i and "GBReference_authors" in i["GBSeq_references"][0]:
                    gblist.append([i['GBSeq_definition'],
                                   i['GBSeq_organism'],
                                   i['GBSeq_update-date'],
                                   ", ".join(i["GBSeq_references"][0]["GBReference_authors"]) +
                                   " " + i["GBSeq_references"][0]["GBReference_title"],
                                   "\",\"".join(i['GBSeq_taxonomy'].split("; "))])
                else:
                    gblist.append([i['GBSeq_definition'],
                                   i['GBSeq_organism'],
                                   i['GBSeq_update-date'],
                                   "no reference",
                                   "\",\"".join(i['GBSeq_taxonomy'].split("; "))])

    for tmp in filelist:
        remove(tmp)  # removing temporary file

except Exception as err:  # removing all temporary files in case of any error
    for i in filelist:
        remove(i)
    raise err


print("Writing annotations into file", f"\"{file}_annotated.{ext}\"")

for i in range(len(idlist)):  # formatting annotation into csv
    annotation = f"\"{sep}\"".join(gblist[i])
    annotated_line = f"{lines[i].strip()}{sep}\"{annotation}\""
    csvout.append(annotated_line)

with open(f"{file}_annotated.{ext}", "w") as out:
    out.write("\n".join(csvout))


print("Done!")