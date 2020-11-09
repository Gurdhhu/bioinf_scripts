# Bioinformatic scripts
Here I will upload some python3 scripts for dealing with nucleotide sequences that I wrote initially for myself but then they appeared to be useful for other people.
## 1. annotate_blast_hits
### Motivation

When you download the results of blast search from NCBI BLAST webpage, the table of blast hits does not contain any information about the matching sequences except their accession number.
It is not convenient to check all the matches manually, especially if you queried hundreds of sequences in a batch.

### Description

This script is intended for taxonomic annotation of **blast** results (blastn, tblastn, blastp or blastx) saved in **Hit Table CSV** format. It uses **efetch** function from **Bio.Entrez** package to get information about accessions from **GenBank Nucleotide (Nuccore)** or **Protein** databases.

**Input:** table in CSV or TSV format (separators: comma, tabulation or semicolon), where GenBank accession numbers are by default in the 2nd column (1 if counted from zero). The column can be changed via optional argument --column, -c (column number counting from zero).

The **output** is a file **"*_annotated.csv"** containing the original table with the following columns added:

* Record name
* Species name
* Date of update
* Reference
* Full taxonomy

### Running the script
This script **"annotate_blast_hits.py"** requires **python3** and package **Biopython** and runs from terminal on Linux.

I did not try to run it on Windows, but probably it should also work.

You can also download an **executable for Linux** that does not require any dependencies. The script with all dependencies was bundled into one file with [PyInstaller](http://www.pyinstaller.org/), so you can just save **"annotate_blast_hits"** and run it.

After downloading the script you should first make it executable:

<code>chmod +x annotate_blast_hits.py</code>

**Usage:**  <code>annotate_blast_hits.py [-h] [-c] csv db email</code>

**Example usage:** <code>./annotate_blast_hits.py example_input.csv n yourname@mail.com -c 1</code>

| arguments | description |
| --- | --- |
| csv | Pathway to csv-formatted Hit Table file with blastn results. Positional argument |
| db | For the output of **nucleotide blast** or **tblastn**, use <code>n</code>. For the output of **protein blast** or **blastx**, use <code>p</code>. Positional argument |
| email | Your e-mail address. It is important when you make many requests to NCBI databases because it helps NCBI to contact you if something goes wrong. Otherwise they can just silently block you. Positional argument |
| -c, --column | Specify the position of the column with GenBank accession numbers, counting from zero. Default value is 1. Optional argument |
| -h, --help | Show help message and exit. Optional argument |


### Possible problems
* If any of GenBank accessions correspond to extremely large sequences, such as whole-genome bacterial sequences or complete chromosome sequences with a length of millions bp, it will run EXTREMELY slowly.

* It was tested with hit tables with up to several thousand accessions. It usually works, but if it's a daytime in USA and NCBI servers are overloaded, the script can stop with a runtime error. In this case, it will try to fetch the data from GenBank once more automatically.

* The script writes some temporary files in the folder containing the hit table and they are normally removed at the end. If the script is terminated manually, some files can stay unremoved.

### Contacts
For any questions and suggestions: ledum_laconicum(at)mail(dot)ru

If anybody will require additional annotation features (e.g., genes), don't hesitate to contact me.




## 2. Additional scripts for VSEARCH pipeline
Several scripts are stored here as a supplementary material for my [publication](https://doi.org/10.1016/j.funeco.2018.11.006)
This is not the best organizational decision from me, but it's too late to change.
### Scripts:
add_silva_taxonomy.py
filter_fasta.py
filter_otutab_and_fasta.py
filter_tsv_by_other_tsv.py
