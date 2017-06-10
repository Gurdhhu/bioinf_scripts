# Public scripts
Here I will upload some python3 scripts that I wrote initially for myself but then they appeared to be useful for other people working with nucleotide sequences.
## 1. annotate_blast_hits
### Motivation

When you download the results of blastn search from NCBI BLAST webpage, the table of blast hits does not contain any information about the matching sequences except their accession number.
It is not convenient to check all the matches manually, especially if you queried hundreds of sequences in a batch.

### Description

This script is intended for taxonomic annotation of blastn results saved in Hit Table CSV format where GenBank accession numbers are in the 4th column.
The output is an annotated CSV file "*_annotated.csv" with the following columns added:

* Record name
* Species name
* Date of update
* Reference
* Full taxonomy

### Running the script
This script **"annotate_blast_hits.py"** requires **python3** and package **Biopython** and runs from terminal on Linux.

I did not try to run it on Windows, but probably it should also work.

You can also download an executable for Linux that does not require any dependencies. The script with all dependencies was bundled into one file with [PyInstaller](http://www.pyinstaller.org/), so you can just save **"annotate_blast_hits"** and run it.



**Usage:**  <code>annotate_blast_hits.py [-h] csv email</code>

**Example usage:** <code>./annotate_blast_hits.py Documents/AlignmentHitTable.csv yourname@mail.com</code>

| arguments | description |
| --- | --- |
| csv | Pathway to csv-formatted Hit Table file with blastn results. Positional argument |
| email | Your e-mail address. Positional argument. It is important when you make many requests to NCBI databases because it helps NCBI to contact you if something goes wrong. Otherwise they can just silently block you |
| -h, --help | show help message and exit. Optional argument |


### Possible problems
* If any of GenBank accessions correspond to extremely large sequences, such as whole-genome bacterial sequences or complete chromosome sequences with a length of millions bp, it will run EXTREMELY slowly.

* The script writes some temporary files that are normally removed at the end. If the script is terminated manually, some files can stay unremoved.

### Contacts
For any questions and suggestions: ledum_laconicum(at)mail(dot)ru

If anybody will require additional annotation features (e.g., genes) or additional databases to be available (e.g., blastp), don't hesitate to contact me.
