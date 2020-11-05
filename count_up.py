#!/usr/bin/env python3

# this is a python script template
# this next line will download the file using curl

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")
    
geneNum = 0 # set a variable: number of genes for Qn2
geneLength = 0 # set a variable: gene length for Qn3
chrLength = 0 # set a variable: chromosome length in Qn4

with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
        if row[2] == 'gene':
            geneNum += 1
            newLength = int(row[4]) - int(row[3])
            geneLength += newLength
    print(f"Qn2: number genes is {geneNum}.")
    print(f"Qn3: total gene length is {geneLength}.")

# Qn4 and Qn5
with gzip.open(fasta, "rt") as f:
    pairs = aspairs(f)
    seqs = dict(pairs)
for k,v in seqs.items():
    chrLength += len(v)
print(f"Qn4: total length of genome is {chrLength}.")
print(f"Qn5: percentage of the genome which is coding {geneLength/chrLength*100:.1f}%(assuming no overlapping genes.)")
