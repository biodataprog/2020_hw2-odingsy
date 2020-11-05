#!/usr/bin/env python3

import os, gzip, itertools

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

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))


numGene_s = 0
len_s = 0
dict_s = {'A':0, 'C':0, 'T':0, 'G':0}

numGene_m = 0
len_m = 0
dict_m = {'A':0, 'C':0, 'T':0, 'G':0}

codon = {}
# creating a dictionary for codon. 
for first in {'A','T','C','G'}:
    for second in {'A','T','C','G'}:
        for third in {'A','T','C','G'}:
            new_codon = first+second+third
            codon[new_codon] = 0
codon_s = codon.copy()
codon_m = codon.copy()

with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        numGene_s += 1
        len_s += len(seqstring)
        for bp in seq[1]:
            dict_s[bp] += 1
        for n in range(0, len(seq[1]), 3):
            codon_s[seq[1][n:n+3]] += 1

GC_s = (dict_s['G'] + dict_s['C'])/sum(dict_s.values())

with gzip.open(file2,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        numGene_m += 1
        len_m += len(seqstring)
        for bp in seq[1]:
            dict_m[bp] += 1
        for n in range(0, len(seq[1]), 3):
            codon_m[seq[1][n:n+3]] += 1

GC_m = (dict_m['G'] + dict_m['C'])/sum(dict_m.values())

print(f"Qn1: total number of genes in Salmonella is {numGene_s} while number of genes in Mycobacterium is {numGene_m}.")
print(f"Qn2: total length of genes in Salmonella is {len_s}bps while total length of genes in Mycobacterium is {len_m}bps.")
print(f"Qn3: GC content in Salmonella is {GC_s*100:.1f}% while GC content in Mycobacterium is {GC_m*100:.1f}%.")
print(f"Qn4: Total number of codon in Salmonella is {len_s/3:.0f} while total number of codon in Mycobacterium is {len_m/3:.0f}.")
print("Qn5:Codon frequency in two species:\nCodon\tSalmonella\tMycobacterium")
for codon in codon:
    print(codon, "\t", codon_s[codon], "\t", codon_m[codon])
