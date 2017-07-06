#!/usr/bin/env python

import os
import sys,re
from Bio import SeqIO

src_dir = "C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/script/temp/"
des_dir = "C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/script/temp/cleaned/"
fasta_list = []
for file in os.listdir(src_dir):
    if file.endswith('fas'):
        fasta_list.append(file)

for fasta in fasta_list:
    print (fasta)
    in_fasta = open(os.path.join(src_dir,fasta))
    out_fasta = open(os.path.join(des_dir,fasta),'a') 
    for line in in_fasta:
        if line.startswith(">"):
            line = line.split("|")
            out_fasta.write(line[0]+"\n")
        else:
            out_fasta.write(line)
    in_fasta.close()
    out_fasta.close()
