import os,re
import sys
from collections import OrderedDict

input_blast = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/hku_s3_cps/hku-s3-cps-gene-blastn.txt'
input_fas = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/hku_s3_cps/hku_sp3.fna'
output_fas = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/hku_s3_cps/hku_sp3_wzh_gene.fna'

fasta = open(input_fas,'U')
fasta_dict = {}
for line in fasta:
    line = line.strip()
    if line =='':
        continue
    if line.startswith('>'):
        header = line.lstrip('>')
        seqname = header.split(' ')[0]
    ##        print seqname
        fasta_dict[seqname]= ''
    else:
        fasta_dict[seqname] += line
##    print(seqname)
fasta.close()

out_seq = []
with open(input_blast) as f:
    lines = f.readlines()
    gene_coordinates = []
    for line in lines:
##        if line.split("\t")[1] in ["SPNOXC_RS01945", "SPNOXC_RS01950", "SPNOXC_RS01960", "SPNOXC_RS01965"]:
        if line.split("\t")[1] == "SPNOXC_RS01925":
            items = line.split("\t")
            gene_coordinates.append([items[6],items[7]])
##            print (gene_coordinates)
            s = min(int(items[6]),int(items[7]))-2000
            e = max(int(items[6]),int(items[7]))
            print(s,e)
            header = '>'+items[0]+'_'+items[1]+':'+str(s) + '-' + str(e) + '\n'
            out_seq.append(header)
            seq= fasta_dict[items[0]][s:e]
            out_seq.append(seq+'\n')
            
with open(output_fas,'w') as e:
    for line in out_seq:
        e.write(line)
            
            
            
