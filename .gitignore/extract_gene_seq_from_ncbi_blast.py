import os,re
import sys
from collections import OrderedDict

input_blast = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/ERR/ERR_blastn_cleaned.txt'
input_fas = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/ERR/ERR_s3_all.fna'
output_fas = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/ERR/ERR_all_cps.fna'


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

d={}
with open(input_blast) as f:
    strain_lines = f.readlines()
    for line in strain_lines:
        d.setdefault(line.split('\t')[0],[]).append(line.split('\t')[1])

out_seq = []
out_contig_info = open('C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/ERR/ERR_contig_info.txt','w')

with open(input_blast) as f:
    strain_lines = f.readlines()
    for key in d:
        contig_coordinates=[]
        contig_info=[]
        for line in strain_lines:
            if line.startswith(key):
                items = line.split('\t')
##                print items
                contig_coordinates.extend([items[6],items[7]])
                contig_info.extend([items[1],items[6],items[7]])
##        print contig,contig_coordinates,len(contig_coordinates)
        contig_coordinates = [int(x) for x in contig_coordinates]
        s = min(contig_coordinates)
        e = max(contig_coordinates)
        header = '>' + key + ':' + str(s) + '-' + str(e) + '\n'
        out_seq.append(header)
        print (header)
        seq= fasta_dict[key][s:e]
        out_seq.append(seq+'\n')
        out_contig_info.write('\n'+key+':' + str(s) + '-' + str(e) + '\n')
        out_contig_info.writelines(["%s\t" %info for info in contig_info])

with open(output_fas,'w') as e:
    for line in out_seq:
        e.write(line)

out_contig_info.close()
