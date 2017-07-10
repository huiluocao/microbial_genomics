import os,sys
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

input_fna='C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/temp/alignment2/ERS1021086-cps.fna.fas'
##strain_list = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/hku_s3_cps/list.txt'
dest = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/cps_sp3/hku_s3_cps/hku_s3_cps_complete/'

with open(input_fna,"rU") as f:
    header = []
    sequences = []
    for record in SeqIO.parse(f,'fasta'):
        header.append(record.id)
        sequence = record.seq[1:]
        sequences.append(str(sequence))
##        print(sequences)
output_fna = open('C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/temp/alignment2/ERS1021086-cps.fna','a')
output_fna.write('>'+header[0].split('|')[0]+'\n')
for seq in sequences:
    output_fna.write(seq)

f.close()
output_fna.close()
