import os,sys
import re
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

src_dir = 'C:/Users/Huiluo/Desktop/HKU-microbiology/pneumococcus/genomes/serotype3_ref_45genomes'

files = os.listdir(src_dir)
print(files)

for file in files:
    with open(os.path.join(src_dir,file)) as f:
        fna = f.read()
        new_header = '>'+file.split('.')[0]+'_'
        fna = fna.replace('>', new_header)
##    print (fna)
        with open(os.path.join(src_dir,file.split('.')[0]+'.fna'),'w') as g:
            g.write(fna)
            
