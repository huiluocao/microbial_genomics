#### plasmid & virus
conda activate genomad
genomad end-to-end --cleanup --splits 8 PA_ZSP7002.fna PA_ZSP7002_genomad /ibex/project/c2170/db/genomad_db/ &


#### prophage
## PHASTEST API
for file in *.fna
do
    a=${file%.fna};
    wget --post-file=$file "https://phastest.ca/phastest_api?contigs=1" -O $a'_phastest'
done

for file in *phastest
do
    a=$(awk -F'"' '{print $4}' $file);
    wget "https://phastest.ca/phastest_api?acc=$a" -O ${file}_results
    echo $a
done

#### phigaro

#### ICEfinder
1) split gbks to individual records:
from Bio import SeqIO
import os,sys
for file in os.listdir('./'):
    if file.endswith ('gbk'):
        with open(file,'r') as f:
            for rec in SeqIO.parse(f,'genbank'):
                with open('/ibex/scratch/projects/c2170/pa/gbks/ICEfinder_linux/gbk/%s.gbk'%rec.id,'w') as g:
                    SeqIO.write(rec,g,'genbank')
2) run ICEfinder                    
module load emboss/6.6.0/gnu-6.4.0
#perl ICEfinder_local.pl bcf_gbk.list
perl ICEfinder_local.pl bcf_gbk_p6.list


#### ISfinder

