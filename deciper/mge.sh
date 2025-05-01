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

#### ICEfinder


#### ISfinder

