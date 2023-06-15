#### run rgi on laptop
conda activate rgi
rgi load --card_json /Users/huiluocao/Desktop/apps/db/card.json --local
rgi database --version --local

for file in *fasta
do
    a=${file%.fasta}
    rgi main --input_sequence $file --output_file $a.rgi --input_type contig --local --clean
done
