# after assembling using CANU, contig can be polished using Illumina reads

module load bwa
module load samtools
module load pilon

bwa index S2_canu/S2.contigs.fasta
bwa mem S2_canu/S2.contigs.fasta hku_reads/rawdata/S2_1.fastq.gz hku_reads/rawdata/S2_2.fastq.gz > S2.sam
samtools view -S S2.sam -b -o S2.bam
samtools sort -o S2.index.bam S2.bam
samtools index S2.index.bam
java -jar $PILON --genome S2_canu/S2.contigs.fasta --frags S2.index.bam 
