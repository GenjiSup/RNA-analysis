#!/bin/bash
declare output="/share/analysis/Carlo/azathioprine"
declare input="/share/analysis/Carlo/trimming/1_Trimmed/"
declare indx="/share/analysis/Carlo/genome.index.star/reference.star"
# declare gtf="Homo_sapiens.GRCh38.100.gtf"

# prepare rsem reference
#rsem-prepare-reference --gtf "/share/analysis/Carlo/genomes/Homo_sapiens.GRCh38.109.gtf" --star -p 16 \
#"/share/analysis/Carlo/genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa" reference.star
 
#rsem: assinging the reads
for fq1 in $input/*_R1*.fastq.gz
do
    bn=$(basename $fq1)
    sample=${bn%_R1*.fastq.gz}
    fq2=$input/${sample}_R2*.fastq.gz
    out=$output/azathioprine${sample}

rsem-calculate-expression --star --star-gzipped-read-file --paired-end -p 16 --no-bam-output $fq1 $fq2 \
  $indx $out
echo $fq1
echo $fq2
echo $out
done
