declare dir="/share/analysis/Carlo/miRNA/DMSO"
declare input="/share/analysis/Carlo/miRNA/DMSO/miRNA_name"
SUFFIX_in="_R1.fastq"

for fq1 in $input/*.fastq
do 
    bn=$(basename $fq1)
    sample=${bn}
    out=$dir/2_mirge2/${sample:0:-${#SUFFIX_in}}

#gunzip $fq1 

miRge2.0 annotate -s $fq1 -d miRBase -pb /share/tools/bowtie-1.1.1 \
 -lib "/share/analysis/Carlo/miRNA" -sp human -cpu 8 -ad illumina

echo $fq1
echo "Annotation done"

done
