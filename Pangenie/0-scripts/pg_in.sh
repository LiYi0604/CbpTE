#!/bin/bash
cd /data2/share/home/liyi/CbpTE/ReseqTE_matrix/Pangenie/genotyping
# source /data2/share/home/liyi/software/miniconda3/bin/activate pan-genome

ref_ge=$1
echo "processing ${ref_ge}"
mkdir ${ref_ge}/indexing

### .bed to .vcf
grep "Chr" /data2/share/home/liyi/CbpTE/TE_Anno/final_bed/${ref_ge}_TE.bed | cut -f 1-3 | bedtools merge -d 1 > ${ref_ge}/indexing/${ref_ge}.pos.bed

awk 'BEGIN {FS=OFS="\t"} {if ($2 <= 1) $2 = 2; print $1, $2-2, $2-1}' ${ref_ge}/indexing/${ref_ge}.pos.bed | \
  bedtools getfasta -fi /data2/share/home/liyi/CbpTE/Genome/fas/${ref_ge}.fasta -bed stdin -tab > ${ref_ge}/indexing/${ref_ge}.ref.bed

awk 'BEGIN {FS=OFS="\t"} {if ($2 <= 1) $2 = 2; print $1, $2-2, $3}' ${ref_ge}/indexing/${ref_ge}.pos.bed | \
  bedtools getfasta -fi /data2/share/home/liyi/CbpTE/Genome/fas/${ref_ge}.fasta -bed stdin -tab > ${ref_ge}/indexing/${ref_ge}.alt.bed

paste ${ref_ge}/indexing/${ref_ge}.pos.bed  ${ref_ge}/indexing/${ref_ge}.ref.bed  ${ref_ge}/indexing/${ref_ge}.alt.bed | \
  awk 'BEGIN {FS=OFS="\t"} {if ($2 <= 1) $2 = 2; print $1, $2-1, ".", $7, $5, ".", "PASS", ".", "GT", "1|1"}' >  ${ref_ge}/indexing/${ref_ge}.vcf

python /data2/share/home/liyi/CbpTE/ReseqTE_matrix/Pangenie/0-scripts/vcf_header.py ${ref_ge} 

cat ${ref_ge}/indexing/${ref_ge}.vcf.header ${ref_ge}/indexing/${ref_ge}.vcf > ${ref_ge}/indexing/${ref_ge}_final.vcf

echo "${ref_ge} finished"
