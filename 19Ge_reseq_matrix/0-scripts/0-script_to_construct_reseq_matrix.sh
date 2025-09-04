#!/bin/bash
cd /data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix

### After completing PanGenie genotyping

### 1. check genotyping quality
for i in 0; do
    echo "$i" >> gt_error
    for file in ../genotyping/*/*.log2; do
        grep -P "Computed kmer abundance peak: $i$" $file > /dev/null && echo "$file" >> gt_error
    done
done
mkdir -p ../check
for reseq_ge in $(sed -n '1,150p' ../0-scripts/reseq_all_list); do
    grep "parameter" ../genotyping/*/*_${reseq_ge}_histogram.histo > ../check/${reseq_ge}.parameter.check
done

### 2. a filteration for genotyping vcf file and merge to annotation bed file 
# get merge detail of reference input file
mkdir -p anno_bed
for ge in C*-*; do
    cut -f 1-3,5,7 /data2/share/home/liyi/CbpTE/TE_Anno/final_bed/${ge}_TE.bed | \
        bedtools merge -d 1 -c 2,3,4,5 -o collapse,collapse,collapse,collapse | \
        awk 'BEGIN{FS=OFS="\t"} $2 = $2-1 {print$0}'> anno_bed/${ge}.anno.merge.bed
done

for reseq in $(cat 0-scripts/list_323); do
    python 0-scripts/coor_switch.e.py $reseq
    echo "$reseq done"
done
# pbs ver
bash 0-scripts/pbs.matrix.sh 1,287   # Cbp: 1,188; Cor: 1,48 ; Cru: 1,54

### 3. population matrix for mapping gt vcf files
python 0-scripts/map.Mix-SC-CC_C162-2.py Cor
python 0-scripts/map.Mix-SC-CC_C162-2.py Cru
python 0-scripts/map.EU-CA_C162-2.py Cor
python 0-scripts/map.EU-CA_C162-2.py Cru
python 0-scripts/map.CA-EU_C162-2.py Cor
python 0-scripts/map.CA-EU_C162-2.py Cru
python 0-scripts/map.cor_C162-2.Cor.py
python 0-scripts/map.cru_C162-2.Cru.py

### 4. mapping gt vcf to population matrix to get a single resequence matrix
for reseq in $(cat 0-scripts/list_323); do
    python 0-scripts/pop2matrix.reseq.ver2.py $reseq
    echo "$reseq done"
done
# pbs ver
bash 0-scripts/pbs.matrix.sh 1,37   # 125,134   # 135,223   # 857,

### 5. performance test for Pangenie and matrix synthesis method
jupyter precision/0-scripts/precision.ipynb

### 6. concatenate single resequence matrixes to resequences matrix
python 0-scripts/matrix_merge.py Cbp
python 0-scripts/matrix_merge.py Cor
python 0-scripts/matrix_merge.py Cru

### 7. convert resequences matrix to VCF format
# 7.1 conver to vcf format
python 0-scripts/reseq_mat_to_vcf.py all
python 0-scripts/reseq_mat_to_vcf.py Cbp
python 0-scripts/reseq_mat_to_vcf.py Cor
python 0-scripts/reseq_mat_to_vcf.py Cru
# 7.2 calculate missing rate
vcftools --vcf filt/reseq_all.vcf --missing-site --out filt/reseq_all.site_missing
vcftools --vcf filt/reseq_Cbp.vcf --missing-site --out filt/reseq_Cbp.site_missing
vcftools --vcf filt/reseq_Cor.vcf --missing-site --out filt/reseq_Cor.site_missing
vcftools --vcf filt/reseq_Cru.vcf --missing-site --out filt/reseq_Cru.site_missing
# 7.3 filter out missing site
python 0-scripts/reseq_vcf_filter.py all
python 0-scripts/reseq_vcf_filter.py Cbp
python 0-scripts/reseq_vcf_filter.py Cor
python 0-scripts/reseq_vcf_filter.py Cru

### 8. statistics
# 8.1 boxplot of TE number in populations
python reseq_all_sample_count.py
# 8.2 figure of missing rate check
jupyter reseq_all_sample_count.ipynb
# 8.3 histogram of TE counts for each sample within each species.
jupyter reseq_mat_site_freq.ipynb

### 9. TE marker in reseq matrix
# 9.1 type2 and type3 marker
python 0-scripts/reseq_marker.py
python 0-scripts/reseq_marker_filtmr.py
# 9.2 type1 marker
python 0-scripts/reseq_type1.filtmr.py Cru
python 0-scripts/reseq_type1.filtmr.py Cor
python 0-scripts/reseq_type1.filtmr.py Cbp
cut -f 1-903 filt/reseq_Cbp.type1.marker > filt/reseq_Cbp.type1.cor.marker
cut -f 1-903 filt/reseq_Cbp.type1.marker.filtmr > filt/reseq_Cbp.type1.cor.marker.filtmr
cut -f 1-6,904-1800 filt/reseq_Cbp.type1.marker > filt/reseq_Cbp.type1.cru.marker
cut -f 1-6,904-1800 filt/reseq_Cbp.type1.marker.filtmr > filt/reseq_Cbp.type1.cru.marker.filtmr

# 9.3 TE marker statistics
jupyter 0-scripts/reseq_type1_test.ipynb




