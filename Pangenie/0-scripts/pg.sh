#PBS -l walltime=100:00:00
#PBS -S /bin/bash
#PBS -k oe

#!/bin/bash
cd /data2/share/home/liyi/CbpTE/ReseqTE_matrix/Pangenie/genotyping
# source /data2/share/home/liyi/software/miniconda3/bin/activate pan-genome

# ref_ge=$1
# echo "processing ${ref_ge}"

### PanGenie-index
/data2/share/home/liyi/software/pangenie/build/src/PanGenie-index \
        -t 20 -k 21 -r /data2/share/home/liyi/CbpTE/Genome/fas/${ref_ge}.fasta \
        -v ${ref_ge}/indexing/${ref_ge}_final.vcf \
        -o ${ref_ge}/indexing/${ref_ge}`` \
        2>> ${ref_ge}/indexing/${ref_ge}.log

### PanGenie
reseq=$(awk -v reseq_ge="$reseq_ge" '$1 == reseq_ge {print $2, $3}' /data2/share/home/liyi/CbpTE/ReseqTE_matrix/Pangenie/0-scripts/matched_all)
file1=$(echo $reseq | awk '{print $1}')
file_type1=$(file -b --mime-type $file1)
if [[ "$file_type1" == "application/x-gzip" ]]; then
    if [[ "$file1" == *.tar.gz ]]; then
        cat_cmd1="tar -xzOf"
    else
        cat_cmd1="zcat"
    fi
elif [[ "$file_type1" == "text/plain" ]]; then
    cat_cmd1="cat"
else
    echo "Unknown file type: $file_type1 for $file1. Exiting."
    exit 1
fi

# PanGenie Command
/data2/share/home/liyi/software/pangenie/build/src/PanGenie \
    -f ${ref_ge}/indexing/${ref_ge} -o ${ref_ge}/${ref_ge}_${reseq_ge} \
    -i <($cat_cmd1 $file1 $file2) -s ${ref_ge}_${reseq_ge} -j 20 -t 20 2>> ${ref_ge}/${ref_ge}_${reseq_ge}.log2

