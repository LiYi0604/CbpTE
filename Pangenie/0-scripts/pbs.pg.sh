#!/bin/bash
cd /data2/share/home/liyi/CbpTE/ReseqTE_matrix/Pangenie

### Genome list


### Generate input files
for ref_ge in C*-*; do
	mkdir -p genotyping/$ref_ge
	bash 0-scripts/pg_in.sh $ref_ge
done

### pangenie-index submit script
# ref_ge=$1
# nodes_name=$1
for ref_ge in C*-*; do
	qsub -N pg.${ref_ge}_index -v ref_ge=${ref_ge} -l nodes=1:ppn=20,mem=20gb 0-scripts/pg.sh
done

### pangenie submit script
line=$1
# nodes_name=$2
for ref_ge in C*-*; do
	for reseq_ge in $(sed -n ${line}'p' 0-scripts/reseq_Cbp_complement_list); do
		qsub -N pg.${ref_ge}_${reseq_ge} -v ref_ge=${ref_ge},reseq_ge=${reseq_ge} -l nodes=1:ppn=10,mem=20gb \
		0-scripts/pg.sh
	done
done

### check completion
for ge in C*-*; do
	python 0-scripts/check_undone.py $ge
done
