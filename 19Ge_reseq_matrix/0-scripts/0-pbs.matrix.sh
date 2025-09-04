#!/bin/bash
cd /data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix

line=$1
# nodes_name=$2
# for reseq in $(sed -n ${line}'p' 0-scripts/list_287); do
# for reseq in C111 C72-3 C90-2 C157-3 C215-4 C261-3; do
# for reseq in $(sed -n ${line}'p' ../Pangenie/0-scripts/reseq_Cbp_complement_list); do
for reseq in $(sed -n ${line}'p' 0-scripts/list_CA); do
    qsub -N $reseq.map -v reseq=${reseq} -l nodes=1:ppn=1 0-scripts/coor_switch.sh
done
