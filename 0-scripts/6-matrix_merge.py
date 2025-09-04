#!/usr/bin/env python3
import os, sys
import pandas as pd

os.chdir(r'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix')
addr_matrix = r'/data2/share/home/liyi/CbpTE/PanTE_matrix'

# 读取基础 matrix
matrix = pd.read_csv(f'{addr_matrix}/matrix_19/panTE.C162-2.matrix.anno', sep='\t', dtype=str, keep_default_na=False)
matrix = matrix.iloc[:, :6]
matrix.columns = ['chr', 'start', 'end', 'id', 'ref', 'class']

def merge_reseq(matrix, list_file, out_file):
    """按 list_file 中的样本依次 merge 并重命名 gt 列，最后输出到 out_file"""
    with open(list_file) as f:
        samples = [l.strip() for l in f if l.strip()]
    merged = matrix
    for s in samples:
        m = pd.read_csv(f'int_file/reseq_map/{s}.mat', sep='\t', dtype=str)
        merged = pd.merge(merged, m, how='left',
                          on=['chr','start','end','id','ref'])
        merged = merged.rename(columns={'gt': s})
        print(f'{s} done')
    merged.to_csv(out_file, sep='\t', index=False)

# 根据输入调用
input_ge = sys.argv[1]
if   input_ge == "Cbp":
    merge_reseq(matrix,
                '../Pangenie/0-scripts/reseq_Cbp_list',
                'filt/reseq_Cbp.mat')
elif input_ge == "Cor":
    merge_reseq(matrix,
                '../Pangenie/0-scripts/reseq_Cor_list',
                'filt/reseq_Cor.mat')
elif input_ge == "Cru":
    merge_reseq(matrix,
                '../Pangenie/0-scripts/reseq_Cru_list',
                'filt/reseq_Cru.mat')
else:
    print('check input genome')
