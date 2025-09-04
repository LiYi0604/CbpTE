import os, sys
import pandas as pd

os.chdir(r'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix/filt')

specie = sys.argv[1]

mr_file = pd.read_csv(f'reseq_{specie}.site_missing.lmiss', sep='\t')
all_mat = pd.read_csv(f'reseq_{specie}.mat', sep='\t')
all_mat1 = all_mat.drop(mr_file[mr_file['F_MISS'] > 0.1].index.tolist(), axis=0)
all_mat1.to_csv(f'reseq_{specie}.mat.filtmr', sep='\t', index=None)