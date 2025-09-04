import os, sys
import pandas as pd

os.chdir(r'/data2/share/home/liyi/CbpTE/PanTE_matrix') 

subge = sys.argv[1]
pop_mat = pd.read_csv(f'matrix_pop/panTE.Mix-SC-CC.{subge}.matrix', sep='\t', dtype=str, keep_default_na=False)
Cru_mat = pd.read_csv(f'matrix_19/panTE.C162-2.{subge}.matrix', sep='\t', dtype=str, keep_default_na=False)

pop_merge_df = pd.merge(pop_mat, Cru_mat, how='outer', on=['chr', 'start', 'end', 'id', 'ref'], suffixes=('_1', '_2'))
ref_merged = pop_merge_df[~pop_merge_df.iloc[:, 5:].isnull().any(axis=1)]
ref_merged1 = ref_merged.iloc[:, :5]
ref_merged1['pop_ref'] = ref_merged1['ref']
ref_merged1['pop_id'] = ref_merged1['id']
ref_merged1.to_csv(f'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix/int_file/panTE.Mix-SC-CC_C162-2.{subge}.map', sep='\t', index=None)