import os, sys
import pandas as pd

os.chdir(r'/data2/share/home/liyi/CbpTE/PanTE_matrix') 

subge = sys.argv[1]
pop_mat = pd.read_csv(f'matrix_pop/panTE.EU-CA.{subge}.matrix', sep='\t', dtype=str, keep_default_na=False)
Cru_mat = pd.read_csv(f'matrix_19/panTE.C162-2.{subge}.matrix', sep='\t', dtype=str, keep_default_na=False)
ge1_list = pop_mat['ref'].unique().tolist()

pop_merge_df = pd.merge(pop_mat, Cru_mat, how='outer', on=['chr', 'start', 'end', 'id', 'ref'], suffixes=('_1', '_2'))
ref_merged = pop_merge_df[~pop_merge_df.iloc[:, 5:].isnull().any(axis=1)]
no_ref_merge = pop_merge_df[pop_merge_df.iloc[:, 10:].isnull().all(axis=1)]
ref_merge_remain = pop_merge_df[pop_merge_df.iloc[:, 5:10].isnull().all(axis=1)]

no_ref_merge_df = pd.DataFrame()
for ge_i in range(0, len(ge1_list)):
    no_ref_merge_id = no_ref_merge[no_ref_merge['ref'] == ge1_list[ge_i]]['id']
    pattern = '(' + '|'.join([f"{id}," for id in no_ref_merge_id]) + ')'
    no_ref_merge_found = ref_merge_remain[ref_merge_remain[f'{ge1_list[ge_i]}_2'].str.contains(pattern)]
    no_ref_merge_found['pop_id'] = no_ref_merge_found[f'{ge1_list[ge_i]}_2'].str.extract(pattern, expand=False).str.rstrip(',')
    no_ref_merge_found['pop_ref'] = ge1_list[ge_i]
    # print(len(no_ref_merge_found))
    no_ref_merge_df = pd.concat([no_ref_merge_df, no_ref_merge_found])
last_col = no_ref_merge_df.pop(no_ref_merge_df.columns[-1])
no_ref_merge_df.insert(5, last_col.name, last_col)
last_col = no_ref_merge_df.pop(no_ref_merge_df.columns[-1])
no_ref_merge_df.insert(6, last_col.name, last_col)

no_ref_merge_df1 = no_ref_merge_df.loc[no_ref_merge_df.iloc[:, :5].drop_duplicates().index, :].iloc[:, :7]
ref_merged1 = ref_merged.iloc[:, :5]
ref_merged1['pop_ref'] = ref_merged1['ref']
ref_merged1['pop_id'] = ref_merged1['id']

merge_df = pd.concat([ref_merged1, no_ref_merge_df1])
merge_df1 = merge_df[~merge_df.iloc[:, :5].duplicated(keep='first')]
merge_df1.to_csv(f'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix/int_file/panTE.EU-CA_C162-2.{subge}.map', sep='\t', index=None)