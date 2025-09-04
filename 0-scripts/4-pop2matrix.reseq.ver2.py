import os
import sys
import pandas as pd

os.chdir(r'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix/')
addr_matrix = r'/data2/share/home/liyi/CbpTE/PanTE_matrix'

input_reseq = sys.argv[1]
pop_list = pd.read_csv('0-scripts/population.txt', sep='\t', dtype=str, header=None, names=['acc', 'pop'])
input_pop = pop_list[pop_list['acc'] == input_reseq]['pop'].values[0]

if input_pop in ['EU', 'AN']:
    input_pops = 'EU-CA'
elif input_pop in ['CA']:
    input_pops = 'CA-EU'
elif input_pop in ['Mix', 'SC', 'CC']:
    input_pops = 'Mix-SC-CC'
elif input_pop == 'Cor' or input_pop == 'Cru':
    input_pops = input_pop.lower()
else:
    raise ValueError(f'chech input pop: {input_pop}')

pop_ges_dict = {'EU-CA': ['C244'],
                'CA-EU': ['C35-2','C37-1','C42-1','C58-2'],
                'Mix-SC-CC': ['C162-2','C98-1','C111','C157-3','C215-4','C261-3','C72-3','C90-2'],
                'cor': ['C13-2','C34-3','C42-2'],
                'cru': ['C879', 'C86', '690']}

matrix_Cor = pd.read_csv(f'{addr_matrix}/matrix_19/panTE.C162-2.Cor.matrix', sep='\t', dtype=str)
matrix_Cru = pd.read_csv(f'{addr_matrix}/matrix_19/panTE.C162-2.Cru.matrix', sep='\t', dtype=str)

if input_pops in ['Mix-SC-CC', 'EU-CA', 'CA-EU']:
    matrix_pop_Cor = pd.read_csv(f'{addr_matrix}/matrix_pop/panTE.{input_pops}.Cor.matrix', sep='\t', dtype=str, usecols=[0,3,4])
    matrix_pop_Cru = pd.read_csv(f'/{addr_matrix}/matrix_pop/panTE.{input_pops}.Cru.matrix', sep='\t', dtype=str, usecols=[0,3,4])
    matrix_map_Cor = pd.read_csv(f'int_file/panTE.{input_pops}_C162-2.Cor.map', sep='\t', dtype=str)
    matrix_map_Cru = pd.read_csv(f'int_file/panTE.{input_pops}_C162-2.Cru.map', sep='\t', dtype=str)
    reseq = input_reseq

    cat_df = pd.DataFrame()
    for ge in pop_ges_dict[input_pops]:
        gt_vcf_df = pd.read_csv(f'../Pangenie/genotyping/{ge}/{ge}_{reseq}_gt.vcf', sep='\t', dtype=str, usecols=[0,3,5])
        merge_ref_Cor_df = pd.merge(matrix_pop_Cor[matrix_pop_Cor['ref'] == f'{ge}.Cor'], gt_vcf_df, on=['chr', 'id'], how='left')
        merge_ref_Cru_df = pd.merge(matrix_pop_Cru[matrix_pop_Cru['ref'] == f'{ge}.Cru'], gt_vcf_df, on=['chr', 'id'], how='left')
        cat_df = pd.concat([cat_df, merge_ref_Cor_df, merge_ref_Cru_df])
    # cat_df.to_csv(f'/data2/share/home/liyi/TEs/Pairwise.new/matrix.new/matrix.map/reseq_test/{reseq}.mat.bed', sep='\t', index=False)
    map_merge_Cor = pd.merge(matrix_map_Cor, cat_df, how='left', left_on=['pop_ref', 'pop_id'], right_on=['ref', 'id'],  suffixes=('', '_map'))
    map_merge_Cru = pd.merge(matrix_map_Cru, cat_df, how='left', left_on=['pop_ref', 'pop_id'], right_on=['ref', 'id'],  suffixes=('', '_map'))
    map_merge_Cor_merge = pd.merge(matrix_Cor, map_merge_Cor, how='left', on=['chr', 'start', 'end', 'id', 'ref']).iloc[:, [0,1,2,3,4,-1]]
    map_merge_Cru_merge = pd.merge(matrix_Cru, map_merge_Cru, how='left', on=['chr', 'start', 'end', 'id', 'ref']).iloc[:, [0,1,2,3,4,-1]]
    map_merge = pd.concat([map_merge_Cor_merge, map_merge_Cru_merge])

    reseq_mat_df = pd.read_csv(f'int_file/reseq_bed/{reseq}.mat.bed', sep='\t', dtype=str)
    rest_mat = pd.merge(map_merge[map_merge['gt'].isnull()].iloc[:, :5], reseq_mat_df, how='left', on=['chr', 'id', 'ref'])
    final_mat = pd.merge(map_merge.iloc[:, :5], pd.concat([map_merge[~map_merge['gt'].isnull()], rest_mat]), how='left', on=['chr', 'start', 'end', 'id', 'ref'])
    final_mat1 = final_mat.replace('0/1', './.').fillna('./.')
    final_mat1['gt'] = final_mat1['gt'].replace({'1/1': '0/0', '0/0': '1/1'})
    final_mat1.to_csv(f'int_file/reseq_map/{reseq}.mat', sep='\t', index=None, na_rep='./.')
    print(f'{reseq} done')
        
elif input_pops in ['cor', 'cru']:
    subge = input_pops.title()
    matrix_pop = pd.read_csv(f'{addr_matrix}/matrix_pop/panTE.{input_pops}.matrix', sep='\t', dtype=str, usecols=[0,3,4])
    matrix_map = pd.read_csv(f'int_file/panTE.{input_pops}_C162-2.{subge}.map', sep='\t', dtype=str)
    reseq = input_reseq
    cat_df = pd.DataFrame()
    for ge in pop_ges_dict[input_pops]:
        gt_vcf_df = pd.read_csv(f'../Pangenie/genotyping/{ge}/{ge}_{reseq}_gt.vcf', sep='\t', dtype=str, usecols=[0,3,5])
        merge_ref_df = pd.merge(matrix_pop[matrix_pop['ref'] == f'{ge}'], gt_vcf_df, on=['chr', 'id'], how='left')
        cat_df = pd.concat([cat_df, merge_ref_df])
    # cat_df.to_csv(f'/data2/share/home/liyi/TEs/Pairwise.new/matrix.new/matrix.map/reseq_test/{reseq}.mat.bed', sep='\t', index=False)
    map_merge_di = pd.merge(matrix_map, cat_df, how='left', left_on=['pop_ref', 'pop_id'], right_on=['ref', 'id'],  suffixes=('', '_map'))
    map_merge_Cor_merge = pd.merge(matrix_Cor, map_merge_di, how='left', on=['chr', 'start', 'end', 'id', 'ref']).iloc[:, [0,1,2,3,4,-1]]
    map_merge_Cru_merge = pd.merge(matrix_Cru, map_merge_di, how='left', on=['chr', 'start', 'end', 'id', 'ref']).iloc[:, [0,1,2,3,4,-1]]
    map_merge = pd.concat([map_merge_Cor_merge, map_merge_Cru_merge])
    
    reseq_mat_df = pd.read_csv(f'int_file/reseq_bed/{reseq}.mat.bed', sep='\t', dtype=str)
    rest_mat = pd.merge(map_merge[map_merge['gt'].isnull()].iloc[:, :5], reseq_mat_df, how='left', on=['chr', 'id', 'ref'])
    final_mat = pd.merge(map_merge.iloc[:, :5], pd.concat([map_merge[~map_merge['gt'].isnull()], rest_mat]), how='left', on=['chr', 'start', 'end', 'id', 'ref'])
    final_mat1 = final_mat.replace('0/1', './.').fillna('./.')
    final_mat1['gt'] = final_mat1['gt'].replace({'1/1': '0/0', '0/0': '1/1'})
    final_mat1.to_csv(f'int_file/reseq_map/{reseq}.mat', sep='\t', index=None, na_rep='./.')
    print(f'{reseq} done')