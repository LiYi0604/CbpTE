import os
import sys
import pandas as pd
import traceback

try:
    os.chdir(r'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix/')
    reseq = sys.argv[1]

    ges = ['C162-2', 'C111', 'C157-3', 'C215-4', 'C244', 'C261-3', 'C35-2', 'C37-1', 'C42-1', 'C58-2', 'C72-3', 'C90-2', 'C98-1', '690', 'C86', 'C879', 'C13-2', 'C34-3', 'C42-2']

    for ge in ges:
        try:
            merge_df = pd.read_csv(f'anno_bed/{ge}.anno.merge.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'start_m', 'end_m', 'id_m', 'class_m'])
            
            vcf_df = pd.read_csv(f'../Pangenie/genotyping/{ge}/{ge}_{reseq}_genotyping.vcf', comment='#', sep='\t', header=None)
            vcf_df = vcf_df.iloc[:, [0, 1, 7, 9]]
            vcf_df.columns = ['chr', 'start', 'qual', 'info']
            vcf_df[['gt', 'gq', 'gl', 'kc']] = vcf_df['info'].str.split(':', expand=True)
            vcf_df = vcf_df.iloc[:, [0, 1, 4, 5]]
            vcf_df = vcf_df[~(vcf_df['gt'] == '.')]

            vcf_merge_df = pd.merge(vcf_df, merge_df, on=['chr', 'start'], how='left')
            vcf_merge_df = vcf_merge_df[~(vcf_merge_df.isnull().any(axis=1))]
            vcf_merge_df = vcf_merge_df[(vcf_merge_df['gt'] != '0/1') & (vcf_merge_df['gq'].astype(int) >= 50)]
            
            multi_df = vcf_merge_df[(vcf_merge_df['start_m'].str.contains(','))]
            columns = ['start_m', 'end_m', 'id_m', 'class_m']
            multi_df.loc[:, columns] = multi_df.loc[:, columns].copy().apply(lambda x: x.str.split(','))
            multi_df_expanded = multi_df.explode(columns)
            
            expanded_vcf_df = pd.concat([vcf_merge_df[~vcf_merge_df['start_m'].str.contains(',')], multi_df.explode(['start_m', 'end_m', 'id_m', 'class_m'])])
            expanded_vcf_df['start_m'] = expanded_vcf_df['start_m'].astype(int) - 1
            expanded_vcf_df['ref'] = ge
            expanded_vcf_df = expanded_vcf_df.iloc[:, [0, 5, 6, 7, 8, 2, 3, 9]]
            expanded_vcf_df.columns = ['chr', 'start', 'end', 'id', 'class', 'gt', 'gq', 'ref']

            expanded_vcf_df.to_csv(f'../Pangenie/genotyping/{ge}/{ge}_{reseq}_gt.vcf', sep='\t', index=False)
            # print(f'{ge}_{reseq}_gt.vcf done')

        except Exception as e:
            print(f"Error processing {ge}_{reseq}_gt.vcf: {str(e)}")
            traceback.print_exc()

except Exception as e:
    print(f"A critical error occurred: {ge}_{reseq}_gt.vcf: {str(e)}")
    traceback.print_exc()
