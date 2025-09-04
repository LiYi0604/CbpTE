import os
import numpy as np
import pandas as pd

os.chdir(r'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix/filt')
addr1 = f'/data2/share/home/liyi/CbpTE/PanTE_matrix/marker_19'

### only for type2 & type3
reseqs = ['Cbp', 'Cor', 'Cru']
for r in reseqs:
    reseq_mat = pd.read_csv(f'reseq_{r}.mat.filtmr', sep='\t', dtype=str)
    if r == 'Cbp':
        reseq_Cor_mat = reseq_mat.iloc[:76891, :]
        reseq_Cru_mat = reseq_mat.iloc[76891:, :]

        types_cor_list = ['type2.cor', 'type3.scor']
        for t in types_cor_list:
            allfilt_mat = pd.read_csv(f'{addr1}/{t}.marker', sep='\t', header=None, dtype=str).iloc[:, :5]
            allfilt_mat.columns = reseq_mat.columns.tolist()[:5]
            allfilt_reseq_mat = pd.merge(allfilt_mat, reseq_Cor_mat, how='left', on=allfilt_mat.columns.tolist())
            allfilt_reseq_mat = allfilt_reseq_mat[~allfilt_reseq_mat.iloc[:, 5:].isnull().all(axis=1)]
            allfilt_reseq_mat.to_csv(f'reseq_{r}.{t}.marker.filtmr', sep='\t', index=None)
            print(f'reseq_{r}.{t}.marker.filtmr done')
        types_cru_list = ['type2.cru', 'type3.scru']
        for t in types_cru_list:
            allfilt_mat = pd.read_csv(f'{addr1}/{t}.marker', sep='\t', header=None, dtype=str).iloc[:, :5]
            allfilt_mat.columns = reseq_mat.columns.tolist()[:5]
            allfilt_reseq_mat = pd.merge(allfilt_mat, reseq_Cru_mat, how='left', on=allfilt_mat.columns.tolist())
            allfilt_reseq_mat = allfilt_reseq_mat[~allfilt_reseq_mat.iloc[:, 5:].isnull().all(axis=1)]
            allfilt_reseq_mat.to_csv(f'reseq_{r}.{t}.marker.filtmr', sep='\t', index=None)
            print(f'reseq_{r}.{t}.marker.filtmr done')

    elif r == 'Cor':
        types_list = ['type2.cor', 'type3.cor']
        for t in types_list:
            allfilt_mat = pd.read_csv(f'{addr1}/{t}.marker', sep='\t', header=None, dtype=str).iloc[:, :5]
            allfilt_mat.columns = reseq_mat.columns.tolist()[:5]
            allfilt_reseq_mat = pd.merge(allfilt_mat, reseq_mat, how='left', on=allfilt_mat.columns.tolist())
            allfilt_reseq_mat = allfilt_reseq_mat[~allfilt_reseq_mat.iloc[:, 5:].isnull().all(axis=1)]
            allfilt_reseq_mat.to_csv(f'reseq_{r}.{t}.marker.filtmr', sep='\t', index=None)
            print(f'reseq_{r}.{t}.marker.filtmr done')

    elif r == 'Cru':
        types_list = ['type2.cru', 'type3.cru']
        for t in types_list:
            allfilt_mat = pd.read_csv(f'{addr1}/{t}.marker', sep='\t', header=None, dtype=str).iloc[:, :5]
            allfilt_mat.columns = reseq_mat.columns.tolist()[:5]
            allfilt_reseq_mat = pd.merge(allfilt_mat, reseq_mat, how='left', on=allfilt_mat.columns.tolist())
            allfilt_reseq_mat = allfilt_reseq_mat[~allfilt_reseq_mat.iloc[:, 5:].isnull().all(axis=1)]
            allfilt_reseq_mat.to_csv(f'reseq_{r}.{t}.marker.filtmr', sep='\t', index=None)
            print(f'reseq_{r}.{t}.marker.filtmr done')