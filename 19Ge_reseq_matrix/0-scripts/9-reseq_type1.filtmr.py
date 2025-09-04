import os, sys
import pandas as pd

os.chdir(r'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix')
addr1 = f'/data2/share/home/liyi/CbpTE/PanTE_matrix/marker_19'

input_ge = sys.argv[1]

def processing_gt(input_df, subge='', ge_override=None):
    # 如果提供了 ge_override，则检查是否在输入数据的列中
    if ge_override:
        if ge_override in input_df.columns:
            ge_id = f'{ge_override}{subge}'
        else:
            if input_ge == "Cru" and (ge_override + ".Cru") in input_df.columns:
                ge_id = ge_override + ".Cru"
            elif input_ge == "Cor" and (ge_override + ".Cor") in input_df.columns:
                ge_id = ge_override + ".Cor"
            else:
                ge_id = f'{ge_override}{subge}'
    else:
        ge_id = f'{pop_ges_dict[input_pop][0]}{subge}'
    
    target_df = input_df.loc[:, [ge_id]]
    
    exploded = (
        target_df
        .assign(id_split=lambda x: x[ge_id].str.split(','))
        .explode('id_split')
        .query("id_split == '0' or id_split.str.contains('TE_|repeat_', na=False)")
        [['id_split'] + target_df.columns.difference(['id']).tolist()]
        .assign(orig_index=lambda df: df.index)
    )
    
    merged = pd.merge(exploded, gt_all_df, left_on='id_split', right_on='id', how='left')
    merged = merged.set_index('orig_index')
    
    result = merged.groupby(merged.index).apply(
        lambda group: group['gt'].dropna().value_counts().index[0]
        if not group['gt'].dropna().empty 
        else ('1/1' if (group['id_split'] == '0').any() else './.')
    )
    
    target_df = target_df.assign(gt=result.reindex(target_df.index, fill_value='./.').values)
    return target_df


# 4. 读取 marker 文件，并对 t1_df 进行复制和预处理
header1 = ['ge', 'chr', 'start', 'end', 'id', 'ref', 
           'C162-2.Cor','C98-1.Cor','C111.Cor','C72-3.Cor','C90-2.Cor','C157-3.Cor','C215-4.Cor','C261-3.Cor','C244.Cor','C58-2.Cor','C42-1.Cor','C37-1.Cor','C35-2.Cor',
           'C13-2', 'C34-3', 'C42-2',
           'C162-2.Cru','C98-1.Cru','C111.Cru','C72-3.Cru','C90-2.Cru','C157-3.Cru','C215-4.Cru','C261-3.Cru','C244.Cru','C58-2.Cru','C42-1.Cru','C37-1.Cru','C35-2.Cru',
           'C879', 'C86', '690']
t1_df = pd.read_csv(f'{addr1}/type1.simp.marker', sep='\t', dtype=str, header=None, names=header1)
t1_df1 = t1_df.copy()
del t1_df

# 5. 对 t1_df1 中每个 ge 列，当值为 '1' 且 ref 列等于该 ge 时，用 id 列替换该位置的值
ge_list = ['C162-2.Cor','C98-1.Cor','C111.Cor','C72-3.Cor','C90-2.Cor','C157-3.Cor','C215-4.Cor','C261-3.Cor','C244.Cor','C58-2.Cor','C42-1.Cor','C37-1.Cor','C35-2.Cor',
           'C13-2', 'C34-3', 'C42-2',
           'C162-2.Cru','C98-1.Cru','C111.Cru','C72-3.Cru','C90-2.Cru','C157-3.Cru','C215-4.Cru','C261-3.Cru','C244.Cru','C58-2.Cru','C42-1.Cru','C37-1.Cru','C35-2.Cru',
           'C879', 'C86', '690']
for col in ge_list:
    t1_df1.loc[(t1_df1[col] == '1') & (t1_df1['ref'] == col), col] = t1_df1.loc[(t1_df1[col] == '1') & (t1_df1['ref'] == col), 'id']
t1_df2 = t1_df1.iloc[:, :6].copy()

# 6. 根据 input_ge 参数读取对应的 reseq 列表
if input_ge == "Cbp":
    reseq_list = open('../Pangenie/0-scripts/reseq_Cbp_list', 'r').read().split('\n')[:-1]
elif input_ge == "Cor":
    reseq_list = open('../Pangenie/0-scripts/reseq_Cor_list', 'r').read().split('\n')[:-1]
elif input_ge == "Cru":
    reseq_list = open('../Pangenie/0-scripts/reseq_Cru_list', 'r').read().split('\n')[:-1]
else:
    print('check input genome')
    sys.exit(1)

# 7. 读取 population 信息，并构造 pop_ges_dict，用于后续根据不同种群选择对应列
pop_list = pd.read_csv('0-scripts/population.txt', sep='\t', dtype=str, header=None, names=['acc', 'pop'])    
pop_ges_dict = {'Cor': ['C13-2', 'C34-3', 'C42-2'], 'Cru': ['C879', 'C86', '690'], 
                'CC': ['C72-3','C90-2','C157-3','C215-4','C261-3'], 'SC': ['C111'], 'Mix': ['C162-2','C98-1'], 
                'CA': ['C58-2','C35-2','C37-1','C42-1'], 'EU': ['C244'], 'AN': ['C244']}

# 8. 对于每个 reseq 样本，处理 gt 数据并合并到 t1_df2 中
for reseq in reseq_list:
    # 根据 reseq_list 中第一个 acc 获取对应种群 input_pop
    input_pop = str(pop_list[pop_list['acc'] == reseq_list[0]]['pop'].iloc[0])
    
    gt_dfs = []
    # 若 input_ge 为 "Cru"，除了读取 pop_ges_dict[input_pop]（Cru组），还读取 pop_ges_dict['EU'] 的 vcf 文件
    if input_ge == "Cru":
        markers_for_vcf = pop_ges_dict[input_pop] + pop_ges_dict['EU']
    else:
        markers_for_vcf = pop_ges_dict[input_pop]
        
    # 循环读取各 marker 对应的 vcf 文件
    for marker in markers_for_vcf:
        vcf_path = f'../Pangenie/genotyping/{marker}/{marker}_{reseq}_gt.vcf'
        try:
            gt_vcf = pd.read_csv(vcf_path, sep='\t')
            gt_vcf = gt_vcf[(gt_vcf['gt'] != '0/1') & (gt_vcf['gq'] > 50)]
            gt_dfs.append(gt_vcf[['id', 'ref', 'gt']])
        except FileNotFoundError:
            print(f"Warning: File {vcf_path} not found")
            continue
    
    # 合并所有读取到的 gt 数据，并统一列名
    gt_all_df = pd.concat(gt_dfs, ignore_index=True)
    gt_all_df.columns = ['id', 'ref', 'gt']
    
    if input_ge == "Cbp":
        # 当 input_ge 为 "Cbp" 时，根据 input_pop 决定联合的 marker 组
        if input_pop == 'Mix':
            markers_for_vcf = pop_ges_dict['Mix'] + pop_ges_dict['SC'] + pop_ges_dict['CC']
        elif input_pop == 'SC':
            markers_for_vcf = pop_ges_dict['SC'] + pop_ges_dict['CC'] + pop_ges_dict['Mix']
        elif input_pop == 'CC':
            markers_for_vcf = pop_ges_dict['CC'] + pop_ges_dict['SC'] + pop_ges_dict['Mix']
        elif input_pop == 'CA':
            markers_for_vcf = pop_ges_dict['CA'] + pop_ges_dict['EU']
        elif input_pop in ['EU', 'AN']:
            markers_for_vcf = pop_ges_dict['EU'] + pop_ges_dict['CA']
        else:
            markers_for_vcf = pop_ges_dict[input_pop]
        
        # 针对 .Cor 分型
        gt_series_cor = pd.Series('./.', index=t1_df1.index)
        for marker in markers_for_vcf:
            # 调用时同时传入 subge 参数和 ge_override，
            # 期望处理的目标列为 marker + ".Cor"
            result_df = processing_gt(t1_df1, subge='.Cor', ge_override=marker)
            this_series = result_df['gt'].replace({'1/1': '0/0', '0/0': '1/1'})
            missing = gt_series_cor == './.'
            gt_series_cor[missing] = this_series[missing]
        t1_df2.loc[:, f'{reseq}_co'] = gt_series_cor
    
        # 针对 .Cru 分型
        gt_series_cru = pd.Series('./.', index=t1_df1.index)
        for marker in markers_for_vcf:
            # 目标列为 marker + ".Cru"
            result_df = processing_gt(t1_df1, subge='.Cru', ge_override=marker)
            this_series = result_df['gt'].replace({'1/1': '0/0', '0/0': '1/1'})
            missing = gt_series_cru == './.'
            gt_series_cru[missing] = this_series[missing]
        t1_df2.loc[:, f'{reseq}_cr'] = gt_series_cru
        
        print(f'{reseq} done')
            
    elif input_ge == "Cor":
        # 当 input_ge 为 "Cor" 时，联合 Cor 组和 CA 组
        markers_for_vcf = pop_ges_dict['Cor'] + pop_ges_dict['CA']
        gt_series = pd.Series('./.', index=t1_df1.index)
        for marker in markers_for_vcf:
            result_df = processing_gt(t1_df1, ge_override=marker)
            this_series = result_df['gt'].replace({'1/1': '0/0', '0/0': '1/1'})
            missing = gt_series == './.'
            gt_series[missing] = this_series[missing]
        t1_df2.loc[:, reseq] = gt_series
        print(f'{reseq} done')
        
    elif input_ge == "Cru":
        # 构造 markers_for_vcf 列表：先使用当前种群，再使用 EU 组
        markers_for_vcf = pop_ges_dict[input_pop] + pop_ges_dict['EU']
        gt_series = pd.Series('./.', index=t1_df1.index)
        for marker in markers_for_vcf:
            result_df = processing_gt(t1_df1, ge_override=marker)
            this_series = result_df['gt'].replace({'1/1': '0/0', '0/0': '1/1'})
            missing = gt_series == './.'
            gt_series[missing] = this_series[missing]
        t1_df2.loc[:, reseq] = gt_series
        print(f'{reseq} done')
        
    else:
        print('check input genome')
        sys.exit(1)

# 9. 根据 input_ge 决定输出格式和过滤条件，并保存最终结果
if input_ge == "Cbp":
    t1_df2_1 = pd.concat([t1_df2.iloc[:, :6], 
                          t1_df2.loc[:, t1_df2.columns.str.contains('_co')],
                          t1_df2.loc[:, t1_df2.columns.str.contains('_cr')]], axis=1)
    t1_df2_1.to_csv(f'filt/reseq_{input_ge}.type1.marker', sep='\t', index=None)
    t1_df3 = t1_df2_1[t1_df2_1.iloc[:, 6:].eq('./.').sum(axis=1) / (t1_df2_1.shape[1]-6) <= 0.1]
    t1_df3.to_csv(f'filt/reseq_{input_ge}.type1.marker.filtmr', sep='\t', index=None)
elif input_ge == "Cor" or input_ge == "Cru":
    t1_df2.to_csv(f'filt/reseq_{input_ge}.type1.marker', sep='\t', index=None)
    t1_df3 = t1_df2[t1_df2.iloc[:, 6:].eq('./.').sum(axis=1) / (t1_df2.shape[1]-6) <= 0.1]
    t1_df3.to_csv(f'filt/reseq_{input_ge}.type1.marker.filtmr', sep='\t', index=None)
else:
    print('check input genome')
    sys.exit(1)