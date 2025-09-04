import os
import sys
import pandas as pd

os.chdir(r'/data2/share/home/liyi/CbpTE/ReseqTE_matrix/19Ge_reseq_matrix/filt')

def dataframe_to_vcf(df):
    vcf_body = []
    for _, row in df.iterrows():
        chrom = row["chr"]
        pos = row["start"]
        te_id = row["id"]
        ref = "A"
        alt = "T"
        qual = "."
        filter_status = "PASS"
        info = f"REF={row['ref']};TE={te_id}"
        format_field = "GT"
        genotypes = "\t".join(row[6:].values)
    
        vcf_line = f"{chrom}\t{pos}\t{te_id}\t{ref}\t{alt}\t{qual}\t{filter_status}\t{info}\t{format_field}\t{genotypes}"
        vcf_body.append(vcf_line)
    
    vcf_content = "\n".join(vcf_body)
    return vcf_content

def main(df, df_name):
    names = {
        'Cbp_df': 'reseq_Cbp',
        'Cor_df': 'reseq_Cor',
        'Cru_df': 'reseq_Cru',
        'all_df': 'reseq_all'
    }

    # df = df.drop(index=90251)
    # df.to_csv(f'{names[df_name]}.matrix', sep='\t', index=False)
    acc = "\t".join(df.columns.tolist()[6:])
    vcf_header = f"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=ID,Number=A,Type=String,Description="Variant IDs.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{acc}
"""
    vcf_content = dataframe_to_vcf(df)
    with open(f"{names[df_name]}.vcf", "w") as f:
        f.write(vcf_header)
        f.write(vcf_content)

specie = sys.argv[1]
df_dict = {}

if specie == 'all':
    cbp_df = pd.read_csv('reseq_Cbp.mat', sep='\t', dtype=str, keep_default_na=False)
    cor_df = pd.read_csv('reseq_Cor.mat', sep='\t', dtype=str, keep_default_na=False)
    cru_df = pd.read_csv('reseq_Cru.mat', sep='\t', dtype=str, keep_default_na=False)

    header = cbp_df.iloc[:, :6]
    all_df = pd.concat([header, cbp_df.iloc[:, 6:], cor_df.iloc[:, 6:], cru_df.iloc[:, 6:]], axis=1)
    all_df.to_csv(f'reseq_all.mat', sep='\t', index=False)
    df_dict['all_df'] = all_df
else:
    df_dict[f'{specie}_df'] = pd.read_csv(f'reseq_{specie}.mat', sep='\t', dtype=str, keep_default_na=False)

if __name__ == "__main__":
    if specie == 'all':
        main(df_dict['all_df'], 'all_df')
    else:
        main(df_dict[f'{specie}_df'], f'{specie}_df')
