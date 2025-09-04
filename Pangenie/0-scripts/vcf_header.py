import os, sys

ge = sys.argv[1]
# fas = sys.argv[2]
infile = f"/data2/share/home/liyi/CbpTE/Genome/fas/{ge}.fasta.fai"
outfile = f"/data2/share/home/liyi/CbpTE/ReseqTE_matrix/Pangenie/genotyping/{ge}/indexing/{ge}.vcf.header"

with open(infile, 'r') as inf1:
    data = inf1.read().strip().split('\n')
inf1.close()

data1 = [d.split("\t")[:2] for d in data if "Chr" in d]
# print(data1)
header1 = f"""##fileformat=VCFv4.1\n##FILTER=<ID=PASS,Description="All filters passed">\n"""
header2 = f"""##INFO=<ID=ID,Number=A,Type=String,Description="Variant IDs.">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n"""
header3 = f"""##reference={ge}\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{ge}\n"""


if os.path.exists(outfile):
    if os.stat(outfile).st_size != 0:
        open(outfile, "w").close()

with open(outfile, 'a') as ouf1:
    ouf1.write(header1)
    for d1 in data1:
        ouf1.write(f"""##contig=<ID={d1[0]},length={d1[-1]}>\n""")
    ouf1.write(header2)
    ouf1.write(header3)
ouf1.close()