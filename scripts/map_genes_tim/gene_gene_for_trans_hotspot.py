import pandas as pd
import pybedtools
from pybedtools import BedTool


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", action="store", help="input gemma file", required=True)
parser.add_argument("-p", "--pvalue", action="store", type= float, help="pvalue threshold", required=True)
parser.add_argument("-o", "--output", action="store", help="output file name", required=True)
args = parser.parse_args()
##usage
"""
python gene_gene_for_trans_hotspot.py -p 0.05 -i columbia_tim.txt -o columbia.out
"""

with open(args.input) as f:
    xx = pd.read_csv(f, sep = '\t')
    # xx[['all_gene_count']] = xx[['all_gene_count']].apply(pd.to_numeric)
    # print xx
    xx_filt = xx[xx.P <= args.pvalue]
    # print xx_filt
    # xx_filt[['chr', 'snp']] = xx_filt['trans_eqtn'].str.split('-', expand=True)
    # # print
    # xx_filt[['snp']] = xx_filt[['snp']].astype('int')
    # xx_filt[['chr']] = xx_filt[['chr']].astype('int')



    xx_filt['start'] = xx_filt['BP'] - 1

    xx_clean = xx_filt[['CHR','start','BP', 'SNP']]

    xx_clean['CHR']=xx_clean['CHR'].astype(str)
    xx_clean['BP']=xx_clean['BP'].astype(str)
    xx_clean['chr_snp_name'] = xx_clean[['CHR','BP','SNP']].apply(lambda x: '-'.join(x), axis=1)

    xx_clean_sort = xx_clean.sort_values(['CHR', 'start'])

    xx_clean_sort_fin = xx_clean_sort[['CHR','start','BP', 'chr_snp_name']]

    best_at_hit = {}
    at_symbol = {}
    at_defline = {}
    at_tf_family = {}
    pt_tf_family = {}

    with open("format_annotation_with_tf.txt") as f:
        f.readline()
        for a in f:
            linea = a.strip().split('\t')
            best_at_hit[linea[0]]=linea[1]
            at_symbol[linea[0]]=linea[2]
            at_defline[linea[0]]=linea[3]
            at_tf_family[linea[0]]=linea[4]
            pt_tf_family[linea[0]]=linea[5]

    trichocarpa3_1_mod = 'trichocarpa_3.0_gene_mod.bed'

    # for a in xx_clean_sort:
        # print a

    x_gene_gene_df = pd.DataFrame(xx_clean_sort).to_csv(index = False, header = False, sep = '\t')

    # print x_gene_gene_df

    x_gene_gene_bedtool = pybedtools.BedTool(x_gene_gene_df, from_string = True)
    slop_gene_gene= x_gene_gene_bedtool.slop(b=10000, g='Ptrichocarpa_444_v3.0_mod.genome')
    combined_df = pd.DataFrame()

    # print x_gene_gene_bedtool

    slop_gene_gene_int = slop_gene_gene.intersect(trichocarpa3_1_mod, wb = True, nonamecheck = True)
    slop_gene_gene_df = pd.read_csv(slop_gene_gene_int.fn, sep = '\t', names = ['a','b','c','d','e','f','g','h','i','j'])
    slop_gene_gene_df = slop_gene_gene_df.drop_duplicates(subset='h', keep='first')

    slop_gene_gene_df['SNP10kb_AtGeneID'] = slop_gene_gene_df.h.map(best_at_hit)
    slop_gene_gene_df['SNP10kb_AtGeneName'] = slop_gene_gene_df.h.map(at_symbol)
    slop_gene_gene_df['SNP10kb_At_Anno'] = slop_gene_gene_df.h.map(at_defline)
    slop_gene_gene_df['SNP10kb_At_TFfamily'] = slop_gene_gene_df.h.map(at_tf_family)
    slop_gene_gene_df['SNP10kb_Pt_TFfamily'] = slop_gene_gene_df.h.map(pt_tf_family)

    # slop_gene_gene_df['Target_gene_At_GeneID'] = slop_gene_gene_df.d.map(best_at_hit)
    # slop_gene_gene_df['Target_gene_At_GeneName'] = slop_gene_gene_df.d.map(at_symbol)
    # slop_gene_gene_df['Target_gene_At_Anno'] = slop_gene_gene_df.d.map(at_defline)
    # slop_gene_gene_df['Target_gene_At_TFfamily'] = slop_gene_gene_df.d.map(at_tf_family)
    # slop_gene_gene_df['Target_gene_Pt_TFfamily'] = slop_gene_gene_df.d.map(pt_tf_family)

    slop_gene_gene_df_fin = slop_gene_gene_df[['d', 'h','SNP10kb_AtGeneID','SNP10kb_AtGeneName','SNP10kb_At_Anno','SNP10kb_At_TFfamily','SNP10kb_Pt_TFfamily']]
    slop_gene_gene_df_fin.columns = ['SNP','SNP10kb_GeneID','SNP10kb_AtGeneID','SNP10kb_AtGeneName','SNP10kb_At_Anno','SNP10kb_At_TFfamily','SNP10kb_Pt_TFfamily']




slop_gene_gene_df_fin.to_csv(args.output, sep = '\t', index = False, header = True)
    # print xx_clean
