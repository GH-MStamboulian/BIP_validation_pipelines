"""
create lod summary table
"""

import pandas as pd
import sys

tables_dir = sys.argv[1]#'/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/bip_3_5_3/lod_report_ALT_200128/'
bip_version = sys.argv[2]#'3.5.3'

table_5_part2_55_genes = 'Table5.LoD.part2.G360CDx.55genes.tsv'
table_5_part2_74_genes = 'Table5.LoD.part2.G360CDx.74genes.tsv'
table_5_55_genes = 'Table5.LoD.G360CDx.55genes.tsv'

variants = ['EGFR.T790M', 'EGFR.L858R', 'KRAS.G12V', 'BRAF.V600E', 'NRAS.Q61R', 'BRCA2.I2944F', 'BRCA1.S1140G', 'Panel-wide.SNV', 'EGFR.p.Glu746_Ala750del',
            'ERBB2.p.Ala775_Gly776insTyrValMetAla', 'EGFR.p.Ala767_Val769dup', 'BRCA2.p.Ser1982fs', 'BRCA1.p.Glu23fs', 'Panel-wide.INDEL', 'MET', 'ERBB2', 'CL_NTRK1_TPM3',
            'CL_RET_CCDC6', 'CL_ROS1_SLC34A2', 'CL_ALK_EML4'
           ]


table_5_part2_55_genes_df = pd.read_csv(tables_dir + table_5_part2_55_genes, sep = '\t')
table_5_part2_74_genes_df = pd.read_csv(tables_dir + table_5_part2_74_genes, sep = '\t')
table_5_55_genes_df = pd.read_csv(tables_dir + table_5_55_genes, sep = '\t')

LoD_sumamry_df = pd.DataFrame(columns = ['Variant', 'Variant_Type', 'LoD_5ng', 'LoD_30ng'])

for i, variant in enumerate(variants):
    if variant in ['BRCA2.I2944F', 'BRCA1.S1140G']:
        variant_type = table_5_part2_74_genes_df[(table_5_part2_74_genes_df['variant'] == variant) & (table_5_part2_74_genes_df['Input'] == '5 ng')]['variant_type'].to_string(header = False, index = False).strip()
        LoD_5ng = float(table_5_part2_74_genes_df[(table_5_part2_74_genes_df['variant'] == variant) & (table_5_part2_74_genes_df['Input'] == '5 ng')]['LoD'].to_string(header = False, index = False).strip())
        LoD_30ng = float(table_5_part2_74_genes_df[(table_5_part2_74_genes_df['variant'] == variant) & (table_5_part2_74_genes_df['Input'] == '30 ng')]['LoD'].to_string(header = False, index = False).strip())
        
    elif variant in ['Panel-wide.SNV', 'Panel-wide.INDEL']:
        variant_type = variant.split('.')[-1]
        LoD_5ng = float(table_5_55_genes_df[(table_5_55_genes_df['variant'] == 'aggregated') & (table_5_55_genes_df['variant_type'] == variant_type) & (table_5_55_genes_df['Input'] == '5 ng')]['LoD'].to_string(header = False, index = False).strip())
        LoD_30ng = float(table_5_55_genes_df[(table_5_55_genes_df['variant'] == 'aggregated') & (table_5_55_genes_df['variant_type'] == variant_type) & (table_5_55_genes_df['Input'] == '30 ng')]['LoD'].to_string(header = False, index = False).strip())        
    
    else:
        variant_type = table_5_part2_55_genes_df[(table_5_part2_55_genes_df['variant'] == variant) & (table_5_part2_55_genes_df['Input'] == '5 ng')]['variant_type'].to_string(header = False, index = False).strip()
        LoD_5ng = float(table_5_part2_55_genes_df[(table_5_part2_55_genes_df['variant'] == variant) & (table_5_part2_55_genes_df['Input'] == '5 ng')]['LoD'].to_string(header = False, index = False).strip())
        LoD_30ng = float(table_5_part2_55_genes_df[(table_5_part2_55_genes_df['variant'] == variant) & (table_5_part2_55_genes_df['Input'] == '30 ng')]['LoD'].to_string(header = False, index = False).strip())
        
    if (LoD_5ng < 0.1):
        LoD_5ng = round(LoD_5ng*100, 1)
    if (LoD_30ng < 0.1):
        LoD_30ng = round(LoD_30ng*100, 1)
    row = [variant, variant_type, LoD_5ng, LoD_30ng]      
    LoD_sumamry_df.loc[i] = row   
    
LoD_sumamry_df.to_csv(tables_dir + 'LoD_summaryTable_bip_'+ str(bip_version)+'.tsv', sep = '\t', index = False)
