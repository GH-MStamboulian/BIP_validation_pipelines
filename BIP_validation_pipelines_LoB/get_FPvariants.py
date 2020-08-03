"""
script to combine and filter variants for false positives
"""

import pandas as pd
import sys

#tables_dir = "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/bip_3_5_0/lob_report_190215/"

tables_dir = sys.argv[1]

variant_types = ["SNV", "INDEL", "CNV", "FUSION"]

filters = ["cfDNA.call", "gDNA.call", "MDL.call"]

cols = ['variant_type', 'chrom', 'position', 'gene', 'transcript_id', 'exon', 'mut_nt', 'mut_aa', 'mut_key', 'cfDNA.call', 'gDNA.call', "MDL.call", 'mut_cnt', 'pool.maf', 'cfDNA.maf', 'gDNA.maf']
FP_variants = pd.DataFrame(columns = cols)

empty_tables = []
row_ctr = 0
for variant_type in variant_types:
    variant_df = pd.read_csv(tables_dir + "Table2.All."+variant_type+".tsv", sep = '\t')    
    if variant_df.shape[0] == 0:
        empty_tables.append(variant_type)
    else:
        for i, row in variant_df.iterrows():
            new_row = [row[col] if col in variant_df.columns else "NaN" for col in cols]
            #new_row.insert(0, variant_type)
            FP_variants.loc[row_ctr] = new_row
            row_ctr += 1;
        
for variant_type in empty_tables:
    new_row = ['NaN']*15
    new_row.insert(0, variant_type)
    FP_variants.loc[row_ctr] = new_row
    row_ctr += 1
    
FP_variants.to_csv(tables_dir + "Table2.All.VARIANTS.tsv", sep = '\t', index = False)

FP_variants_filtered = FP_variants[(FP_variants['cfDNA.call'] == False) & (FP_variants['gDNA.call'] == False) & (FP_variants['MDL.call'] == False)]
remaining_variants = list(set(variant_types).difference(set(list(FP_variants_filtered['variant_type']))))
row_ctr = len(FP_variants_filtered)
for variant_type in remaining_variants:
    new_row = ['NaN']*15
    new_row.insert(0, variant_type)
    FP_variants_filtered.loc[row_ctr] = new_row
    row_ctr += 1

FP_variants_filtered = FP_variants_filtered.drop(columns = ['cfDNA.call', 'gDNA.call', 'MDL.call'])
FP_variants_filtered.to_csv(tables_dir + "Table2.All.VARIANTS.filtered.tsv", sep = '\t', index = False)
