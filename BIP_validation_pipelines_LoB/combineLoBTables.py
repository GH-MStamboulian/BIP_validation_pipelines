"""
script that will take FP variants from both bip versions and will merge and combine them to one table
"""

import pandas as pd
import sys


#table_1_f = "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/bip_3_5_0/lob_report_190215/Table2.All.VARIANTS.filtered.tsv"
#bip_version_1 = "3.5.0"
#table_2_f = "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/bip_3_5_3/lob_report_190215/Table2.All.VARIANTS.filtered.tsv"
#bip_version_2 = "3.5.3"
#fileName = '/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/LoB1_bip_' + bip_version_1 + '_LoB2_bip_' + bip_version_2 + '.pdf'

table_1_f = sys.argv[1]
bip_version_1 = sys.argv[2]
table_2_f = sys.argv[3]
bip_version_2 = sys.argv[4]
out_dir = sys.argv[5]

fileName = out_dir + 'LoB1_bip_' + bip_version_1 + '_LoB2_bip_' + bip_version_2 + '.pdf'
out_fname = out_dir + 'LoB1_bip_' + bip_version_1 + '_LoB2_bip_' + bip_version_2 + '.tsv'

cols = ['variant_type', 'chrom', 'position', 'gene', 'transcript_id', 'exon', 'mut_nt', 'mut_aa', 'mut_key', 'mut_cnt', 'pool.maf', 'cfDNA.maf', 'gDNA.maf']

table_1_df = pd.read_csv(table_1_f, sep = '\t')
table_2_df = pd.read_csv(table_2_f, sep = '\t')

merged_df = pd.merge(table_1_df, table_2_df, how = "inner", on = cols)

merged_df["bip_" + bip_version_1] = [True]*len(merged_df)
merged_df["bip_" + bip_version_2] = [True]*len(merged_df)

table_1_only_rows = list()
for i, row in table_1_df.iterrows():
    if row["mut_aa"] not in list(merged_df["mut_aa"]) and row["mut_key"] not in list(merged_df["mut_key"]):
        row = list(row)
        row.extend([True, False])
        table_1_only_rows.append(row)
    
    
table_2_only_rows = list()
for i, row in table_2_df.iterrows():
    if row["mut_aa"] not in list(merged_df["mut_aa"]) and row["mut_key"] not in list(merged_df["mut_key"]):
        row = list(row)
        row.extend([True, False])
        table_1_only_rows.append(row)

row_ctr = len(merged_df)
if (table_1_only_rows):
    for row in table_1_only_rows:
        merged_df.loc[row_ctr] = row
        row_ctr += 1

if (table_2_only_rows):
    for row in table_2_only_rows:
        merged_df.loc[row_ctr] = row
        row_ctr += 1

merged_df = merged_df.sort_values(by = "variant_type")   

merged_df.to_csv(out_fname, sep = '\t', index = False)

data = list()
data.append(list(merged_df.columns))

for i, row in merged_df.iterrows():
    data.append(list(merged_df.loc[i]))

from reportlab.platypus import SimpleDocTemplate, Paragraph
from reportlab.lib.pagesizes import letter
from reportlab.platypus import Table
from reportlab.platypus import TableStyle
from reportlab.lib import colors
from reportlab.lib.units import mm, inch
from reportlab.lib.styles import getSampleStyleSheet


styles = getSampleStyleSheet()
styleN = styles['Heading1']

pagesize = (400 * mm, 210 * mm)

pdf = SimpleDocTemplate(
    fileName,
    pagesize=pagesize
)

table = Table(data)

# add style


style = TableStyle([
    #('BACKGROUND', (0,0), (3,0), colors.green),
    ('TEXTCOLOR',(0,0),(-1,0),colors.black),

    ('ALIGN',(0,0),(-1,-1),'CENTER'),

    ('FONTNAME', (0,0), (-1,0), 'Courier'),
    ('FONTSIZE', (0,0), (-1,0), 12),
    
    ('BOTTOMPADDING', (0,0), (-1,0), 12),

    ('BACKGROUND',(0,1),(-1,-1),colors.beige),
])
table.setStyle(style)


# 3) Add borders
ts = TableStyle(
    [
    ('BOX',(0,0),(-1,-1),2,colors.black),

    #('LINEBEFORE',(2,1),(2,-1),2,colors.red),
    #('LINEABOVE',(0,2),(-1,2),2,colors.green),

    ('GRID',(0,0),(-1,-1),2,colors.black),
    ]
)
table.setStyle(ts)

P = Paragraph("Table summarizing the comparison of LoB (False positive) variants between two BIP versions", styleN)

elems = []
elems.append(P)
elems.append(table)

pdf.build(elems)
