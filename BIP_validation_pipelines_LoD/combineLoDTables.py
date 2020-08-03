"""
construct comparison table from BIP
"""

import pandas as pd
import sys

out1_dir = sys.argv[1]#'/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/LoD_comparisons/out_1_bip_3.5.3/'
out2_dir = sys.argv[2]#'/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/LoD_comparisons/out_2_bip_3.5.3/'

version1 = sys.argv[3]#'3.5.3'
version2 = sys.argv[4]#'3.5.3'

out_dir = sys.argv[5]

df1 = pd.read_csv(out1_dir + 'LoD_summaryTable_bip_' + str(version1)+ '.tsv', sep = '\t')
df2 = pd.read_csv(out2_dir + 'LoD_summaryTable_bip_' + str(version2)+ '.tsv', sep = '\t')

merged_df = pd.DataFrame(columns = ['Variant', 'Variant_Type', 'LoD_5ng_bip:' + str(version1), 'LoD_5ng_bip:' + str(version2), 'LoD_30ng_bip:' + str(version1), 'LoD_30ng_bip:' + str(version2)])

for (i, row1), (j, row2) in zip(df1.iterrows(), df2.iterrows()):
    variant = row1['Variant']
    variant_type = row1['Variant_Type']
    LoD_5ng_1 = row1['LoD_5ng']
    LoD_5ng_2 = row2['LoD_5ng']
    LoD_30ng_1 = row1['LoD_30ng']    
    LoD_30ng_2 = row2['LoD_30ng']
    merged_df.loc[i] = [variant, variant_type, LoD_5ng_1, LoD_5ng_2, LoD_30ng_1, LoD_30ng_2]

merged_df.to_csv(out_dir + "/" + 'LoD1_bip_'+str(version1)+'_LoD2_bip_'+str(version2)+'.tsv', sep = '\t', index = False)

fileName = out_dir + "/" + 'LoD1_bip_'+str(version1)+'_LoD2_bip_'+str(version2)+'.pdf'

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

pagesize = (320 * mm, 210 * mm)

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

P = Paragraph("Table summarizing the comparison of LoD values between two BIP:" + str(version1) + " & BIP:"+str(version2), styleN)

elems = []
elems.append(P)
elems.append(table)


pdf.build(elems)
