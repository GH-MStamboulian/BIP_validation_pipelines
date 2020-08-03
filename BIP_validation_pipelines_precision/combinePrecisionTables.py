"""
script that will take the Positive Percent Agreement from both bip versions and will merge and combine them to one table
"""

import pandas as pd
import sys

out_1_dir = sys.argv[1]
bip_version_1 = sys.argv[2]
out_2_dir = sys.argv[3]
bip_version_2 = sys.argv[4]
out_dir = sys.argv[5]



table_1_f = out_1_dir + "Table7_variant_class_level_PPA.tsv"
table_2_f = out_2_dir + "Table7_variant_class_level_PPA.tsv"

table_3_f = out_1_dir + "Table7_variant_class_level_FalseNegatives_PPA.tsv"
table_4_f = out_2_dir + "Table7_variant_class_level_FalseNegatives_PPA.tsv"

table_5_f = out_1_dir + "Table12_within_condition_PPA.tsv"
table_6_f = out_2_dir + "Table12_within_condition_PPA.tsv"

table_7_f = out_1_dir + "Table10_sample_level_NPA.tsv"
table_8_f = out_2_dir + "Table10_sample_level_NPA.tsv"

fileName = out_dir + 'precision1_bip_' + bip_version_1 + '_precision2_bip_' + bip_version_2 + '.pdf'

"""
table_1_f = sys.argv[1]
bip_version_1 = sys.argv[2]
table_2_f = sys.argv[3]
bip_version_2 = sys.argv[4]
out_dir = sys.argv[5]
"""

############################
##  variant class level PPAs
############################

table_1_df = pd.read_csv(table_1_f, sep = '\t')
table_2_df = pd.read_csv(table_2_f, sep = '\t')


table_1_df.columns = ["variant_class", "bip:"+bip_version_1+"_+ve_calls", "expected_calls", "bip:"+bip_version_1+"_PPA", "bip:"+bip_version_1+"_CI_0.95_lower", "bip:"+bip_version_1+"_CI_0.95_upper"]
table_2_df.columns = ["variant_class", "bip:"+bip_version_2+"_+ve_calls", "expected_calls", "bip:"+bip_version_2+"_PPA", "bip:"+bip_version_2+"_CI_0.95_lower", "bip:"+bip_version_2+"_CI_0.95_upper"]

table_variant_level_ppa_df = pd.DataFrame(columns = ["variant_class", "bip:"+bip_version_1+"_+ve_calls", "bip:"+bip_version_2+"_+ve_calls", "expected_calls", "bip:"+bip_version_1+"_PPA", "bip:"+bip_version_2+"_PPA", "bip:"+bip_version_1+"_CI_0.95_lower", "bip:"+bip_version_1+"_CI_0.95_upper", "bip:"+bip_version_2+"_CI_0.95_lower", "bip:"+bip_version_2+"_CI_0.95_upper"])
table_variant_level_ppa_df["variant_class"] = table_1_df["variant_class"]
table_variant_level_ppa_df["bip:"+bip_version_1+"_+ve_calls"] = table_1_df["bip:"+bip_version_1+"_+ve_calls"]
table_variant_level_ppa_df["bip:"+bip_version_2+"_+ve_calls"] = table_2_df["bip:"+bip_version_2+"_+ve_calls"]
table_variant_level_ppa_df["expected_calls"] = table_1_df["expected_calls"]
table_variant_level_ppa_df["bip:"+bip_version_1+"_PPA"] = table_1_df["bip:"+bip_version_1+"_PPA"]
table_variant_level_ppa_df["bip:"+bip_version_2+"_PPA"] = table_2_df["bip:"+bip_version_2+"_PPA"]
table_variant_level_ppa_df["bip:"+bip_version_1+"_CI_0.95_lower"] = table_1_df["bip:"+bip_version_1+"_CI_0.95_lower"]
table_variant_level_ppa_df["bip:"+bip_version_1+"_CI_0.95_upper"] = table_1_df["bip:"+bip_version_1+"_CI_0.95_upper"]
table_variant_level_ppa_df["bip:"+bip_version_2+"_CI_0.95_lower"] = table_2_df["bip:"+bip_version_2+"_CI_0.95_lower"]
table_variant_level_ppa_df["bip:"+bip_version_2+"_CI_0.95_upper"] = table_2_df["bip:"+bip_version_2+"_CI_0.95_upper"]

table_variant_level_ppa_df.to_csv(out_dir + 'variantLevelPPA_precision1_bip_' + bip_version_1 + '_precision2_bip_' + bip_version_2 + '.pdf')

#######################################
##  variant class level False Negatives
#######################################

table_3_df = pd.read_csv(table_3_f, sep = '\t')
table_4_df = pd.read_csv(table_4_f, sep = '\t')

table_variant_level_False_Negatives_df = pd.DataFrame(columns = ["variant_class", "variant", "bip:"+bip_version_1+"_FN", "bip:"+bip_version_2+"_FN"])

all_variant_class_variant_tuples = list()
for i, row in table_3_df.iterrows():
    all_variant_class_variant_tuples.append((row["variant_class"], row["variant"]))
    
for i, row in table_4_df.iterrows():
    tup = (row["variant_class"], row["variant"])
    if tup not in all_variant_class_variant_tuples:
        all_variant_class_variant_tuples.append(tup)
        
table_variant_level_False_Negatives_df["variant_class"] = [item[0] for item in all_variant_class_variant_tuples]       
table_variant_level_False_Negatives_df["variant"] = [item[1] for item in all_variant_class_variant_tuples]

table_variant_level_False_Negatives_df["bip:"+bip_version_1+"_FN"] = [int(table_3_df[table_3_df["variant"] == var]["number of false negatives"]) if len(table_3_df[table_3_df["variant"] == var]) > 0 else "NaN" for var in  [item[1] for item in all_variant_class_variant_tuples]]
table_variant_level_False_Negatives_df["bip:"+bip_version_2+"_FN"] = [int(table_3_df[table_3_df["variant"] == var]["number of false negatives"]) if len(table_4_df[table_3_df["variant"] == var]) > 0 else "NaN" for var in  [item[1] for item in all_variant_class_variant_tuples]]

table_variant_level_False_Negatives_df.to_csv(out_dir + 'variantLevelFalseNegatives_precision1_bip_' + bip_version_1 + '_precision2_bip_' + bip_version_2 + '.pdf')

########################
##  within condition PPA
########################

table_5_df = pd.read_csv(table_5_f, sep = '\t')
table_6_df = pd.read_csv(table_6_f, sep = '\t')

table_within_condition_PPA_df = pd.DataFrame(columns = ["Condition", "bip:" + bip_version_1 + "_n_+ve_calls", "bip:" + bip_version_2 + "_n_+ve_calls", "n_expected_calls", "bip:"+bip_version_1+"_PPA", "bip:"+bip_version_1+"_PPA", "bip:"+bip_version_1+"_CI_0.95_lower", "bip:"+bip_version_1+"_CI_0.95_upper", "bip:"+bip_version_2+"_CI_0.95_lower", "bip:"+bip_version_2+"_CI_0.95_upper"])

table_within_condition_PPA_df["Condition"] = table_5_df["Condition"]
table_within_condition_PPA_df["bip:" + bip_version_1 + "_n_+ve_calls"] = table_5_df["number of positive calls"]
table_within_condition_PPA_df["bip:" + bip_version_2 + "_n_+ve_calls"] = table_6_df["number of positive calls"]
table_within_condition_PPA_df["n_expected_calls"] = table_5_df["total number of expected calls"]
table_within_condition_PPA_df["bip:" + bip_version_1 + "_PPA"] = table_5_df["number of positive calls"]
table_within_condition_PPA_df["bip:" + bip_version_2 + "_PPA"] = table_6_df["number of positive calls"]
table_within_condition_PPA_df["bip:"+bip_version_1+"_CI_0.95_lower"] = table_5_df["CI_0.95_lower"]
table_within_condition_PPA_df["bip:"+bip_version_1+"_CI_0.95_upper"] = table_5_df["CI_0.95_upper"]
table_within_condition_PPA_df["bip:"+bip_version_2+"_CI_0.95_lower"] = table_6_df["CI_0.95_lower"]
table_within_condition_PPA_df["bip:"+bip_version_2+"_CI_0.95_upper"] = table_6_df["CI_0.95_upper"]

table_within_condition_PPA_df.to_csv(out_dir + 'withinConditionPPA_precision1_bip_' + bip_version_1 + '_precision2_bip_' + bip_version_2 + '.pdf')

####################
##  sample level NPA
####################

table_7_df = pd.read_csv(table_7_f, sep = '\t')
table_8_df = pd.read_csv(table_8_f, sep = '\t')

table_sample_level_NPA_df = pd.DataFrame(columns = ["N_samples", "bip:"+bip_version_1+ "N_samples_with_FP", "bip:"+bip_version_2+ "N_samples_with_FP", "bip:"+bip_version_1+ "_sample_level_NPA", "bip:"+bip_version_2+ "_sample_level_NPA", "95%_CI"])

table_sample_level_NPA_df["N_samples"] = table_7_df["Number_samples"]
table_sample_level_NPA_df["bip:"+bip_version_1+ "N_samples_with_FP"] = table_7_df["N_samples_with_False_Positives"]
table_sample_level_NPA_df["bip:"+bip_version_2+ "N_samples_with_FP"] = table_8_df["N_samples_with_False_Positives"]
table_sample_level_NPA_df["bip:"+bip_version_1+ "_sample_level_NPA"] = table_7_df["sample_level_NPA"]
table_sample_level_NPA_df["bip:"+bip_version_2+ "_sample_level_NPA"] = table_8_df["sample_level_NPA"]
table_sample_level_NPA_df["95%_CI"] = table_7_df["95%_CI"]

table_sample_level_NPA_df.to_csv(out_dir + 'sampleLevelNPA_precision1_bip_' + bip_version_1 + '_precision2_bip_' + bip_version_2 + '.pdf')

################
##  build report
################

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.pagesizes import letter
from reportlab.platypus import Table
from reportlab.platypus import TableStyle
from reportlab.lib import colors
from reportlab.lib.units import mm, inch
from reportlab.lib.styles import getSampleStyleSheet



styles = getSampleStyleSheet()
styleN = styles['Heading1']

pagesize = (600 * mm, 400 * mm)

pdf = SimpleDocTemplate(
    fileName,
    pagesize=pagesize
)



def construct_table(data):
    
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
    
    return table

data = list()
data.append(list(table_variant_level_ppa_df.columns))

for i, row in table_variant_level_ppa_df.iterrows():
    data.append(list(table_variant_level_ppa_df.loc[i]))
    
table = construct_table(data)

P = Paragraph("Table summarizing the comparison of precision (Positive Percent Agreement) across variant classes between two BIP versions", styleN)

elems = []
elems.append(P)
elems.append(table)


data = list()
data.append(list(table_variant_level_False_Negatives_df.columns))

for i, row in table_variant_level_False_Negatives_df.iterrows():
    data.append(list(table_variant_level_False_Negatives_df.loc[i]))
    
table = construct_table(data)

P = Paragraph("Table summarizing the comparison of varriant\n class level False Negatives across variant\n classes between two BIP versions", styleN)

P2 = Spacer(1, 1.25*inch)

elems.append(P2)
elems.append(P)
elems.append(table)


data = list()
data.append(list(table_within_condtion_PPA_df.columns))

for i, row in table_within_condtion_PPA_df.iterrows():
    data.append(list(table_within_condtion_PPA_df.loc[i]))
    
table = construct_table(data)

P = Paragraph("Table summarizing the comparison of within experimental condtion precision between two BIP versions", styleN)

P2 = Spacer(1, 1.25*inch)

elems.append(P2)
elems.append(P)
elems.append(table)


data = list()
data.append(list(table_sample_level_NPA.columns))

for i, row in table_sample_level_NPA.iterrows():
    data.append(list(table_sample_level_NPA.loc[i]))
    
table = construct_table(data)

P = Paragraph("Table summarizing the comparison of number of samples with False Postives between two BIP versions", styleN)

P2 = Spacer(1, 1.25*inch)

elems.append(P2)
elems.append(P)
elems.append(table)


pdf.build(elems)
