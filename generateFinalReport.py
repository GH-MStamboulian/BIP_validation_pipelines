"""
Final report generator
"""

import os
import sys
import pandas as pd
import sys

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
from reportlab.lib.pagesizes import letter
from reportlab.platypus import Table
from reportlab.platypus import TableStyle
from reportlab.lib import colors
from reportlab.lib.units import mm, inch
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle


version1 = sys.argv[1]
version2 = sys.argv[2]

out_dir = sys.argv[3]

data_dir = sys.argv[4]

fileName = out_dir + 'BIP_validation_report.pdf'

accuracy_folder = "accuracy_comparison/"
precision_folder = "precision_comparison/"
sensitivity_folder = "sensitivity_comparison/"
specificity_folder = "specificity_comparison/"

def round_col(col):
    return [round(float(item), 2) for item in list(col)]


def construct_table(data):
    
    table = Table(data)

    # add style
    style = TableStyle([
        #('BACKGROUND', (0,0), (3,0), colors.green),
        ('TEXTCOLOR',(0,0),(-1,0),colors.black),

        ('ALIGN',(0,0),(-1,-1),'CENTER'),
        ('FONTNAME', (0,0), (-1, 0), 'Courier-Bold'),
        ('FONTSIZE', (0,0), (-1, 0), 6),
        ('FONTNAME', (0,1), (-1,-1), 'Courier'),
        ('FONTSIZE', (0,1), (-1,-1), 6),
        #('FONTSIZE', (0,4), (-1,0), 6),
        ('BOTTOMPADDING', (0,0), (-1,0), 8),

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

def construct_table_2(data):
    
    table = Table(data, hAlign='LEFT')

    # add style
    style = TableStyle([
        #('BACKGROUND', (0,0), (3,0), colors.green),
        ('TEXTCOLOR',(0,0),(-1,0),colors.black),

        ('ALIGN',(0,0),(-1,-1),'CENTER'),
        ('FONTNAME', (0,0), (-1, 0), 'Courier-Bold'),
        ('FONTSIZE', (0,0), (-1, 0), 6),
        ('FONTNAME', (0,1), (-1,-1), 'Courier'),
        ('FONTSIZE', (0,1), (-1,-1), 6),
        #('FONTSIZE', (0,4), (-1,0), 6),
        ('BOTTOMPADDING', (0,0), (-1,0), 8),

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

height, width = A4

styles = getSampleStyleSheet()
heading1_styleN = styles['Heading1']

style_title = ParagraphStyle(
        name='Normal',
        fontSize=24,
        leftIndent = 100,
    )

style_header = ParagraphStyle(
        name='Normal',
        fontSize=18,        
    )

style_header_3 = ParagraphStyle(
        name='Normal',
        fontSize=14,        
    )

style_header_4 = ParagraphStyle(
        name='Normal',
        fontSize=12.5,        
    )

style_text = ParagraphStyle(
        name='Normal',
        fontSize=12,        
    )

style_text_1 = ParagraphStyle(
        name='Normal',
        fontSize=10,        
    )

style_footer = ParagraphStyle(
        name='Normal',
        fontSize=6,        
    )

pagesize = (height, width)
elems = []

pdf = SimpleDocTemplate(
    fileName,
    pagesize=pagesize
)

P = Paragraph("BIP Validation Pipelines",  style = style_title)
elems.append(P)

P = Paragraph("The purpose of this report is to summarize results obtained from the analysis of BIP performance metrics, and comparing it across the different"
              "versions of BIP releases. In particular we are assessing four performance measures for each of the two BIP versions that are being compared. "
              "The measures assessed in this analysis are: Accuracy, Precision, Sensitivity and Specificity.", 
              style = style_text)
S2 = Spacer(1, 1.25*inch)
S1 = Spacer(1, 0.7*inch)
S0_5 = Spacer(1, 0.4*inch)
S0 = Spacer(1, 0.2*inch)
elems.append(S2)
elems.append(P)
elems.append(S1)
P = Paragraph("Table Of Contents", style = style_header)
elems.append(P)
elems.append(S0)
P = Paragraph('Accuracy', style = style_text, bulletText='-')
elems.append(P)
P = Paragraph('Precision', style = style_text, bulletText='-')
elems.append(P)
P = Paragraph('Sensitivity', style = style_text, bulletText='-')
elems.append(P)
P = Paragraph('Specificity', style = style_text, bulletText='-')
elems.append(P)
P = Paragraph('Appendix', style = style_text, bulletText='-')
elems.append(P)
elems.append(PageBreak())

############
##  accuracy
############
P = Paragraph("I. ACCURACY", style = style_header)
elems.append(P)
elems.append(S1)
P = Paragraph("In this section we are assessing the analytical accuracy of Guardant360 CDx (G360 CDx). To quantify accuracy for G360, we compare results obtained"
              "to an external orthognal validator conducted by MD Anderson Cancer Center, where we treat the latter as the 'ground truth', and thus allowing us to"
              "quantify and compare the analytical accuracy across the entire pane and the reportable range. For the purpose of this study three sample collections were"
              "made. Aliquots from these samples were taken and tested once by Guardant360 CDx and another by LBP70 (MD Anderson). Accuracy in this context was quantified"
              "as Positive Percent Agreement (PPA) and Negative Percent Agreement (NPA), between the two methods. PPA here is defined as the proportion of variants detected"
              "by G360 CDx, in comparison to the total number of variants in the samples (presented by LBP70). Similarly NPA is defined as the proportion of the variants "
              "not detected by G360 CDx in comparison to the total number of variants not detected by LBP70. We quantify and compare PPA and NPA over all variants between the"
              "two bip versions in comparison."
              , 
              style = style_text)

elems.append(P)



accuracy_dir = out_dir + accuracy_folder

accuracy_PPA_table_f = accuracy_dir + 'accuracy_PPA_merged_table_1_bip_'+version1+'_2_bip_'+version2+'.tsv'
accuracy_PPA_table_df = pd.read_csv(accuracy_PPA_table_f, sep = '\t')

accuracy_PPA_table_df.columns = ['variant \ncategory', 'variant\n type', 'acceptable\n PPA',
       'bip:3.5.2\nCDx+/LBP70+', 'bip:3.5.3\nCDx+/LBP70+', 'bip:3.5.2\nCDx-/LBP70+',
       'bip:3.5.3\nCDx-/LBP70+', 'bip:3.5.2\n PPA', 'bip:3.5.3\n PPA', 'LLCI',
       'ULCI']

data = list()
data.append(list(accuracy_PPA_table_df.columns))
for i, row in accuracy_PPA_table_df.iterrows():
    data.append(list(accuracy_PPA_table_df.loc[i]))
    
table = construct_table(data)
#elems.append(S0_5)
P = Paragraph("Accuracy PPA", style_header_3)
elems.append(S0)
elems.append(P)
elems.append(S0)
elems.append(table)

P = Paragraph("acceptable PPA:minimum acceptance criteria for PPA calculated in accordance to the report 'D-000111'; CDx / LBP70:variants detected by G360 method and MD Anderson"
              "+ for detected, - for not detected; PPA:Positive Percent Agreement; LLCI:Lower level of 95% confidence interval; ULCI:Upper level of 95% confidence intervel"
              ,
              style = style_footer)
elems.append(P)

accuracy_NPA_table_f = accuracy_dir + 'accuracy_NPA_merged_table_1_bip_'+version1+'_2_bip_'+version2+'.tsv'
accuracy_NPA_table_df = pd.read_csv(accuracy_NPA_table_f, sep = '\t')

accuracy_NPA_table_df.columns = ['variant \ncategory', 'variant \ntype', 'acceptable \nNPA',
       'bip:'+version1+'\nCDx+/LBP70-', 'bip:'+version2+'\nCDx+/LBP70-', 'bip:'+version1+'\nCDx-/LBP70-',
       'bip:'+version2+'\nCDx-/LBP70-', 'bip:'+version1+' \nNPA', 'bip:'+version2+' \nNPA', 'LLCI',
       'ULCI']



data = list()
data.append(list(accuracy_NPA_table_df.columns))
for i, row in accuracy_NPA_table_df.iterrows():
    data.append(list(accuracy_NPA_table_df.loc[i]))
    
table = construct_table(data)
#elems.append(S0_5)
P = Paragraph("Accuracy NPA", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)
P = Paragraph("acceptable NPA:minimum acceptance criteria for NPA calculated in accordance to the report 'D-000111'; NPA:Positive Percent Agreement; "
              ,
              style = style_footer)
elems.append(P)


P = Paragraph("Accuracy comparison conclusions", style_header_3)
elems.append(S0)
elems.append(P)
P = Paragraph("Accuracy PPA", style_header_4)
elems.append(S0)
elems.append(P)

accuracy_ppa_conclusions = ""

variants = ["SNV", "Indel", "CNA", "Fusion"]
variant_types = ["Clinically Relevant", "Panel-wide"]

i=0
for variant_type in variant_types:
    for variant in variants:        
        version_1 = accuracy_PPA_table_df[(accuracy_PPA_table_df['variant \ncategory'] == variant) & (accuracy_PPA_table_df['variant\n type'] == variant_type)]
        version_2 = accuracy_PPA_table_df[(accuracy_PPA_table_df['variant \ncategory'] == variant) & (accuracy_PPA_table_df['variant\n type'] == variant_type)]
        if len(version_1['bip:'+version1+'\nCDx-/LBP70+']) == 0:
            missed_by_version1 = 0
        else:
            missed_by_version1 = int(version_1['bip:'+version1+'\nCDx-/LBP70+'])
        if len(version_2['bip:'+version2+'\nCDx-/LBP70+']) == 0:
            missed_by_version2 = 0
        else:
            missed_by_version2 = int(version_2['bip:'+version2+'\nCDx-/LBP70+'])
        
        accuracy_ppa_conclusions += str(i+1) + ")There were "+str(missed_by_version1) +" " + variant_type + " " +variant + " missed by BIP version " + version1 + " and " +str(missed_by_version2) +" " + variant_type + " " + variant + " missed by BIP version " + version2 +"<br/><br/>"
        i+=1 
        
P = Paragraph(accuracy_ppa_conclusions, style = style_text_1)
elems.append(P)


elems.append(S0)
P = Paragraph("Accuracy NPA", style_header_4)
elems.append(S0)
elems.append(P)

accuracy_npa_conclusions = ""

variants = ["SNV", "Indel", "CNA", "Fusion"]
variant_types = ["Clinically Relevant", "Panel-wide"]

i = 0
for variant_type in variant_types:
    for variant in variants:        
        version_1 = accuracy_NPA_table_df[(accuracy_NPA_table_df['variant \ncategory'] == variant) & (accuracy_NPA_table_df['variant \ntype'] == variant_type)]
        version_2 = accuracy_NPA_table_df[(accuracy_NPA_table_df['variant \ncategory'] == variant) & (accuracy_NPA_table_df['variant \ntype'] == variant_type)]        
        if len(version_1['bip:'+version1+'\nCDx+/LBP70-']) == 0:
            missed_by_version1 = 0
        else:
            missed_by_version1 = int(version_1['bip:'+version1+'\nCDx+/LBP70-'])
        if len(version_2['bip:'+version2+'\nCDx+/LBP70-']) == 0:
            missed_by_version2 = 0
        else:
            missed_by_version2 = int(version_2['bip:'+version2+'\nCDx+/LBP70-'])
        
        accuracy_npa_conclusions += str(i+1) + ")There were "+str(missed_by_version1) +" " + variant_type + " " +variant + " false +ve calls by BIP version " + version1 + " and " +str(missed_by_version2) +" " + variant_type + " " + variant + " false +ve calls by BIP version " + version2 +"<br/><br/>"
        i+=1
P = Paragraph(accuracy_npa_conclusions, style = style_text_1)
elems.append(P)
        


elems.append(PageBreak())

#############
##  precision
#############

P = Paragraph("II. PRECISION", style = style_header)
elems.append(P)
elems.append(S1)
P = Paragraph("In this section we are assessing the analytical precision of Guardant360 CDx (G360 CDx). Precision of an assay is defined as the closeness of agreement between"
              "measured qualitative output or quantitative estimates obtained in replicate testing. To assess precision under repeatibilty and within site conditions"
              "experiments were performed over multiple replicates of cell free DNA (cfDNA) sample pools within the same run (repeatability) and different runs (days, instruments)"
              "etc, assessing within site precision. All cdDNA samples contained a set of targetted variants to assess the analytical precision (7 SNVs, 5 Indels, 2 CNAs and 4 Fusions)"
              "Three distinct sample pools were prepared using cfDNA each containing the set of 18 targetted variants. Each of these pools were sequenced by an orthogonal"
              "validation method perfromed by MD Andersion Cancer Center (LBP70 assay), alongside with Guardant Health G360 CDx sequencing. Class level Positive Percentage Agreement (PPA)"
              "for each of the 4 variant classes were quantified accross all conditions. Sample level Negative Percentage Agreement (NPA, False Positives) was quantified accross 90 replicates,"
              "which were asessed over previously defined 32 variants. Finally within-condtion PPA accross all targeted variants were also quantified. There were a total of 3 conditions"
              "that were assessed in this experiment.", 
              style = style_text)

elems.append(P)

precision_dir = out_dir + precision_folder

precision_variantLevelPPA_table_f = precision_dir + 'variantLevelPPA_precision1_bip_' + version1 + '_precision2_bip_' + version2 + '.tsv'
precision_variantLevelPPA_table_df = pd.read_csv(precision_variantLevelPPA_table_f, sep = ',')
del precision_variantLevelPPA_table_df['Unnamed: 0']

precision_variantLevelPPA_table_df.columns = ['variant\nclass', 'bip:'+version1+'\n+ve calls', 'bip:'+version2+'\n+ve calls',
       'expected\ncalls', 'bip:'+version1+'\nPPA', 'bip:'+version2+'\nPPA',
       'bip:'+version1+'\nCI 0.95 lower', 'bip:'+version1+'\nCI 0.95 upper',
       'bip:'+version2+'\nCI 0.95 lower', 'bip:'+version2+'\nCI 0.95 upper']

precision_variantLevelPPA_table_df['bip:'+version1+'\nPPA'] = round_col(precision_variantLevelPPA_table_df['bip:'+version1+'\nPPA'])
precision_variantLevelPPA_table_df['bip:'+version2+'\nPPA'] = round_col(precision_variantLevelPPA_table_df['bip:'+version2+'\nPPA'])
precision_variantLevelPPA_table_df['bip:'+version1+'\nCI 0.95 lower'] = round_col(precision_variantLevelPPA_table_df['bip:'+version1+'\nCI 0.95 lower'])
precision_variantLevelPPA_table_df['bip:'+version1+'\nCI 0.95 upper'] = round_col(precision_variantLevelPPA_table_df['bip:'+version1+'\nCI 0.95 upper'])
precision_variantLevelPPA_table_df['bip:'+version2+'\nCI 0.95 lower'] = round_col(precision_variantLevelPPA_table_df['bip:'+version2+'\nCI 0.95 lower'])
precision_variantLevelPPA_table_df['bip:'+version2+'\nCI 0.95 upper'] = round_col(precision_variantLevelPPA_table_df['bip:'+version2+'\nCI 0.95 upper'])

data = list()
data.append(list(precision_variantLevelPPA_table_df.columns))
for i, row in precision_variantLevelPPA_table_df.iterrows():
    data.append(list(precision_variantLevelPPA_table_df.loc[i]))
    
table = construct_table(data)
#elems.append(S0_5)
P = Paragraph("Precision PPA (Variant level)", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)
P = Paragraph("+ve calls:number of variants detected by the method; expected calls:number of variants expected to be found in a sample"
              ,
              style = style_footer)
elems.append(P)


precision_variantLevelFalseNegatives_table_f = precision_dir + 'variantLevelFalseNegatives_precision1_bip_'+version1+'_precision2_bip_'+version2+'.tsv'
precision_variantLevelFalseNegatives_table_df = pd.read_csv(precision_variantLevelFalseNegatives_table_f, sep = ',')
del precision_variantLevelFalseNegatives_table_df['Unnamed: 0']
precision_variantLevelFalseNegatives_table_df.columns = ['variant\nclass', 'variant', 'bip:'+version1+'\nFN', 'bip:'+version2+'\nFN']
data = list()
data.append(list(precision_variantLevelFalseNegatives_table_df.columns))
for i, row in precision_variantLevelFalseNegatives_table_df.iterrows():
    data.append(list(precision_variantLevelFalseNegatives_table_df.loc[i]))
    
table = construct_table(data)
#elems.append(S0_5)
P = Paragraph("Variants missed by BIP", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)
P = Paragraph("FN:number of true variants in the sample failed to be classified as positive by the method"
              ,
              style = style_footer)
elems.append(P)
elems.append(PageBreak())

precision_withinCondtionPPA_table_f = precision_dir + 'withinConditionPPA_precision1_bip_' + version1 + '_precision2_bip_' + version2 + '.tsv'
precision_withinCondtionPPA_table_df = pd.read_csv(precision_withinCondtionPPA_table_f, sep = ',')
del precision_withinCondtionPPA_table_df['Unnamed: 0']

precision_withinCondtionPPA_table_df.columns = ['Condition', 'bip:'+version1+'\n+ve calls', 'bip:'+version2+'\n+ve calls',
       'expected\ncalls', 'bip:'+version1+'\nPPA', 'bip:'+version2+'\nPPA',
       'bip:'+version1+'\nCI 0.95 lower', 'bip:'+version1+'\nCI 0.95 upper',
       'bip:'+version2+'\nCI 0.95 lower', 'bip:'+version2+'\nCI 0.95 upper']

precision_withinCondtionPPA_table_df['bip:'+version1+'\nCI 0.95 lower'] = round_col(precision_withinCondtionPPA_table_df['bip:'+version1+'\nCI 0.95 lower'])
precision_withinCondtionPPA_table_df['bip:'+version1+'\nCI 0.95 upper'] = round_col(precision_withinCondtionPPA_table_df['bip:'+version1+'\nCI 0.95 upper'])
precision_withinCondtionPPA_table_df['bip:'+version2+'\nCI 0.95 lower'] = round_col(precision_withinCondtionPPA_table_df['bip:'+version2+'\nCI 0.95 lower'])
precision_withinCondtionPPA_table_df['bip:'+version2+'\nCI 0.95 upper'] = round_col(precision_withinCondtionPPA_table_df['bip:'+version2+'\nCI 0.95 upper'])

data = list()
data.append(list(precision_withinCondtionPPA_table_df.columns))
for i, row in precision_withinCondtionPPA_table_df.iterrows():
    data.append(list(precision_withinCondtionPPA_table_df.loc[i]))
    
table = construct_table(data)
#elems.append(S0_5)
P = Paragraph("Precision PPA (Across different conditions)", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)


precision_sampleLevelNPA_table_f = precision_dir + 'sampleLevelNPA_precision1_bip_' + version1 + '_precision2_bip_' + version2 + '.tsv'
precision_sampleLevelNPA_table_df = pd.read_csv(precision_sampleLevelNPA_table_f, sep = ',')
del precision_sampleLevelNPA_table_df['Unnamed: 0']

precision_sampleLevelNPA_table_df.columns = ['N samples', 'bip:' + version1 + '\nN samples with FP', 'bip:' + version2 + '\nN samples with FP',
       'bip:' + version1 + '\nsample level NPA', 'bip:' + version2 + '\nsample level NPA', '95% CI']

data = list()
data.append(list(precision_sampleLevelNPA_table_df.columns))
for i, row in precision_sampleLevelNPA_table_df.iterrows():
    data.append(list(precision_sampleLevelNPA_table_df.loc[i]))
    
table = construct_table(data)
#elems.append(S0_5)
P = Paragraph("Precision NPA (Across samples)", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)
P = Paragraph("N samples with FP:Sample level False Postives, number of samples failed to be considered as positive"
              ,
              style = style_footer)
elems.append(P)


P = Paragraph("Precision comparison conclusions", style_header_3)
elems.append(S0)
elems.append(P)
P = Paragraph("Precision variant level PPA", style_header_4)
elems.append(S0)
elems.append(P)

precision_ppa_conclusions = ""

variants = ["SNV", "INDEL", "CNV", "FUSION"]

for i, variant in enumerate(variants):
    version_1 = precision_variantLevelPPA_table_df[precision_variantLevelPPA_table_df['variant\nclass'] == variant]
    version_2 = precision_variantLevelPPA_table_df[precision_variantLevelPPA_table_df['variant\nclass'] == variant]
    expected_calls = int(precision_variantLevelPPA_table_df[precision_variantLevelPPA_table_df['variant\nclass'] == variant]['expected\ncalls'])
    if len(version_1['bip:'+version1+'\n+ve calls']) == 0:
        missed_by_version1 = 0
    else:
        missed_by_version1 = int(version_1['bip:'+version1+'\n+ve calls'])
    if len(version_2['bip:'+version2+'\n+ve calls']) == 0:
        missed_by_version2 = 0
    else:
        missed_by_version2 = int(version_2['bip:'+version2+'\n+ve calls'])

    precision_ppa_conclusions += str(i+1) + ")There were "+str(missed_by_version1) +" " + variant + " calls by BIP version " + version1 + " and " +str(missed_by_version2) +" " + variant + " calls by BIP version " + version2 +" from a total of " +str(expected_calls)+ " expected calls<br/><br/>"

P = Paragraph(precision_ppa_conclusions, style = style_text_1)
elems.append(P)

elems.append(S0)
P = Paragraph("Precision cross conditions PPA", style_header_4)
elems.append(S0)
elems.append(P)
precision_ppa_conclusions = ""

variants = ["PC 1", "PC 2", "PC 3"]

for i, variant in enumerate(variants):
    version_1 = precision_withinCondtionPPA_table_df[precision_withinCondtionPPA_table_df['Condition'] == variant]
    version_2 = precision_withinCondtionPPA_table_df[precision_withinCondtionPPA_table_df['Condition'] == variant]
    expected_calls = int(precision_withinCondtionPPA_table_df[precision_withinCondtionPPA_table_df['Condition'] == variant]['expected\ncalls'])
    if len(version_1['bip:'+version1+'\n+ve calls']) == 0:
        missed_by_version1 = 0
    else:
        missed_by_version1 = int(version_1['bip:'+version1+'\n+ve calls'])
    if len(version_2['bip:'+version2+'\n+ve calls']) == 0:
        missed_by_version2 = 0
    else:
        missed_by_version2 = int(version_2['bip:'+version2+'\n+ve calls'])

    precision_ppa_conclusions += str(i+1) + ")There were "+str(missed_by_version1) + " calls by BIP version " + version1 + " and " +str(missed_by_version2) +" calls by BIP version " + version2 +" from a total of " +str(expected_calls)+ " expected calls in condition " +variant+ "<br/><br/>"

P = Paragraph(precision_ppa_conclusions, style = style_text_1)
elems.append(P)



elems.append(PageBreak())

###############
##  sensitivity
###############

P = Paragraph("III. SENSITIVITY", style = style_header)

elems.append(P)
elems.append(S1)
P = Paragraph("In this section we are assessing the analytical Sensitivity, which in other words is defined as the Limit of Detection (LoD) for Guardant360 CDx (G360 CDx)."
              "Sensitivity in this context is defined as the concentration of the analyte (number of mutant molecules) in the sample, which could be measured as the Mutant"
              "Allele Frequency (MAF), such that the variant of interest would be detected by BIP in 95% or more of the experiments. Sample pools were created by extracting"
              "cell free DNA identified for each of the individual targetted variants. These pools are mixed together as separate variant pools for each of the individual"
              "variants for the purpose of increasing the sample masses (mutant molecules). These variant pools then were mixed with wild type pools with different volumens"
              "to generate the targetted dilution levels to test the LoD. A total of 168 clinical samples were used to create these pools. Two sample pools were maintained"
              "throughout the experiment, one starting at 5ng and the other at 30ng concentrations, 5 dilution levels were created to test for LoD. Regression models were used"
              "to estimate the targetted MAF values for detecting variants at least 95% of the time if there were 3 or more detection points available, otherwise the lowest MAF "
              "estimate where detection is at least 95%, was used as the target MAF value. LoD was quantified over a total of 36 variants (17 unique variants once starting from"
              "5ng samples and another time from 30ng). "              
              , 
              style = style_text)

elems.append(P)

sensitivity_dir = out_dir + sensitivity_folder

sensitivity_variantLoD_table_f = sensitivity_dir + 'LoD1_bip_' + version1 + '_LoD2_bip_' + version2 + '.tsv'
sensitivity_variantLoD_table_df = pd.read_csv(sensitivity_variantLoD_table_f, sep = '\t')

sensitivity_variantLoD_table_df.columns = ['Variant', 'Variant\nType', 'LoD 5ng\nbip:' + version1, 'LoD 5ng\nbip:' + version2,
       'LoD 30ng\nbip:' + version1, 'LoD 30ng\nbip:' + version2]


data = list()
data.append(list(sensitivity_variantLoD_table_df.columns))
for i, row in sensitivity_variantLoD_table_df.iterrows():
    data.append(list(sensitivity_variantLoD_table_df.loc[i]))
    
table = construct_table(data)
#elems.append(S0_5)
P = Paragraph("Sensitivity (Limit of Detection)", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)
P = Paragraph("LoD 5ng:Limit of Detection in 5ng concetration samples; LoD 30ng:Limit of Detection in 30ng concentration samples"
              ,
              style = style_footer)
elems.append(P)


P = Paragraph("Sensitivity comparison conclusions", style_header_3)
elems.append(S0)
elems.append(P)
P = Paragraph("Limit of Detection", style_header_4)
elems.append(S0)
elems.append(P)

LoD_conclusions = ""


for i, lod_row in sensitivity_variantLoD_table_df.iterrows():
    variant = lod_row['Variant']
    variant_type = lod_row['Variant\nType']
    new_line = False
    if lod_row['LoD 5ng\nbip:' + version1] != lod_row['LoD 5ng\nbip:' + version2]:
        LoD_conclusions += str(i+1)+ ") The LoD values for the " + variant_type + " " + variant + " at 5ng differs between bip:" + version1 + " LoD:" + str(lod_row['LoD 5ng\nbip:' + version1]) + " and bip:" + version2 + " LoD:" + str(lod_row['LoD 5ng\nbip:' + version2])
        new_line = True
    if lod_row['LoD 30ng\nbip:' + version1] != lod_row['LoD 30ng\nbip:' + version2]:    
        if not new_line:
            LoD_conclusions +=str(i+1)+ ") "
        LoD_conclusions += "The LoD values for the " + variant_type + " " + variant + " at 30ng differs between bip:" + version1 + " LoD:" + str(lod_row['LoD 30ng\nbip:' + version1]) + " and bip:" + version2 + " LoD:" + str(lod_row['LoD 30ng\nbip:' + version2])
        new_line = True
    if new_line:
        LoD_conclusions += "<br/>"
    
if LoD_conclusions == "":
    LoD_conclusions = "LoD values do not differ between the two bip versions in comparsion (bip:"+ str(version1)+ " & bip:" + str(version2) +") accross the targeted " + str(len(sensitivity_variantLoD_table_df)) + " variants"
        
P = Paragraph(LoD_conclusions, style = style_text_1)
elems.append(P)



elems.append(PageBreak())



###############
##  specificity
###############

P = Paragraph("IV. SPECIFICITY", style = style_header)

elems.append(P)
elems.append(S1)

P = Paragraph("In this section we are assessing the analytical Specificity, which in other words is defined as the Limit of Blank (LoB) for Guardant360 CDx (G360 CDx). "              
              "Specificity in context of Guardant360 CDx is defined as the highest measurement result that is likely to be observed in a blank sample (i.e. sample not "
              "containing the variant). Limit of Blank here is quantified as the fraction of pisitively identified variants in blank samples (analyticial False Positives)."
              "For the purpose of this experiment  cfDNA from 20 healthy donoros were extracted (with no reported cancer cases). Two orthogonal methods were used as validation"
              "over the 20 samples used, an external MD Anderson NGS-based comparsion (LBP70), as well as guardant health sequencing was performed over these samples. All variants"
              "detected by BIP (any version) were first filtered for germline variants. At any point a variant reported by BIP was considered to be a false positive if it is "
              "not verified by any of the orthogonal tests used, and hence analytical false positive is quantified. Variants verified by any of the orthogonal methods were ."
              "instead considered as analytical true positive. Analytical False Positive variants were reported for the BIP versions in comparison."
              , 
              style = style_text)

elems.append(P)

specificity_dir = out_dir + specificity_folder

specificity_variantLoB_table_f = specificity_dir + 'LoB1_bip_' + version1 + '_LoB2_bip_' + version2 + '.tsv'
specificity_variantLoB_table_df = pd.read_csv(specificity_variantLoB_table_f, sep = '\t')

specificity_variantLoB_table_df.columns = ['variant\ntype', 'chrom', 'position', 'gene', 'transcript\nid', 'exon',
       'mut\nnt', 'mut\naa', 'mut\nkey', 'mut\ncnt', 'pool.maf', 'cfDNA.maf',
       'gDNA.maf', 'detected\nbip:' + version1, 'detected\nbip:' + version2]

#del specificity_variantLoB_table_df['transcript\nid']
del specificity_variantLoB_table_df['cfDNA.maf']
del specificity_variantLoB_table_df['gDNA.maf']

data = list()
data.append(list(specificity_variantLoB_table_df.columns))
for i, row in specificity_variantLoB_table_df.iterrows():
    data.append(list(specificity_variantLoB_table_df.loc[i]))
    
table = construct_table(data)
#elems.append(S0_5)
P = Paragraph("Specificity (Limit of Blank, False Positives)", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)
P = Paragraph("Chrom:the chromosome number for which the variant is found; position:starting coordinate of the variant on the chromosome; transcript id:the id of the transcript"
              "where the variant is on; exon:the exon number where the variant is located; mut nt: nucelotide level mutation; mut aa:Amino acid level mutation; mut cnt: "
              "number of mutant molecules in the sample; pool,maf:mutant allele frequency in the pooled sample; not detected:variant detected or not within that sample."
              ,
              style = style_footer)
elems.append(P)


P = Paragraph("Specificity comparison conclusions", style_header_3)
elems.append(S0)
elems.append(P)
P = Paragraph("Limit of Blank", style_header_4)
elems.append(S0)
elems.append(P)

LoB_conclusions = ""

version_1_variant_lob_dic, version_2_variant_lob_dic = {'CNV':0, 'FUSION':0, 'INDEL':0, 'SNV':0}, {'CNV':0, 'FUSION':0, 'INDEL':0, 'SNV':0}

for i, lob_row in specificity_variantLoB_table_df.iterrows():
    if str(lob_row['gene']) != 'nan':        
        variant_type = lob_row['variant\ntype']
        version_1_detected = lob_row['detected\nbip:' + version1]
        version_2_detected = lob_row['detected\nbip:' + version2]
        if version_1_detected:
            version_1_variant_lob_dic[variant_type] += 1
        if version_2_detected:
            version_2_variant_lob_dic[variant_type] += 1    
for i, k in enumerate(version_1_variant_lob_dic):
    LoB_conclusions += str(i+1) + ") there were a total of "+str(version_1_variant_lob_dic[k]) + " " + k + " False detection cases for bip:" + version1 +", and a total of "+str(version_2_variant_lob_dic[k]) + " " + k + " False detection cases for bip:" + version2 + "<br/>" 

P = Paragraph(LoB_conclusions, style = style_text_1)
elems.append(P)


elems.append(PageBreak())

############
##  appendix
############

specificity_samples_f = data_dir + "specificity/" + version1+"/G360/sample_list.csv"
specificity_flowcells = list(set(list(pd.read_csv(specificity_samples_f, sep = ',')['runid'])))

sensitivity_samples_f = data_dir + "sensitivity/" + version1+"/sample_list.181204.csv"
sensitivity_flowcells = list(set(list(pd.read_csv(sensitivity_samples_f, sep = ',')['runid'])))

precision_samples_f = data_dir + "precision/" + version1+"/LineData_variant.csv"
precision_flowcells = list(set(list(pd.read_csv(precision_samples_f, sep = ',')['runid'])))

accuracy_samples_f = data_dir + "accuracy/" + version1+"/line_data.xlsx"
accuracy_samples_df = pd.read_excel(accuracy_samples_f, sheet_name = "Sample Manifest")
accuracy_flowcells = list(set(list(accuracy_samples_df[accuracy_samples_df['Processing Lab'] == 'cdx']['Flowcell ID'])))

P = Paragraph("IV. APPENDIX", style = style_header)
elems.append(P)

data = list()
data.append(["Specificity Flowcells"])
data.extend([[item] for item in specificity_flowcells])
    
table = construct_table_2(data)

P = Paragraph("Specificity Flowcells list", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)


data = list()
data.append(["Sensitivity Flowcells"])
data.extend([[item] for item in sensitivity_flowcells])
    
table = construct_table_2(data)
elems.append(S0)
P = Paragraph("Sensitivity Flowcells list", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)


data = list()
data.append(["Precision Flowcells"])
data.extend([[item] for item in precision_flowcells])
    
table = construct_table_2(data)
elems.append(S0)
P = Paragraph("Precision Flowcells list", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)


data = list()
data.append(["Accuracy Flowcells"])
data.extend([[item] for item in accuracy_flowcells])
    
table = construct_table_2(data)
elems.append(S0)
P = Paragraph("accuracy Flowcells list", style_header_3)
elems.append(S0_5)
elems.append(P)
elems.append(S0)
elems.append(table)


pdf.build(elems)
