COLUMN_KEYS = {
    'rm_reportable': 'True if a SNV or Indel is in the RRR',
    'In RRR': 'In Restricted Reportable Range',
    'Cancer Type': 'The diagnosed cancer type',
    'Patient ID': 'Identifier for donor',
    'Sample ID': 'Identifier the sample in each Streck BCT',
    'Flowcell ID': 'Identifier for the flowcell on which the sample was processed',
    'Amino acid change': 'Variant amino acid change',
    'Genomic position': 'Chromosome coordinate',
    'Gene': 'Gene',
    'Sequencing Data Folder': 'Folder containing BIP output',
    'Raw Sequencing Data Folder': 'Folder containing raw sequencing data',
    'Biobank Source': 'The biobank from which the plasma aliquot was selected',
    'Panel': 'The version of the panel with which the sample was sequenced',
    'Processing Lab': 'The lab at which the sample was processed',
    'Sample Collection': 'Either Sample Collection 1, 2, or 3 as desribed in the study protocol',
    'Detection': 'Integer encoding if the variant was not detected (0), detected (1), or ' \
                 'reported (2)',
    'homopolymer': 'A boolean indicating if the variant is adjacent to a 5bp homopolymer stretch',
    'long_indel': 'A boolean indicating if the variant has an insertion or deletion >30bp',
    'exclude_reason': 'The information provided to Guardant about MDL QC failures',
    'exclude_reason_status': 'The PASS/FAIL status for overall QC non-CDx-processed sampless',
    'total': 'A summary PASS/FAIL column for all in-process QC, sequencing and cfDNA' \
             'quantitations',
    'splice_effect': 'Splice Effect',
    'chrom': 'Chromosome',
    'del_size': 'Size of deletion (if any)',
    'panel_GH2.10': 'boolean indicating if this variant is designed to be detected by this panel',
    'panel_GH2.10.1': 'boolean indicating if this variant is designed to be detected by this' \
                      ' panel',
    'panel_GH2.11': 'boolean indicating if this variant is designed to be detected by this panel',
    'panel_MDACC_v1.0': 'boolean indicating if this variant is designed to be detected by this' \
                        ' panel',

    'mdl': 'The original call by LBP70',
    'cdx': 'The original call by G360 CDx',
    'ldt': 'The original call by G360 LDT',
    'mdl_reviewed': 'The reviewed call for LBP70',
    'cdx_reviewed': 'The reviewed call for G360 CDx',
    'ldt_reviewed': 'The reviewed call for G360 LDT',

    'runid': 'Identifier for the flowcell on which the sample was processed',
    'status': 'the PASS/FAIL/REVIEW status of the QC metric',
    'value': 'the measured value for the QC metric',
    'unit': 'the unit of measure for the QC metric',
    'run_sample_id': 'Identifier for each sample',
    'threshold': 'the threshold for the QC metric',
    'metric': 'the name of the QC metric',
    'category': 'the category of the QC metric (sample, control, or flowcell)',
    'verbose_name': 'verbose name for the QC metric',
    'operator': 'the operator to be used in conjuction with the QC threshold',
    'decimal_places': 'the number of decimal places to keep for the QC metric',

    'flowcell_aio_controls': 'A flowcell-level QC',
    'flowcell_auto_qc_total': 'A flowcell-level QC',
    'flowcell_cluster_density': 'A flowcell-level QC',
    'flowcell_clusters_passing_filter': 'A flowcell-level QC',
    'flowcell_phasing_1': 'A flowcell-level QC',
    'flowcell_phasing_2': 'A flowcell-level QC',
    'flowcell_phasing_index': 'A flowcell-level QC',
    'flowcell_prephasing_1': 'A flowcell-level QC',
    'flowcell_prephasing_2': 'A flowcell-level QC',
    'flowcell_prephasing_index': 'A flowcell-level QC',
    'flowcell_qscore_1': 'A flowcell-level QC',
    'flowcell_qscore_2': 'A flowcell-level QC',
    'flowcell_qscore_index': 'A flowcell-level QC',
    'flowcell_aio_controls_status': 'A flowcell-level QC',
    'flowcell_auto_qc_total_status': 'A flowcell-level QC',
    'flowcell_cluster_density_status': 'A flowcell-level QC',
    'flowcell_clusters_passing_filter_status': 'A flowcell-level QC',
    'flowcell_phasing_1_status': 'A flowcell-level QC',
    'flowcell_phasing_2_status': 'A flowcell-level QC',
    'flowcell_phasing_index_status': 'A flowcell-level QC',
    'flowcell_prephasing_1_status': 'A flowcell-level QC',
    'flowcell_prephasing_2_status': 'A flowcell-level QC',
    'flowcell_prephasing_index_status': 'A flowcell-level QC',
    'flowcell_qscore_1_status': 'A flowcell-level QC',
    'flowcell_qscore_2_status': 'A flowcell-level QC',
    'flowcell_qscore_index_status': 'A flowcell-level QC',
    'flowcell_autoqc_total': 'A flowcell-level QC',
    'flowcell_autoqc_total_status': 'A flowcell-level QC',

    'aiocontrol_false_positive_indel': 'A control-level QC',
    'aiocontrol_false_positive_snv': 'A control-level QC',
    'aiocontrol_sensitivity_cnv': 'A control-level QC',
    'aiocontrol_sensitivity_fusion': 'A control-level QC',
    'aiocontrol_sensitivity_indel': 'A control-level QC',
    'aiocontrol_sensitivity_snv': 'A control-level QC',
    'sample_contamination_pct': 'A control-level QC',
    'aiocontrol_false_positive_indel_status': 'A control-level QC',
    'aiocontrol_false_positive_snv_status': 'A control-level QC',
    'aiocontrol_sensitivity_cnv_status': 'A control-level QC',
    'aiocontrol_sensitivity_fusion_status': 'A control-level QC',
    'aiocontrol_sensitivity_indel_status': 'A control-level QC',
    'aiocontrol_sensitivity_snv_status': 'A control-level QC',
    'sample_contamination_pct_status': 'A control-level QC',

    'pathogenic_germline_num': 'A sample-level QC',
    'sample_autoqc_total': 'A sample-level QC',
    'sample_contamination_pct': 'A sample-level QC',
    'sample_coverage_exceptions': 'A sample-level QC',
    'sample_female_chry_molecules': 'A sample-level QC',
    'sample_gc_bias': 'A sample-level QC',
    'sample_gender_status_mismatch': 'A sample-level QC',
    'sample_germline_contamination': 'A sample-level QC',
    'sample_non_singleton_families': 'A sample-level QC',
    'sample_on_target_rate': 'A sample-level QC',
    'pathogenic_germline_num_status': 'A sample-level QC',
    'sample_autoqc_total_status': 'A sample-level QC',
    'sample_contamination_pct_status': 'A sample-level QC',
    'sample_coverage_exceptions_status': 'A sample-level QC',
    'sample_female_chry_molecules_status': 'A sample-level QC',
    'sample_gc_bias_status': 'A sample-level QC',
    'sample_gender_status_mismatch_status': 'A sample-level QC',
    'sample_germline_contamination_status': 'A sample-level QC',
    'sample_non_singleton_families_status': 'A sample-level QC',
    'sample_on_target_rate_status': 'A sample-level QC',
    'xtr_quant': 'The yield of cfDNA from the extraction process (5ng). NA=not applicable.',
    'Extraction Yield': 'The yield of cfDNA from the extraction process (5ng). NA=not applicable.',
    'en_quant': 'The molarity of the prepared library (nM). NA=not applicable.',
    'Enrichment Molarity': 'The molarity of the prepared library (nM). NA=not applicable.',
    'xtr_quant_status': 'The Pass/Fail status of xtr_quant.',
    'en_quant_status': 'The Pass/Fail status of en_quant',

    'Call': '0=detected but not reported, 1=reported (for SNVs, indels, and fusions), '
            '1=amplification without information to identify it as focal or whole-chromosome-arm, '
            '2=focal amplification, 3=whole-chromosome-arm amplification (for CNVs)',
    'Alteration Type': 'Four alteration types: SNV, Indel, CNV, Fusion',
    'Clinical Variant Class': 'The clinically relevant category according to D-000039',
    'Unique Variant Identifier': 'A concatenation of loci information to uniquely identify '
                                 'variant',
    'MAF/CN': 'the minor allele frequency for SNVs, indels, and fusions; the copy number for CNAs',
    'Mutant Molecules': 'the number of mutant molecules supporting the called variant',
    'Total Molecules': 'the total number of molecules at variant site',
    'Nucleotide change': 'the nucelotide changes for the called variant',
    'QC Note': 'Additional notes about sample processing. A blank cell indicates there were no '
               'additional QC flags'
}
