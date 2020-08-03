suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
library(rstudioapi)
require(optparse)
options(repr.plot.width=6, repr.plot.height=5)

if (!require("binom")) install.packages("binom")

`%notin%` <- Negate(`%in%`)

option_list <- list(
  make_option(c("-d", "--data_dir"),
              type = "character",
              default = NULL,
              help = "directory for the list of samples and variants tables of interest for this study",
              metavar = "character"),
  make_option(c("-o", "--out_dir"),
              type = "character",
              default = NULL,
              help = "output directory to store the tables",
              metavar = "character"),
  make_option(c("-v", "--version"),
              type = "character",
              default = NULL,
              help = "specify BIP version to process the data",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


#bip_version <- "3.5.0"
bip_version <- opt$version
#saveDir = output_dir <- "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/precision_comparison/test_out_3/"
saveDir = output_dir <- opt$out_dir
if(!file.exists(saveDir)) dir.create(saveDir, recursive = TRUE)

#outdir <- "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/precision_data/bip_3_5_3/"
#outdir <- "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/precision_data/bip_3_5_2/"
#outdir <- "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/precision_data/bip_3_5_0/"
outdir = data_dir <- opt$data_dir

sample <- read.table(file.path(outdir, "sample_sheet.txt"), header=T, sep="\t")

colnames(sample) <- c("Date", "OWNER", "INSTRUMENT", "FLOWCELL_ID", "Study.Name", "SAMPLE_ID", "Historical.Sample.ID.Notes", "Sample_input_ng", "Condition", "Batch", "Flowcell_ID", "fc_loc", "fc_id", "sequencer", "pool")

sanity_check <- function(df){
  sc <- TRUE
  ncombination <- length(unique(df$Condition))
  print(paste0("total number of combinations is ", ncombination))
  replicates_num <- df %>% group_by(Batch, pool) %>% summarize(replicates = n())
  nbatch <- length(unique(replicates_num$Batch))
  n_replicates = nbatch*5
  print(paste0("total number of batch is ", nbatch))
  if(nbatch != 6){sc <- FALSE}
  pool_num <- replicates_num %>% group_by(Batch) %>% summarize(npool = n())
  for(i in 1:nrow(pool_num)){
    print(paste0("Batch ", as.character(pool_num$Batch[i]), " has ", pool_num$npool[i], " pools"))
    if(pool_num$npool[i] != 3){sc <- FALSE}
  }
  if(! all.equal(replicates_num$replicates, rep(5, 18))){
    sc <- FALSE
  } else {print(paste0("each pool in each of ", nbatch, " batches has 5 replicates, making a total of ", nbatch*5, " replicates for each pool."))}
  return(list("sc" = sc, "n_replicates" = n_replicates))
}


sc_res <- sanity_check(sample)
sc_res$sc
n_replicates <- sc_res$n_replicates

autoqc_sample <- function(fc_loc){
  autoqc_sample <- fread(file.path(fc_loc, "autoqc_sample_qc.hdr.tsv"), header=T, sep="\t")
  return(autoqc_sample)
}

fc_loc_unique <- unique(sample$fc_loc)
fc_loc_unique <- fc_loc_unique[!is.na(fc_loc_unique)]
length(fc_loc_unique)

auto_qc <- lapply(fc_loc_unique, autoqc_sample)
auto_qc_all <- do.call(rbind, auto_qc)

flowcell_qc <- auto_qc_all %>% filter(category == "flowcell") %>% select(-run_sample_id)

control_qc <- auto_qc_all %>% filter(category == "control")

sample_qc <- auto_qc_all %>% filter(category == "sample" & run_sample_id %in% sample$SAMPLE_ID)

write.csv(flowcell_qc, file=file.path(outdir, "LineData_fc_QC.csv"), row.names=F, quote=F)
write.csv(control_qc, file=file.path(outdir, "LineData_control_QC.csv"), row.names=F, quote=F)
write.csv(sample_qc, file=file.path(outdir, "LineData_sample_QC.csv"), row.names=F, quote=F)

flowcell_qc %>% filter(status == "FAIL")
control_qc %>% filter(status == "FAIL")

sample %>% filter(is.na(fc_loc))

sample_qc %>% filter(status == "FAIL" & metric != "sample_autoqc_total")

nrow(sample_qc %>% filter(status == "REVIEW" & metric != "sample_autoqc_total" & ! verbose_name %in% c("Gender Status Mismatch",
                                                                                                       "Germline contamination"
)) %>% select(run_sample_id) %>% distinct())

nrow(sample)

sample$QC_status <- "Pass"

sample_manifest <- sample %>% select(-Date, -OWNER, -INSTRUMENT, -Study.Name, -Historical.Sample.ID.Notes, -Flowcell_ID)
write.csv(sample_manifest, file=file.path(outdir, "LineData_sample_manifest.csv"), row.names=F)

####################
## VARIANT LEVEL PPA
####################

snv_pos_sample <- function(fc_loc, sample){
  snv_call_sample <- read.table(file.path(fc_loc, paste0(sample, ".snv_call.hdr.tsv")), header=T, sep="\t")
  snv_pos_sample <- snv_call_sample %>% filter(! grepl("not reportable", variant_comment))
  return(snv_pos_sample)
}

indel_pos_sample <- function(fc_loc, sample){
  indel_call_sample <- read.table(file.path(fc_loc, paste0(sample, ".indel_call.hdr.tsv")), header=T, sep="\t")
  indel_pos_sample <- indel_call_sample %>% filter(! (grepl("oncogene fs var_ds <= 1", variant_comment) | grepl("single molecule indel present", variant_comment)) & ! grepl("not reportable", variant_comment))
  return(indel_pos_sample)
}

cnv_pos_sample <- function(fc_loc, sample){
  cnv_call_sample <- read.table(file.path(fc_loc, paste0(sample, ".cnv_call.hdr.tsv")), header=T, sep="\t")
  cnv_pos_sample <- cnv_call_sample %>% filter(call > 0)
  return(cnv_pos_sample)
}

fusion_pos_sample <- function(fc_loc, sample){
  fusion_call_sample <- read.table(file.path(fc_loc, paste0(sample, ".fusion_call.hdr.tsv")), header=T, sep="\t")
  fusion_pos_sample <- fusion_call_sample %>% filter(call == 1)
  return(fusion_pos_sample)
}


snv_pos <- apply(sample, 1, function(s) snv_pos_sample(s[12], s[6]))
snv_pos_all <- do.call(rbind, snv_pos)
indel_pos <- apply(sample, 1, function(s) indel_pos_sample(s[12], s[6]))
indel_pos_all <- do.call(rbind, indel_pos)
cnv_pos <- apply(sample, 1, function(s) cnv_pos_sample(s[12], s[6]))
cnv_pos_all <- do.call(rbind, cnv_pos)
fusion_pos <- apply(sample, 1, function(s) fusion_pos_sample(s[12], s[6]))
fusion_pos_all <- do.call(rbind, fusion_pos)


snv_pos_all_class <- snv_pos_all %>% left_join(sample, by=c("run_sample_id" = "SAMPLE_ID"))
indel_pos_all_class <- indel_pos_all %>% left_join(sample, by=c("run_sample_id" = "SAMPLE_ID"))
cnv_pos_all_class <- cnv_pos_all %>% left_join(sample, by=c("run_sample_id" = "SAMPLE_ID"))
fusion_pos_all_class <- fusion_pos_all %>% left_join(sample, by=c("run_sample_id" = "SAMPLE_ID"))


expected_pos_snv <- c("EGFR.L858R.Pool1", "EGFR.T790M.Pool1", "KRAS.G12V.Pool1", "NRAS.Q61R.Pool1", "BRAF.V600E.Pool1", "BRCA1.S1140G.Pool2", "BRCA2.I2944F.Pool2")
expected_pos_indel <- c("EGFR.E746_A750del.Pool1", "ERBB2.A775_G776insYVMA.Pool1","EGFR.A767_V769dup.Pool2", "BRCA1.E23fs.Pool2", "BRCA2.S1982fs.Pool2")
expected_pos_cnv <- c("ERBB2.Pool1", "MET.Pool1")
expected_pos_fusion <- c("RET.Pool1", "ALK.Pool2", "ROS1.Pool2", "NTRK1.Pool3")

expected_neg_snv <- c("BRCA1.S1140G.Pool1", "BRCA2.I2944F.Pool1", "EGFR.L858R.Pool2", "EGFR.T790M.Pool2", "NRAS.Q61R.Pool2", "EGFR.L858R.Pool3", "EGFR.T790M.Pool3", "KRAS.G12V.Pool3", "NRAS.Q61R.Pool3", "BRAF.V600E.Pool3", "BRCA1.S1140G.Pool3", "BRCA2.I2944F.Pool3")
expected_neg_indel <- c("EGFR.A767_V769dup.Pool1", "BRCA1.E23fs.Pool1", "EGFR.E746_A750del.Pool2", "EGFR.E746_A750del.Pool3", "ERBB2.A775_G776insYVMA.Pool3", "EGFR.A767_V769dup.Pool3", "BRCA1.E23fs.Pool3", "BRCA2.S1982fs.Pool3")
expected_neg_cnv <- c("ERBB2.Pool2", "MET.Pool2", "ERBB2.Pool3", "MET.Pool3")
expected_neg_fusion <- c("RET.Pool2", "RET.Pool3", "ALK.Pool1", "ALK.Pool3", "ROS1.Pool1", "ROS1.Pool3", "NTRK1.Pool1", "NTRK1.Pool2")


snv_pos_all_class$var_key <- paste0(snv_pos_all_class$gene, ".", snv_pos_all_class$mut_aa, ".", snv_pos_all_class$pool)
indel_pos_all_class$var_key <- paste0(indel_pos_all_class$gene, ".", indel_pos_all_class$mut_aa_short, ".", indel_pos_all_class$pool)
cnv_pos_all_class$var_key <- paste0(cnv_pos_all_class$gene, ".", cnv_pos_all_class$pool)
fusion_pos_all_class$var_key <- paste0(fusion_pos_all_class$gene_a, ".", fusion_pos_all_class$pool)

## check alterative TP EGFR exon 19 deletion according to D-000319
indel_pos_all_class %>% filter(chrom == "7" & position == "55242466" & mut_nt == "GAATTAAGAGAAGCAA>G")

snv_pos_all_class %>% filter(chrom == "7" & position == "55242466" & mut_nt == "G>A")


indel_pos_all_class$var_key[indel_pos_all_class$run_sample_id == "A014107801" & 
                              indel_pos_all_class$chrom == "7" &
                              indel_pos_all_class$position == "55242466" &
                              indel_pos_all_class$mut_nt == "GAATTAAGAGAAGCAA>G"
                            ] <- "EGFR.E746_A750del.Pool1"
indel_pos_all_class$percentage[indel_pos_all_class$run_sample_id == "A014107801" & 
                                 indel_pos_all_class$chrom == "7" &
                                 indel_pos_all_class$position == "55242466" &
                                 indel_pos_all_class$mut_nt == "GAATTAAGAGAAGCAA>G"
                               ] <- NA
indel_pos_all_class$mut_aa_short[indel_pos_all_class$run_sample_id == "A014107801" & 
                                   indel_pos_all_class$chrom == "7" &
                                   indel_pos_all_class$position == "55242466" &
                                   indel_pos_all_class$mut_nt == "GAATTAAGAGAAGCAA>G"
                                 ] <- "E746_A750del"
indel_pos_all_class$mut_aa[indel_pos_all_class$run_sample_id == "A014107801" & 
                             indel_pos_all_class$chrom == "7" &
                             indel_pos_all_class$position == "55242466" &
                             indel_pos_all_class$mut_nt == "GAATTAAGAGAAGCAA>G"
                           ] <- "p.Glu746_Ala750del"

snv_expected_pos_neg_all <- snv_pos_all_class %>% filter(var_key %in% c(expected_pos_snv, expected_neg_snv))
indel_expected_pos_neg_all <- indel_pos_all_class %>% filter(var_key %in% c(expected_pos_indel, expected_neg_indel))
cnv_expected_pos_neg_all <- cnv_pos_all_class %>% filter(var_key %in% c(expected_pos_cnv, expected_neg_cnv))
fusion_expected_pos_neg_all <- fusion_pos_all_class %>% filter(var_key %in% c(expected_pos_fusion, expected_neg_fusion))


snv_expected_pos_neg_all %>% group_by(var_key) %>% summarize(mean_MAF = mean(percentage))

indel_expected_pos_neg_all %>% group_by(var_key) %>% summarize(mean_MAF = mean(percentage, na.rm = T))

cnv_expected_pos_neg_all %>% group_by(var_key) %>% summarize(mean_CNV = mean(copy_number, na.rm = T))

fusion_expected_pos_neg_all %>% group_by(var_key) %>% summarize(mean_MAF = mean(percentage, na.rm = T))

snv_ppa_df <- snv_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_snv) %>% 
  group_by(var_key) %>% 
  summarize(number_pos_var = n())
snv_ppa_df

indel_ppa_df <- indel_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_indel) %>% 
  group_by(var_key) %>% 
  summarize(number_pos_var = n())
indel_ppa_df


cnv_ppa_df <- cnv_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_cnv) %>% 
  group_by(var_key) %>% 
  summarize(number_pos_var = n())
cnv_ppa_df



fusion_ppa_df <- fusion_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_fusion) %>% 
  group_by(var_key) %>% 
  summarize(number_pos_var = n())
fusion_ppa_df


#for false negative calls (i.e. expected but not reported) all variants

missed_variants_list <- list()

snv_ppa_df_df <- data.frame(var_key = snv_ppa_df[1]$var_key, number_pos_var = snv_ppa_df[2]$number_pos_var)
for (i in seq_len(nrow(snv_ppa_df_df))){
  if(snv_ppa_df_df[i, ]$number_pos_var < n_replicates){
    n_missed <- n_replicates - snv_ppa_df_df[i, ]$number_pos_var
    variant <- as.character(snv_ppa_df_df$var_key[[i]])
    row <- list('FUSION', variant, n_missed)
    missed_variants_list <- append(missed_variants_list, list(row))
  }
}

indel_ppa_df_df <- data.frame(var_key = indel_ppa_df[1]$var_key, number_pos_var = indel_ppa_df[2]$number_pos_var)
for (i in seq_len(nrow(indel_ppa_df_df))){
  if(indel_ppa_df_df[i, ]$number_pos_var < n_replicates){
    n_missed <- n_replicates - indel_ppa_df_df[i, ]$number_pos_var
    variant <- as.character(indel_ppa_df_df$var_key[[i]])
    row <- list('INDEL', variant, n_missed)
    missed_variants_list <- append(missed_variants_list, list(row))
  }
}


cnv_ppa_df_df <- data.frame(var_key = cnv_ppa_df[1]$var_key, number_pos_var = cnv_ppa_df[2]$number_pos_var)
for (i in seq_len(nrow(cnv_ppa_df_df))){
  if(cnv_ppa_df_df[i, ]$number_pos_var < n_replicates){
    n_missed <- n_replicates - cnv_ppa_df_df[i, ]$number_pos_var
    variant <- as.character(cnv_ppa_df_df$var_key[[i]])
    row <- list('CNA', variant, n_missed)
    missed_variants_list <- append(missed_variants_list, list(row))
  }
}


fusion_ppa_df_df <- data.frame(var_key = fusion_ppa_df[1]$var_key, number_pos_var = fusion_ppa_df[2]$number_pos_var)
for (i in seq_len(nrow(fusion_ppa_df_df))){
  if(fusion_ppa_df_df[i, ]$number_pos_var < n_replicates){
    n_missed <- n_replicates - fusion_ppa_df_df[i, ]$number_pos_var
    variant <- as.character(fusion_ppa_df_df$var_key[[i]])
    row <- list('FUSION', variant, n_missed)
    missed_variants_list <- append(missed_variants_list, list(row))
  }
}


n_pos_SNV <- sum(snv_ppa_df$number_pos_var)
n_pos_INDEL <- sum(indel_ppa_df$number_pos_var)
n_pos_CNV <- sum(cnv_ppa_df$number_pos_var)
n_pos_FUSION <- sum(fusion_ppa_df$number_pos_var)

total_expectPos_SNV <- length(expected_pos_snv)*30
total_expectPos_INDEL <- length(expected_pos_indel)*30
total_expectPos_CNV <- length(expected_pos_cnv)*30
total_expectPos_FUSION <- length(expected_pos_fusion)*30

print (paste0("detected ", n_pos_SNV, " SNVs from ", total_expectPos_SNV, " samples."))
print (paste0("detected ", n_pos_INDEL, " Indels from ", total_expectPos_INDEL, " samples."))
print (paste0("detected ", n_pos_CNV, " CNVs from ", total_expectPos_CNV, " samples."))
print (paste0("detected ", n_pos_FUSION, " fusions from ", total_expectPos_FUSION, " samples."))

snv_ppa <- n_pos_SNV / total_expectPos_SNV
indel_ppa <- n_pos_INDEL / total_expectPos_INDEL
cnv_ppa <- n_pos_CNV / total_expectPos_CNV
fusion_ppa <- n_pos_FUSION / total_expectPos_FUSION


library(binom)
CI_0.95_lower <- function(n_pos, n_total){
  CI_0.95_lower = binom.confint(n_pos, n_total, conf.level = 0.95, methods = "exact")[1,5]
  return(CI_0.95_lower)
} 
CI_0.95_upper <- function(n_pos, n_total){
  CI_0.95_upper = binom.confint(n_pos, n_total, conf.level = 0.95, methods = "exact")[1,6]
  return(CI_0.95_upper)
}

CI_0.95_lower_snv_ppa <- CI_0.95_lower(n_pos_SNV, total_expectPos_SNV)
CI_0.95_upper_snv_ppa <- CI_0.95_upper(n_pos_SNV, total_expectPos_SNV)


CI_0.95_lower_indel_ppa <- CI_0.95_lower(n_pos_INDEL, total_expectPos_INDEL)
CI_0.95_upper_indel_ppa <- CI_0.95_upper(n_pos_INDEL, total_expectPos_INDEL)


CI_0.95_lower_cnv_ppa <- CI_0.95_lower(n_pos_CNV, total_expectPos_CNV)
CI_0.95_upper_cnv_ppa <- CI_0.95_upper(n_pos_CNV, total_expectPos_CNV)


CI_0.95_lower_fusion_ppa <- CI_0.95_lower(n_pos_FUSION, total_expectPos_FUSION)
CI_0.95_upper_fusion_ppa <- CI_0.95_upper(n_pos_FUSION, total_expectPos_FUSION)

var_ppa <- data.frame(c("SNV", "INDEL", "CNV", "FUSION"), c(n_pos_SNV, n_pos_INDEL, n_pos_CNV, n_pos_FUSION), c(total_expectPos_SNV, total_expectPos_INDEL, total_expectPos_CNV, total_expectPos_FUSION), c(snv_ppa, indel_ppa, cnv_ppa, fusion_ppa), c(CI_0.95_lower_snv_ppa, CI_0.95_lower_indel_ppa, CI_0.95_lower_cnv_ppa, CI_0.95_lower_fusion_ppa), c(CI_0.95_upper_snv_ppa, CI_0.95_upper_indel_ppa, CI_0.95_upper_cnv_ppa, CI_0.95_upper_fusion_ppa))
colnames(var_ppa) <- c("variant_class", "number of positive calls", "total number of expected calls", "PPA", "CI_0.95_lower", "CI_0.95_upper")
var_ppa


write.table(var_ppa, file = file.path(output_dir, "Table7_variant_class_level_PPA.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

false_negatives_df <- do.call(rbind, missed_variants_list)
colnames(false_negatives_df) <- c("variant_class", "variant", "number of false negatives")
row.names(false_negatives_df) <- 1:nrow(false_negatives_df)
false_negatives_df

write.table(false_negatives_df, file = file.path(output_dir, "Table7_variant_class_level_FalseNegatives_PPA.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)



#for false negative calls (i.e. expected but not reported) rm_reportable only

if(bip_version == "3.5.3"){
  snv_expected_pos_neg_all_rm_reportable <- snv_expected_pos_neg_all %>% filter(rm_reportable == 1)
  indel_expected_pos_neg_all_rm_reportable <- indel_expected_pos_neg_all %>% filter(rm_reportable == 1)
  cnv_expected_pos_neg_all_rm_reportable <- cnv_expected_pos_neg_all %>% filter(call == 2 | call == 3)
  fusion_expected_pos_neg_all_rm_reportable <- fusion_expected_pos_neg_all %>% filter(call == 1)
  
  
  snv_ppa_df_rm_reportable <- snv_expected_pos_neg_all_rm_reportable %>% 
    filter(var_key %in% expected_pos_snv) %>% 
    group_by(var_key) %>% 
    summarize(number_pos_var = n())
  snv_ppa_df_rm_reportable
  
  
  indel_ppa_df_rm_reportable <- indel_expected_pos_neg_all_rm_reportable %>% 
    filter(var_key %in% expected_pos_indel) %>% 
    group_by(var_key) %>% 
    summarize(number_pos_var = n())
  indel_ppa_df_rm_reportable
  
  
  cnv_ppa_df_rm_reportable <- cnv_expected_pos_neg_all_rm_reportable %>% 
    filter(var_key %in% expected_pos_cnv) %>% 
    group_by(var_key) %>% 
    summarize(number_pos_var = n())
  cnv_ppa_df_rm_reportable
  
  fusion_ppa_df_rm_reportable <- fusion_expected_pos_neg_all_rm_reportable %>% 
    filter(var_key %in% expected_pos_fusion) %>% 
    group_by(var_key) %>% 
    summarize(number_pos_var = n())
  fusion_ppa_df_rm_reportable
  
  
  
  missed_variants_list <- list()
  
  snv_ppa_df_rm_reportable_df <- data.frame(var_key = snv_ppa_df_rm_reportable[1]$var_key, number_pos_var = snv_ppa_df_rm_reportable[2]$number_pos_var)
  for (i in seq_len(nrow(snv_ppa_df_rm_reportable_df))){
    if(snv_ppa_df_rm_reportable_df[i, ]$number_pos_var < n_replicates){
      n_missed <- n_replicates - snv_ppa_df_rm_reportable_df[i, ]$number_pos_var
      variant <- as.character(snv_ppa_df_rm_reportable_df$var_key[[i]])
      row <- list('FUSION', variant, n_missed)
      missed_variants_list <- append(missed_variants_list, list(row))
    }
  }
    
  indel_ppa_df_rm_reportable_df <- data.frame(var_key = indel_ppa_df_rm_reportable[1]$var_key, number_pos_var = indel_ppa_df_rm_reportable[2]$number_pos_var)
  for (i in seq_len(nrow(indel_ppa_df_rm_reportable_df))){
    if(indel_ppa_df_rm_reportable_df[i, ]$number_pos_var < n_replicates){
      n_missed <- n_replicates - indel_ppa_df_rm_reportable_df[i, ]$number_pos_var
      variant <- as.character(indel_ppa_df_rm_reportable_df$var_key[[i]])
      row <- list('INDEL', variant, n_missed)
      missed_variants_list <- append(missed_variants_list, list(row))
    }
  }
  
  
  cnv_ppa_df_rm_reportable_df <- data.frame(var_key = cnv_ppa_df_rm_reportable[1]$var_key, number_pos_var = cnv_ppa_df_rm_reportable[2]$number_pos_var)
  for (i in seq_len(nrow(cnv_ppa_df_rm_reportable_df))){
    if(cnv_ppa_df_rm_reportable_df[i, ]$number_pos_var < n_replicates){
      n_missed <- n_replicates - cnv_ppa_df_rm_reportable_df[i, ]$number_pos_var
      variant <- as.character(cnv_ppa_df_rm_reportable_df$var_key[[i]])
      row <- list('CNA', variant, n_missed)
      missed_variants_list <- append(missed_variants_list, list(row))
    }
  }
  
  
  fusion_ppa_df_rm_reportable_df <- data.frame(var_key = fusion_ppa_df_rm_reportable[1]$var_key, number_pos_var = fusion_ppa_df_rm_reportable[2]$number_pos_var)
  for (i in seq_len(nrow(fusion_ppa_df_rm_reportable_df))){
    if(fusion_ppa_df_rm_reportable_df[i, ]$number_pos_var < n_replicates){
      n_missed <- n_replicates - fusion_ppa_df_rm_reportable_df[i, ]$number_pos_var
      variant <- as.character(fusion_ppa_df_rm_reportable_df$var_key[[i]])
      row <- list('FUSION', variant, n_missed)
      missed_variants_list <- append(missed_variants_list, list(row))
    }
  }
  
  
  
  n_pos_SNV_rm_reportable <- sum(snv_ppa_df_rm_reportable$number_pos_var)
  n_pos_INDEL_rm_reportable <- sum(indel_ppa_df_rm_reportable$number_pos_var)
  n_pos_CNV_rm_reportable <- sum(cnv_ppa_df_rm_reportable$number_pos_var)
  n_pos_FUSION_rm_reportable <- sum(fusion_ppa_df_rm_reportable$number_pos_var)
  
  total_expectPos_SNV_rm_reportable <- (length(expected_pos_snv) - 2)*30 # expected positive doesn't include BRCA1 and BRCA2 SNV as rm_reportable is not expected to be 1
  total_expectPos_INDEL_rm_reportable <- length(expected_pos_indel)*30
  total_expectPos_CNV_rm_reportable <- length(expected_pos_cnv)*30
  total_expectPos_FUSION_rm_reportable <- length(expected_pos_fusion)*30
  
  print (paste0("detected ", n_pos_SNV_rm_reportable, " SNVs from ", total_expectPos_SNV_rm_reportable, " samples."))
  print (paste0("detected ", n_pos_INDEL_rm_reportable, " Indels from ", total_expectPos_INDEL_rm_reportable, " samples."))
  print (paste0("detected ", n_pos_CNV_rm_reportable, " CNVs from ", total_expectPos_CNV_rm_reportable, " samples."))
  print (paste0("detected ", n_pos_FUSION_rm_reportable, " fusions from ", total_expectPos_FUSION_rm_reportable, " samples."))
  
  snv_ppa_rm_reportable <- n_pos_SNV_rm_reportable / total_expectPos_SNV_rm_reportable
  indel_ppa_rm_reportable <- n_pos_INDEL_rm_reportable / total_expectPos_INDEL_rm_reportable
  cnv_ppa_rm_reportable <- n_pos_CNV_rm_reportable / total_expectPos_CNV_rm_reportable
  fusion_ppa_rm_reportable <- n_pos_FUSION_rm_reportable / total_expectPos_FUSION_rm_reportable
  
  
  
  CI_0.95_lower_snv_ppa_rm_reportable <- CI_0.95_lower(n_pos_SNV_rm_reportable, total_expectPos_SNV_rm_reportable)
  CI_0.95_upper_snv_ppa_rm_reportable <- CI_0.95_upper(n_pos_SNV_rm_reportable, total_expectPos_SNV_rm_reportable)
  
  
  CI_0.95_lower_indel_ppa_rm_reportable <- CI_0.95_lower(n_pos_INDEL_rm_reportable, total_expectPos_INDEL_rm_reportable)
  CI_0.95_upper_indel_ppa_rm_reportable <- CI_0.95_upper(n_pos_INDEL_rm_reportable, total_expectPos_INDEL_rm_reportable)
  
  
  CI_0.95_lower_cnv_ppa_rm_reportable <- CI_0.95_lower(n_pos_CNV_rm_reportable, total_expectPos_CNV_rm_reportable)
  CI_0.95_upper_cnv_ppa_rm_reportable <- CI_0.95_upper(n_pos_CNV_rm_reportable, total_expectPos_CNV_rm_reportable)
  
  
  CI_0.95_lower_fusion_ppa_rm_reportable <- CI_0.95_lower(n_pos_FUSION_rm_reportable, total_expectPos_FUSION_rm_reportable)
  CI_0.95_upper_fusion_ppa_rm_reportable <- CI_0.95_upper(n_pos_FUSION_rm_reportable, total_expectPos_FUSION_rm_reportable)
  
  
  var_ppa_rm_reportable <- data.frame(c("SNV", "INDEL", "CNV", "FUSION"), c(n_pos_SNV_rm_reportable, n_pos_INDEL_rm_reportable, n_pos_CNV_rm_reportable, n_pos_FUSION_rm_reportable), c(total_expectPos_SNV_rm_reportable, total_expectPos_INDEL_rm_reportable, total_expectPos_CNV_rm_reportable, total_expectPos_FUSION_rm_reportable), c(snv_ppa_rm_reportable, indel_ppa_rm_reportable, cnv_ppa_rm_reportable, fusion_ppa_rm_reportable), c(CI_0.95_lower_snv_ppa_rm_reportable, CI_0.95_lower_indel_ppa_rm_reportable, CI_0.95_lower_cnv_ppa_rm_reportable, CI_0.95_lower_fusion_ppa_rm_reportable), c(CI_0.95_upper_snv_ppa_rm_reportable, CI_0.95_upper_indel_ppa_rm_reportable, CI_0.95_upper_cnv_ppa_rm_reportable, CI_0.95_upper_fusion_ppa_rm_reportable))
  colnames(var_ppa_rm_reportable) <- c("variant_class", "number of positive calls", "total number of expected calls", "PPA", "CI_0.95_lower", "CI_0.95_upper")
  var_ppa_rm_reportable
  
  write.table(var_ppa_rm_reportable, file = file.path(output_dir, "Table7_variant_class_level_PPA_rm_reportable.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  false_negatives_df <- do.call(rbind, missed_variants_list)
  colnames(false_negatives_df) <- c("variant_class", "variant", "number of false negatives")
  row.names(false_negatives_df) <- 1:nrow(false_negatives_df)
  false_negatives_df
  
  write.table(false_negatives_df, file = file.path(output_dir, "Table7_variant_class_level_FalseNegatives_PPA_rm_reportable.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
}

###################

###################
## sample level NPA
###################

##all 74 genes no rm_reportable filtering

sample_neg_snv <- unique(snv_expected_pos_neg_all$run_sample_id[snv_expected_pos_neg_all$var_key %in% expected_neg_snv])
sample_neg_indel <- unique(indel_expected_pos_neg_all$run_sample_id[indel_expected_pos_neg_all$var_key %in% expected_neg_indel])
sample_neg_cnv <- unique(cnv_expected_pos_neg_all$run_sample_id[cnv_expected_pos_neg_all$var_key %in% expected_neg_cnv])
sample_neg_fusion <- unique(fusion_expected_pos_neg_all$run_sample_id[fusion_expected_pos_neg_all$var_key %in% expected_neg_fusion])

sample_neg <- unique(c(sample_neg_snv, sample_neg_indel, sample_neg_cnv, sample_neg_fusion))

n_neg_sample <- length(sample_neg)

n_neg_sample

sample_npa <- 1 - n_neg_sample / 90
sample_npa

sample_npa_df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(sample_npa_df) <- c("Number_samples", "N_samples_with_False_Positives", "sample_level_NPA", "95%_CI")
sample_npa_df[nrow(sample_npa_df) +1, ] = c(90, n_neg_sample, (sample_npa*100), "[96%-100%]")
sample_npa_df

write.table(sample_npa_df, file = file.path(output_dir, "Table10_sample_level_NPA.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


nrow(indel_pos_all_class %>% filter(chrom == "7" & position == "55242466" & mut_nt == "GAATTAAGAGAAGCAA>G" & pool %in% c("Pool2", "Pool3")))

## limit to rm_reportable = 1
## rm_reportables only filtered out above

if(bip_version == "3.5.3"){
  sample_neg_snv_rm_reportable <- unique(snv_expected_pos_neg_all_rm_reportable$run_sample_id[snv_expected_pos_neg_all_rm_reportable$var_key %in% expected_neg_snv])
  sample_neg_indel_rm_reportable <- unique(indel_expected_pos_neg_all_rm_reportable$run_sample_id[indel_expected_pos_neg_all_rm_reportable$var_key %in% expected_neg_indel])
  sample_neg_cnv_rm_reportable <- unique(cnv_expected_pos_neg_all_rm_reportable$run_sample_id[cnv_expected_pos_neg_all_rm_reportable$var_key %in% expected_neg_cnv])
  sample_neg_fusion_rm_reportable <- unique(fusion_expected_pos_neg_all_rm_reportable$run_sample_id[fusion_expected_pos_neg_all_rm_reportable$var_key %in% expected_neg_fusion])
  
  sample_neg_rm_reportable <- unique(c(sample_neg_snv_rm_reportable, sample_neg_indel_rm_reportable, sample_neg_cnv_rm_reportable, sample_neg_fusion_rm_reportable))
  
  n_neg_sample_rm_reportable <- length(sample_neg_rm_reportable)
  
  n_neg_sample_rm_reportable
  
  sample_npa_rm_reportable <- 1 - n_neg_sample_rm_reportable / 90
  sample_npa_rm_reportable
  
  sample_npa_rm_reportable_df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(sample_npa_rm_reportable_df) <- c("Number_samples", "N_samples_with_False_Positives", "sample_level_NPA", "95%_CI")
  sample_npa_rm_reportable_df[nrow(sample_npa_rm_reportable_df) +1, ] = c(90, n_neg_sample_rm_reportable, (sample_npa_rm_reportable*100), "[96%-100%]")
  sample_npa_rm_reportable_df
  
  write.table(sample_npa_rm_reportable_df, file = file.path(output_dir, "Table11_sample_level_NPA_rm_reportable.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

  nrow(indel_expected_pos_neg_all_rm_reportable %>% filter(chrom == "7" & position == "55242466" & mut_nt == "GAATTAAGAGAAGCAA>G" & pool %in% c("Pool2", "Pool3")))

}
#false positives for EGFR 7:55242466:GAATTAAGAGAAGCAA>G + 7:55242466:G>A in Pool 2 and Pool 3


#######################
## WITHIN-CONDITION PPA
#######################

##all 74 genes no rm_reportable filtering
snv_concord <- snv_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_snv) %>% 
  group_by(var_key, pool, Condition) %>%
  summarize(n_sample_corc = n())

indel_concord <- indel_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_indel) %>% 
  group_by(var_key, pool, Condition) %>%
  summarize(n_sample_corc = n())

cnv_concord <- cnv_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_cnv) %>% 
  group_by(var_key, pool, Condition) %>%
  summarize(n_sample_corc = n())

fusion_concord <- fusion_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_fusion) %>% 
  group_by(var_key, pool, Condition) %>%
  summarize(n_sample_corc = n())

fusion_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_fusion & gene_a == 'ALK') %>%
  group_by(gene_a, gene_b, downstream_gene) %>% summarize(n = n())

fusion_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_fusion & gene_a == 'NTRK1') %>%
  group_by(gene_a, gene_b, downstream_gene) %>% summarize(n = n())

fusion_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_fusion & gene_a == 'RET') %>%
  group_by(gene_a, gene_b, downstream_gene) %>% summarize(n = n())

fusion_expected_pos_neg_all %>% 
  filter(var_key %in% expected_pos_fusion & gene_a == 'ROS1') %>%
  group_by(gene_a, gene_b, downstream_gene) %>% summarize(n = n())

concord <- rbind(snv_concord, indel_concord, cnv_concord, fusion_concord)

concord_combination1 <- concord %>% filter(Condition == 'Precision Combination 1')

n_pos_samp_condition1 <- sum(concord_combination1$n_sample_corc)
n_total_samp_condition1 <- 18*10

combination1_ppa <- n_pos_samp_condition1 / n_total_samp_condition1

concord_combination2 <- concord %>% filter(Condition == 'Precision Combination 2')

n_pos_samp_condition2 <- sum(concord_combination2$n_sample_corc)
n_total_samp_condition2 <- 18*10

combination2_ppa <- n_pos_samp_condition2 / n_total_samp_condition2

concord_combination3 <- concord %>% filter(Condition == 'Precision Combination 3')

n_pos_samp_condition3 <- sum(concord_combination3$n_sample_corc)
n_total_samp_condition3 <- 18*10

combination3_ppa <- n_pos_samp_condition3 / n_total_samp_condition3

CI_0.95_lower_combination1_ppa <- CI_0.95_lower(n_pos_samp_condition1, n_total_samp_condition1)
CI_0.95_upper_combination1_ppa <- CI_0.95_upper(n_pos_samp_condition1, n_total_samp_condition1)


CI_0.95_lower_combination2_ppa <- CI_0.95_lower(n_pos_samp_condition2, n_total_samp_condition2)
CI_0.95_upper_combination2_ppa <- CI_0.95_upper(n_pos_samp_condition2, n_total_samp_condition2)


CI_0.95_lower_combination3_ppa <- CI_0.95_lower(n_pos_samp_condition3, n_total_samp_condition3)
CI_0.95_upper_combination3_ppa <- CI_0.95_upper(n_pos_samp_condition3, n_total_samp_condition3)

condition_ppa <- data.frame(c("PC 1", "PC 2", "PC 3"), c(n_pos_samp_condition1, n_pos_samp_condition2, n_pos_samp_condition3), c(n_total_samp_condition1, n_total_samp_condition2, n_total_samp_condition3), c(combination1_ppa, combination2_ppa, combination3_ppa), c(CI_0.95_lower_combination1_ppa, CI_0.95_lower_combination2_ppa, CI_0.95_lower_combination3_ppa), c(CI_0.95_upper_combination1_ppa, CI_0.95_upper_combination2_ppa, CI_0.95_upper_combination3_ppa))
colnames(condition_ppa) <- c("Condition", "number of positive calls", "total number of expected calls", "PPA", "CI_0.95_lower", "CI_0.95_upper")
condition_ppa

write.table(condition_ppa, file = file.path(output_dir, "Table12_within_condition_PPA.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


## limit to rm_reportable = 1

if(bip_version == "3.5.3"){
  snv_concord_rm_reportable <- snv_expected_pos_neg_all_rm_reportable %>% 
    filter(var_key %in% expected_pos_snv) %>% 
    group_by(var_key, pool, Condition) %>%
    summarize(n_sample_corc = n())
  
  indel_concord_rm_reportable <- indel_expected_pos_neg_all_rm_reportable %>% 
    filter(var_key %in% expected_pos_indel) %>% 
    group_by(var_key, pool, Condition) %>%
    summarize(n_sample_corc = n())
  
  cnv_concord_rm_reportable <- cnv_expected_pos_neg_all_rm_reportable %>% 
    filter(var_key %in% expected_pos_cnv) %>% 
    group_by(var_key, pool, Condition) %>%
    summarize(n_sample_corc = n())
  
  fusion_concord_rm_reportable <- fusion_expected_pos_neg_all_rm_reportable %>% 
    filter(var_key %in% expected_pos_fusion) %>% 
    group_by(var_key, pool, Condition) %>%
    summarize(n_sample_corc = n())
  
  concord_rm_reportable <- rbind(snv_concord_rm_reportable, indel_concord_rm_reportable, cnv_concord_rm_reportable, fusion_concord_rm_reportable)
  
  concord_combination1_rm_reportable <- concord_rm_reportable %>% filter(Condition == 'Precision Combination 1')
  
  n_pos_samp_condition1_rm_reportable <- sum(concord_combination1_rm_reportable$n_sample_corc)
  n_total_samp_condition1_rm_reportable <- (18 - 2)*10
  
  combination1_ppa_rm_reportable <- n_pos_samp_condition1_rm_reportable / n_total_samp_condition1_rm_reportable
  
  concord_combination2_rm_reportable <- concord_rm_reportable %>% filter(Condition == 'Precision Combination 2')
  
  n_pos_samp_condition2_rm_reportable <- sum(concord_combination2_rm_reportable$n_sample_corc)
  n_total_samp_condition2_rm_reportable <- (18 - 2)*10
  
  combination2_ppa_rm_reportable <- n_pos_samp_condition2_rm_reportable / n_total_samp_condition2_rm_reportable
  
  concord_combination3_rm_reportable <- concord_rm_reportable %>% filter(Condition == 'Precision Combination 3')
  
  n_pos_samp_condition3_rm_reportable <- sum(concord_combination3_rm_reportable$n_sample_corc)
  n_total_samp_condition3_rm_reportable <- (18 - 2)*10
  
  combination3_ppa_rm_reportable <- n_pos_samp_condition3_rm_reportable / n_total_samp_condition3_rm_reportable
  
  CI_0.95_lower_combination1_ppa_rm_reportable <- CI_0.95_lower(n_pos_samp_condition1_rm_reportable, n_total_samp_condition1_rm_reportable)
  CI_0.95_upper_combination1_ppa_rm_reportable <- CI_0.95_upper(n_pos_samp_condition1_rm_reportable, n_total_samp_condition1_rm_reportable)
  
  
  CI_0.95_lower_combination2_ppa_rm_reportable <- CI_0.95_lower(n_pos_samp_condition2_rm_reportable, n_total_samp_condition2_rm_reportable)
  CI_0.95_upper_combination2_ppa_rm_reportable <- CI_0.95_upper(n_pos_samp_condition2_rm_reportable, n_total_samp_condition2_rm_reportable)
  
  
  CI_0.95_lower_combination3_ppa_rm_reportable <- CI_0.95_lower(n_pos_samp_condition3_rm_reportable, n_total_samp_condition3_rm_reportable)
  CI_0.95_upper_combination3_ppa_rm_reportable <- CI_0.95_upper(n_pos_samp_condition3_rm_reportable, n_total_samp_condition3_rm_reportable)
  
  condition_ppa_rm_reportable <- data.frame(c("PC 1", "PC 2", "PC 3"), c(n_pos_samp_condition1_rm_reportable, n_pos_samp_condition2_rm_reportable, n_pos_samp_condition3_rm_reportable), c(n_total_samp_condition1_rm_reportable, n_total_samp_condition2_rm_reportable, n_total_samp_condition3_rm_reportable), c(combination1_ppa_rm_reportable, combination2_ppa_rm_reportable, combination3_ppa_rm_reportable), c(CI_0.95_lower_combination1_ppa_rm_reportable, CI_0.95_lower_combination2_ppa_rm_reportable, CI_0.95_lower_combination3_ppa_rm_reportable), c(CI_0.95_upper_combination1_ppa_rm_reportable, CI_0.95_upper_combination2_ppa_rm_reportable, CI_0.95_upper_combination3_ppa_rm_reportable))
  colnames(condition_ppa_rm_reportable) <- c("Condition", "number of positive calls", "total number of expected calls", "PPA", "CI_0.95_lower", "CI_0.95_upper")
  condition_ppa_rm_reportable
  
  write.table(condition_ppa_rm_reportable, file = file.path(output_dir, "Table13_within_condition_PPA_rm_reportable.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
}  
