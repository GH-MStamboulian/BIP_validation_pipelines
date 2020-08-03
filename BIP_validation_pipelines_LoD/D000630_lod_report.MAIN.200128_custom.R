library(rstudioapi)
require(optparse)
setwd(getwd())

# run LoD
option_list <- list(
  make_option(c("-i", "--input_dir"),
              type = "character",
              default = NULL,
              help = "directory containing the bip outputs for all the flowcells involved in this study",
              metavar = "character"),
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
              metavar = "character"),
  make_option(c("-a", "--abs_dir"),
              type = "character",
              default = NULL,
              help = "specify the absolute directory for the scripts",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


resdir = saveDir <- opt$out_dir
dataDir <- opt$data_dir
bip_output = repo <- opt$input_dir#directory for the bip output files as inputs for these scripts
bip_version <-opt$version
scripts_dir <- opt$abs_dir

#bip_output = repo <- "/ghds/ivd/flowcentral/"
#resdir = saveDir <- "/home/mstamboulian/compare_R_versions/R_3.6.2_bip_3.5.3/"
#dataDir <- "/home/mstamboulian/D000630_data_bip_3_5_3/"
#bip_version <- "3.5.3"
#scripts_dir <- getwd()


if(!file.exists(saveDir)) dir.create(saveDir, recursive = TRUE)

library(ghutils)
library(dplyr)
library(ghcnv)
library(ggplot2)
library(data.table)
library(jsonlite)

cat("Done loading dependencies...")
getwd()


#######################
##  Read annotation
######################

type <- read.csv(file.path(dataDir, "sample_list.181204.csv"), header=T, stringsAsFactors=F)
type$Condition = type$Dilution = type$Sample_Description4
type$run_sample_id = type$SAMPLE_ID
type$ng_seq = type$Sample_Description2
type$Type = "poolA" #first set everything to pool A
type$Type[type$Sample_Description5 == "Pool2"] = "poolB" #then set for the rows with Pool2 as PoolB

type = type[grepl("Pool", type$Historical.Sample.ID.Notes),]#removes all the control samples from the list


variants.5ng = read.csv(paste(dataDir , "variants_LoD_5ng_v3.2.csv", sep = '/'), header=T, stringsAsFactors=F)
variants.30ng = read.csv(paste(dataDir , "variants_LoD_30ng_v3.2.csv", sep = '/'), header=T, stringsAsFactors=F)
targeted.maf = targeted.maf.init = read.csv(file=file.path(dataDir,  "targeted_maf.190404.csv"), header=T, stringsAsFactors=F)

uvPoolA = variants.5ng$target_variant[variants.5ng$Pool == "Pool1"]
uvPoolB = variants.5ng$target_variant[variants.5ng$Pool == "Pool2"]

table(type$Dilution, type$Type, type$ng_seq)

#########################
##  Load data from repo
#########################

if (bip_version == '3.5.3')
{
  parameterSetJson <- "/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/repositories/ghpipeline/parameter_sets/GH2.11.json"
  probesBed <- "/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/repositories/ghpipeline/parameter_sets/G360/v2.11/GH2.11_probes.bed"
  bip.version = "3.5.3-0-g8857b98"
}else if (bip_version == '3.5.0' | bip_version == '3.5')
{
  parameterSetJson <- "/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/repositories/ghpipeline/parameter_sets/GH2.11.json"
  probesBed <- "/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/repositories/ghpipeline/parameter_sets/G360/v2.11/GH2.11_probes.bed"
  bip.version = "3.5-0-g15db6a0"
}

#change the if to 0 if the data is allready loaded before and created as an R object then you can just load the R object,
#make this if statement smarter, check whether or not the R object exists, if so then don't enter the if statement
if(!file.exists(file.path(resdir, "work.RData")) ) {
  algoVersion = "v3.5"
  parameterSet <- fromJSON(parameterSetJson)
  probes = getCnvProbes(probeFile=probesBed, algoVersion = algoVersion)
  params = getCnvParams (jsonFile = parameterSetJson, algoVersion = algoVersion)
  
  tables = c("snv_call", "indel_call", "fusion_call", "cnv_call", "sample_coverage", "ghcnv_qc", "qc_on_target", "gh_sample")
  files.ext = c(".snv_call.hdr.tsv", ".indel_call.hdr.tsv", ".fusion_call.hdr.tsv", ".cnv_call.hdr.tsv", ".coverage.hdr.tsv", 
                ".ghcnv_qc.hdr.tsv", ".on_target_db.hdr.tsv", ".gh_sample_db.hdr.tsv")
  data = list()
  
  run_sample_id = type$run_sample_id
  runid = type$runid
  folderid = type$Folder
  for(n in files.ext) data[[tables[which(n == files.ext)]]] = 
    import_data_from_repo(run_sample_id, folderid, repo=repo, params, probes, 
                          what = n, sep="\t", version = bip.version, mc.cores=8)
  names(data) = tables
  data = lapply(tables, function(n) data[[n]] = rbindlist(data[[n]]))
  names(data) = tables
  
  snv = data$snv_call
  snv$mut_key = paste(snv$gene, snv$mut_aa, sep = ".")
  # dedup flowcell
  key = paste(snv$run_sample_id, snv$runid, snv$gene, snv$position, snv$mut_aa, sep="_")
  snv = snv[match(unique(key), key),]
  print(paste("Pulled ", nrow(snv), " SNVs for ", length(unique(snv$run_sample_id)), " samples"))
  
  indel = data$indel_call
  indel$mut_key = paste(indel$gene, indel$mut_aa, sep = ".")
  # dedup flowcell
  key = paste(indel$run_sample_id, indel$runid, indel$mut_key, indel$position, sep="_")
  indel = indel[match(unique(key), key),]
  print(paste("Pulled ", nrow(indel), " indels for ", length(unique(indel$run_sample_id)), " samples"))
  
  fusion = data$fusion_call
  fusion$mut_key = paste(fusion$gene_a, fusion$pos_a, sep = ".")
  key = paste(fusion$run_sample_id, fusion$runid, fusion$mut_key, fusion$gene_b, fusion$pos_b, sep="_")
  fusion = fusion[match(unique(key), key),]
  print(paste("Pulled ", nrow(fusion), " fusions for ", length(unique(fusion$run_sample_id)), " samples"))
  
  cnv = data$cnv_call
  cnv$copy_number = 2*cnv$mean + 2;
  cnv$mut_key = paste(cnv$gene, "CNV", sep = ".")
  key = paste(cnv$run_sample_id, cnv$runid, cnv$gene, sep="_")
  cnv = cnv[match(unique(key), key),]
  print(paste("Pulled ", nrow(cnv), " CNVs for ", length(unique(cnv$run_sample_id)), " samples"))
  
  nsc = data$sample_coverage
  key = paste(nsc$run_sample_id, nsc$runid, nsc$chrom, nsc$start_pos, nsc$end_pos, sep="_")
  nsc = nsc[match(unique(key), key),]
  print(paste("Pulled ", nrow(nsc), " NSCs for ", length(unique(nsc$run_sample_id)), " samples"))
  
  autoqc = lapply(file.path(repo, unique(runid), "autoqc_sample_qc.hdr.tsv"), function(f) fread(f))
  autoqc = do.call(rbind, autoqc)
  flqc = autoqc %>% filter(category == "flowcell")
  control = autoqc %>% filter(category == "control")
  autoqc = autoqc %>% filter(run_sample_id %in% type$run_sample_id)
  print(paste("Pulled ", nrow(autoqc), " autoqc for ", length(unique(autoqc$run_sample_id)), " samples"))
  
  #pull ghcnv qc...
  ghcnv_qc = data$ghcnv_qc
  key = paste(ghcnv_qc$run_sample_id, ghcnv_qc$runid, sep="_")
  ghcnv_qc = ghcnv_qc[match(unique(key), key),]
  print(paste("Pulled ", nrow(ghcnv_qc), " CNV QCs for ", length(unique(ghcnv_qc$run_sample_id)), " samples"))
  
  #pull on_target qc...
  ont = data$qc_on_target
  key = paste(ont$run_sample_id, ont$runid, sep="_")
  ont = ont[match(unique(key), key),]
  print(paste("Pulled ", nrow(ont), " on-target qc for ", length(unique(ont$run_sample_id)), " samples"))
  
  #pull gh_sample qc...
  gh_sample = data$gh_sample
  key = paste(gh_sample$run_sample_id, gh_sample$runid, sep="_")
  gh_sample = gh_sample[match(unique(key), key),]
  print(paste("Pulled ", nrow(gh_sample), " gh_sample qc for ", length(unique(gh_sample$run_sample_id)), " samples"))
  
  # QC metrics
  anno = ont %>% left_join(ghcnv_qc, by=c("run_sample_id", "runid")) %>% left_join(gh_sample, by=c("run_sample_id", "runid"))
  
  save(snv, indel, fusion, cnv, nsc, autoqc, flqc, control, anno, file=file.path(resdir, paste("work", "RData", sep=".")))
}
load(file.path(resdir, "work.RData"))

type = type[match(anno$run_sample_id, type$run_sample_id),]
table(type$Dilution, type$Type, type$ng_seq)


###############################
##  Table 6&7: TRUE MAF values
##  LBP70 and G360CDx methods
###############################


#Table7.targeted.maf.level.G360CDx.tsv Needed
source(paste(scripts_dir, "D000630_lod_report.targeted_maf.200128.R", sep = '/')) #checked all the files have same values as original run

source(paste(scripts_dir,"D000630_lod_report.Detection.74genes.200128.R", sep = '/')) #checked all the files have same values as original run

source(paste(scripts_dir, "D000630_lod_report.Detection.55genes.200128.R", sep = '/')) #checked all the files have same values as original run

#Table5.LoD.part2.G360CDx.74genes.tsv  => requires Table7.targeted.maf.G360CDx.RData by D000630_lod_report.targeted_maf.200128.R
source(paste(scripts_dir, "D000630_lod_report.LoD.G360CDX.74genes.200128.R", sep = '/')) #values differ

#Table5.LoD.part2.G360CDx.55genes.tsv & needed
source(paste(scripts_dir, "D000630_lod_report.LoD.G360CDX.55genes.200128.R", sep = '/'))

