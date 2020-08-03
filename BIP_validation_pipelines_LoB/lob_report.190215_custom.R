#also extract some support information as well, such as MAF and mutant molecule counts
#repeat with 3.5.2 as well, 

library(rstudioapi)
require(optparse)

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
g360.repo = bip_output = repo <- opt$input_dir#directory for the bip output files as inputs for these scripts
bip_version <-opt$version
abs_dir <- opt$abs_dir


setwd(abs_dir)
print(getwd())


# Code location
gitDir = "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/repositories"
bipDir = file.path(gitDir, "ghpipeline") #file.path's delimiter is a / autimatically so no need to end your strings with "/"
localDir = file.path(gitDir, "data_science")
#setwd("/host_home/git/data_science/02_DEVELOPMENT/181010_IVD_LoB")

## Data location
#dataDir = "/ghds/groups/ivd_study_data/av/181010_IVD_LoB/data"
#vcfDir = "/ghds/groups/ivd_study_data/av/181010_IVD_LoB/data/vcf"
#dataDir = "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES//lob_data/181010_data_bip_3_5_0"
G360Dir = file.path(dataDir, "G360") #use either this
WESDir = file.path(dataDir, "WES")
LBP70Dir = file.path(dataDir, "LBP70") #or this for reference to filter out the variants first instead of WESDir

#g360.repo = "/ghds/ivd/flowcentral"

#bip_version = "3.5.0"

##  Results location
#saveDir = resdir = "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/bip_3_5_0/lob_report_190215"
if(!file.exists(saveDir)) dir.create(saveDir, recursive = TRUE)

###################
##  Dependencies
###################

library(VariantAnnotation)
library(GenomicRanges)
library(RColorBrewer)
library(MASS)
library(readr)
library(ggplot2)
library(ghutils)
library(ghcnv)
library(data.table)
library(dplyr)
library(jsonlite)
source(file.path(localDir, "ghutils/R/prepareData_GH.R"))

########################
##  Specify BIP version 
##  and parameters
########################

if (bip_version == '3.5.3')
{
  parameterSetJson <- file.path(bipDir, "parameter_sets/GH2.11.json")
  probesBed <- file.path(bipDir, "parameter_sets/G360/v2.11/GH2.11_probes.bed")
  bip.version = "3.5.3-0-g8857b98"
}else if (bip_version == '3.5.0' | bip_version == '3.5')
{
  parameterSetJson <- file.path(bipDir, "parameter_sets/GH2.11.json")
  probesBed <- file.path(bipDir, "parameter_sets/G360/v2.11/GH2.11_probes.bed")
  bip.version = "3.5-0-g15db6a0"
}else if (bip_version == '3.5.2')
{
  parameterSetJson <- file.path(bipDir, "parameter_sets/GH2.11.json")
  probesBed <- file.path(bipDir, "parameter_sets/G360/v2.11/GH2.11_probes.bed")
  bip.version = "3.5.2-rc2-0-g9e05762"
}

########################
##  Load samples
########################

#G360 samples to folders maps
type.s2f <- read.csv(file.path(G360Dir, "sample2folder.csv"), header=T, stringsAsFactors=F)

# G360 neat cfDNA
type.cfDNA <- read.csv(file.path(G360Dir, "sample_list_neat.csv"), header=T, stringsAsFactors=F)
colnames(type.cfDNA)[ncol(type.cfDNA)] = "Donor"
type.cfDNA$run_sample_id = type.cfDNA$SAMPLE_ID
type.cfDNA = type.cfDNA %>% filter(!grepl("Control", Sample_Description3))
#type.cfDNA$Donor[1:9] = paste("Donor", 1:9, sep="_")
type.cfDNA$Study.Name = "cfDNA LoB"
type.cfDNA$Pool = "LoBPoolA"
type.cfDNA$Pool[type.cfDNA$Donor %in% paste("Donor", 11:20, sep="_")] = "LoBPoolB"

# G360 neat gDNA
type.gDNA <- read.csv(file.path(G360Dir, "sample_list_neat_gDNA.csv"), header=T, stringsAsFactors=F)
colnames(type.gDNA)[ncol(type.gDNA)] = "Donor"
type.gDNA$run_sample_id = type.gDNA$SAMPLE_ID
type.gDNA = type.gDNA %>% filter(!grepl("Control", Sample_Description3))
type.gDNA$Donor = sapply(type.gDNA$Sample_Description3, function(x) gsub(" ", "_", x))
type.gDNA$Pool = "LoBPoolA"
type.gDNA$Pool[type.gDNA$Donor %in% paste("Donor", 11:20, sep="_")] = "LoBPoolB"

cols = intersect(colnames(type.gDNA), colnames(type.cfDNA))
type = rbind(type.cfDNA[,cols], type.gDNA[,cols])

# G360 Pools
type.pools = read.csv(file=file.path(G360Dir, "sample_list.csv"), header=T, stringsAsFactors=F)
type.pools$run_sample_id = type.pools$SAMPLE_ID
type.pools$Study.Name = "pools LoB"
type.pools$Pool = sapply(type.pools$Historical.Sample.ID.Notes, function(x) strsplit(x, "_")[[1]][1])
#type$runid = type$FLOWCELL_ID
type.pools$Donor = "pool"
type.pools = type.pools[,colnames(type)]
type = rbind(type, type.pools)

type = type %>% filter(Pool %in% c("LoBPoolA", "LoBPoolB"))

get_FolderIDs <- function(runid) {
  folderid = c()
  for (run_id in runid)
    folderid <- c(folderid,  type.s2f[which(type.s2f$runid == run_id), ]$Folder)
  
  return(folderid)
}

#########################
##  Load data from repo
##  bip outputs
#########################


if( !file.exists(file.path(resdir, "work.file.RData"))) {
  repo = g360.repo
  algoVersion = "v3.5"
  parameterSet <- fromJSON(parameterSetJson)
  probes = getCnvProbes(probeFile=probesBed, algoVersion = algoVersion)
  params = getCnvParams (jsonFile = parameterSetJson, algoVersion = algoVersion)
  
  tables = c("snv_call", "indel_call", "fusion_call", "cnv_call", "sample_coverage", "ghcnv_qc", "qc_on_target", "gh_sample")
  files.ext = c(".snv_call.hdr.tsv", ".indel_call.hdr.tsv", ".fusion_call.hdr.tsv", ".cnv_call.hdr.tsv", ".coverage.hdr.tsv", 
                ".ghcnv_qc.hdr.tsv", ".on_target_db.hdr.tsv", ".gh_sample_db.hdr.tsv")
  data = list()
  #repo = "/ghds/ivd/flowcentral"
  #repo = "/ghds/groups/BIP_release_testing/output/3.5.1/BIP3.5.1_rc4_av_LoD/pipeline"
  #repo = "/ghds/groups/BIP_release_testing/output/3.5.1/BIP3.5.1_rc4_av_lob/pipeline"
  run_sample_id = type$run_sample_id
  runid = type$runid
  
  folderid <- get_FolderIDs(runid)
  
  for(n in files.ext) data[[tables[which(n == files.ext)]]] = 
    import_data_from_repo(run_sample_id, folderid, repo='', params, probes, 
                          what = n, sep="\t", version = bip.version, mc.cores=8)
  names(data) = tables
  data = lapply(tables, function(n) data[[n]] = rbindlist(data[[n]]))
  names(data) = tables
  
  snv = data$snv_call
  snv$mut_key = paste(snv$gene, snv$mut_aa, sep = ".")
  # dedup flowcell
  snv$mut_key = paste(snv$gene, snv$mut_aa, sep = ".")
  snv$mut_key1 = paste(paste0("chr", snv$chrom), snv$position, snv$mut_nt, sep = ".")
  snv$mut_key2 = paste(snv$Pool, snv$chrom, snv$position, snv$mut_nt, sep="_")
  key = paste(snv$run_sample_id, snv$runid, snv$gene, snv$position, snv$mut_aa, sep="_")
  snv = snv[match(unique(key), key),]
  print(paste("Pulled ", nrow(snv), " SNVs for ", length(unique(snv$run_sample_id)), " samples"))
  
  indel = data$indel_call
  indel$mut_key = paste(indel$gene, indel$mut_aa, sep = ".")
  indel$mut_key1 = paste(paste0("chr", indel$chrom), indel$position, indel$mut_nt, sep = ".")
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
  anno = anno[,match(unique(colnames(anno)), colnames(anno))]
  
  save(snv, indel, fusion, cnv, nsc, autoqc, flqc, control, anno, 
       file=file.path(resdir, paste("work.file","RData", sep=".")))
}
load(file.path(resdir, "work.file.RData"))

type = type[match(anno$run_sample_id, type$run_sample_id),]
table(type$Donor, type$Study.Name)

#################
##  Ambry data
#################

if(!file.exists(file.path(resdir, "work.ref.RData"))) { #this bit needs to be commented out later
  ##  Ambry vcf files
  sample.info = fread(file=file.path(WESDir, "sample_list_ambry_180830.csv"), stringsAsFactors=F)
  sample.info$sample_id = sapply(sample.info$`Ambry Accessions`, function(x) gsub("-", "_", x))
  file_name <- function(f) strsplit(f, "/")[[1]][length(strsplit(f, "/")[[1]])] 
  vcf.files = file.path(WESDir, "vcf", dir(file.path(WESDir, "vcf")))
  sample_ids = sapply(vcf.files, function(f) gsub(".vcf", "", file_name(f)))
  snv.ref = mclapply(vcf.files, function(f) 
    VCF2df(readVcf(f), sample_id = sample_ids[which(f==vcf.files)], flowcell_id = "Ambry"), 
    mc.cores = 10)
  names(snv.ref) = sample_ids
  for(i in 1:length(snv.ref)) snv.ref[[i]]$mut_key2 = paste(sample.info$Pool[i], 
                                                            sapply(snv.ref[[i]]$chrom, function(x) gsub("chr", "", x)),
                                                            snv.ref[[i]]$position, snv.ref[[i]]$mut_nt, sep="_")
  for(i in 1:length(snv.ref)) snv.ref[[i]]$mut_key1 = paste(snv.ref[[i]]$chrom, 
                                                            snv.ref[[i]]$position, snv.ref[[i]]$mut_nt, sep=".")
  #snv.ref$mut_key = paste(snv.ref$gene, snv.ref$mut_aa, sep=".")
  
  # Ambry coverage files
  coverage.files = file.path(WESDir, "coverage", dir(file.path(WESDir, "coverage")))
  coverage.ref = mclapply(coverage.files, function(f) read.csv(file=f, skip=5, header=T, stringsAsFactors=F),
                          mc.cores = 10)
  targeted.regions.ref = coverage.ref[[1]] %>% dplyr::select(c("Gene", "Isoform", "Coding.exon.CDS.", "CDS.Length")) %>%
    mutate(exon = gsub("CDS", "", Coding.exon.CDS. ))
  save(snv.ref, coverage.ref, targeted.regions.ref, file=file.path(saveDir,  paste("work.ref", "RData", sep=".")))
}
load(file.path(resdir, "work.ref.RData"))


if( !file.exists(file.path(resdir, "work.MDL.RData") )) {
  # MDL data
  repo = LBP70Dir
  files = dir(repo)
  files = files[grepl("^SM", files)]
  type.MDL = fread(file.path(repo, "LOD cfDNA GH2MDL 10252018 Key.csv"))
  type.MDL$run_sample_id = sapply(type.MDL$`MDL Accn Number`, function(x) gsub("-", "", x))
  type.MDL = type.MDL %>% filter(!grepl("POOL", `GH Accn Number`))
  #type.MDL$Donor = type$Donor[match(type.MDL$`GH Accn Number`, type$Historical.Sample.ID.Notes)]
  type.MDL$Donor = type$Donor[sapply(type.MDL$`GH Accn Number`, function(x) 
    which(x == type$'Historical.Sample.ID.Notes' & 
            type$Study.Name == "cfDNA LoB"))]
  
  snv.MDL = lapply(type.MDL$run_sample_id, function(n) {f = files[grepl(n, files) & grepl("snv_call.hdr.tsv", files)];fread(file.path(repo, f))})
  snv.MDL = do.call(rbind, snv.MDL)
  snv.MDL$mut_key = paste(snv.MDL$gene, snv.MDL$mut_aa, sep = ".")
  indel.MDL = lapply(type.MDL$run_sample_id, function(n) {f = files[grepl(n, files) & grepl("indel_call.hdr.tsv", files)];fread(file.path(repo, f))})
  indel.MDL = do.call(rbind, indel.MDL)
  indel.MDL$mut_key = paste(indel.MDL$gene, indel.MDL$mut_aa, sep = ".")
  cnv.MDL = lapply(type.MDL$run_sample_id, function(n) {f = files[grepl(n, files) & grepl("cnv_call.hdr.tsv", files)];fread(file.path(repo, f))})
  cnv.MDL = do.call(rbind, cnv.MDL)
  cnv.MDL$mut_key = paste(cnv.MDL$gene, "CNV", sep = ".")
  fusion.MDL = lapply(type.MDL$run_sample_id, function(n) {f = files[grepl(n, files) & grepl("fusion_call.hdr.tsv", files)];fread(file.path(repo, f))})
  fusion.MDL = do.call(rbind, fusion.MDL)
  fusion.MDL$mut_key = paste(fusion.MDL$gene_a, fusion.MDL$pos_a, sep = ".")
  ghcnv_qc.MDL = lapply(type.MDL$run_sample_id, function(n) {f = files[grepl(n, files) & grepl("ghcnv_qc.hdr.tsv", files)];fread(file.path(repo, f))})
  ghcnv_qc.MDL = do.call(rbind, ghcnv_qc.MDL)
  gh_board.MDL = lapply(type.MDL$run_sample_id, function(n) {f = files[grepl(n, files) & grepl("gh_board.hdr.tsv", files)];fread(file.path(repo, f))})
  gh_board.MDL = do.call(rbind, gh_board.MDL)
  nsc.MDL = lapply(type.MDL$run_sample_id, function(n) {f = files[grepl(n, files) & grepl("coverage.hdr.tsv", files)];fread(file.path(repo, f))})
  nsc.MDL = do.call(rbind, nsc.MDL)
  sample_qc.MDL = fread(file.path(repo, "autoqc_sample_qc.hdr.tsv"))
  anno.MDL = ghcnv_qc.MDL
  save(type.MDL, snv.MDL, indel.MDL, cnv.MDL, fusion.MDL, nsc.MDL, anno.MDL, sample_qc.MDL, 
       file=file.path(saveDir, paste("work.MDL", "RData", sep=".")))
}
load(file.path(resdir, "work.MDL.RData"))

###############################
##  Table 1: QC metrics
###############################

##  Appendix 1.A: sequencing QC metrics
runid = unique(type$runid[type$Study.Name == "pools LoB"])
m = "flowcell_cluster_density"
zz = flqc %>% filter(metric == m)
ClusterDensity = zz$value[match(runid, zz$runid)]
m = "flowcell_clusters_passing_filter"
zz = flqc %>% filter(metric == m)
ClusterPassingFilter = zz$value[match(runid, zz$runid)]
metrics = c("flowcell_phasing_1", "flowcell_phasing_2", "flowcell_prephasing_1", "flowcell_prephasing_2")
zz = sapply(metrics, function(m) {zz = flqc %>% filter(metric == m);zz$value[match(runid, zz$runid)]})
PhasingAndPrePhasingReads1And2 = apply(zz, 1, function(x) ifelse(all(x < 0.01), "<0.01%", ">0.01%") )
m = "flowcell_qscore_1"
zz = flqc %>% filter(metric == m)
QScoreRead1 = zz$value[match(runid, zz$runid)]
m = "flowcell_qscore_2"
zz = flqc %>% filter(metric == m)
QScoreRead2 = zz$value[match(runid, zz$runid)]

tab = data.frame(runid = runid,
                 ClusterDensity = ClusterDensity,
                 ClusterPassingFilter = ClusterPassingFilter,
                 PhasingAndPrePhasingReads1And2 = PhasingAndPrePhasingReads1And2,
                 QScoreRead1 = QScoreRead1,
                 QScoreRead2 = QScoreRead2)
tab
write.table(tab, file=file.path(resdir, "Appendix1.A.tsv"), sep="\t", col.names=T, row.names=T, quote=F )


summaryQC = table(autoqc$metric, autoqc$status)
summaryQC
write.table(summaryQC, file=file.path(resdir, "Appendix1.B.tsv"), sep="\t", col.names=T, row.names=T, quote=F)

#####################
##  Control samples
#####################

##  controls nsc
tables = c( "sample_coverage")
files.ext = c(".coverage.hdr.tsv")
data.control = list()
repo = "/ghds/ivd/flowcentral" #get control results from the main flowcentral repo
#repo = "/ghds/groups/BIP_release_testing/output/3.5.1/BIP3.5.1_rc4_av_lob/pipeline" ##these folders do not have BIP results
run_sample_id = unique(control$run_sample_id)
runid = control$runid[match(run_sample_id, control$run_sample_id)]

folderid <- get_FolderIDs(runid)

for(n in files.ext) data.control[[tables[which(n == files.ext)]]] = 
  import_data_from_repo(run_sample_id, folderid, repo='', params, probes, 
                        what = n, sep="\t", version = "3.5", mc.cores=8)
names(data.control) = tables
data.control = lapply(tables, function(n) data.control[[n]] = rbindlist(data.control[[n]]))
names(data.control) = tables
nsc.control = data.control$sample_coverage
key = paste(nsc.control$run_sample_id, 
            nsc.control$runid, 
            nsc.control$chrom, 
            nsc.control$start_pos, 
            nsc.control$end_pos, sep="_")
nsc.control = nsc.control[match(unique(key), key),]
print(paste("Pulled ", nrow(nsc.control), " NSCs for ", length(unique(nsc.control$run_sample_id)), " samples"))

##  Appendix 2: CONTROL PERFORMANCE CHARACTERISTICS
run_sample_id = unique(control$run_sample_id)
runid = control$runid[match(run_sample_id, control$run_sample_id)]

m = "sample_contamination_pct"
zz = control %>% filter(metric == m)
FamilyContamination = zz$value[match(runid, zz$runid)]

NonSingletonCoverage = tapply(nsc.control$median_coverage, nsc.control$run_sample_id, median)
NonSingletonCoverage = NonSingletonCoverage[match(run_sample_id, names(NonSingletonCoverage))]

metrics = c("aiocontrol_sensitivity_snv", 
            "aiocontrol_sensitivity_indel", 
            "aiocontrol_sensitivity_cnv", 
            "aiocontrol_sensitivity_fusion")
zz = sapply(metrics, function(m) {zz = control %>% filter(metric == m);zz$value[match(runid, zz$runid)]})
SensitivityForTargetedVariants = apply(zz, 1, function(x) ifelse(all(x == 100), "100%", "<100%") )

metrics = c("aiocontrol_false_positive_indel", "aiocontrol_false_positive_snv")
zz = sapply(metrics, function(m) {zz = control %>% filter(metric == m);zz$value[match(runid, zz$runid)]})
SpecificityForTargetedVariants = apply(zz, 1, function(x) ifelse(all(x == 0), "100%", "<100%") )

tab = data.frame(run_sample_id = run_sample_id,
                 runid = runid,
                 NonSingletonCoverage = NonSingletonCoverage,
                 FamilyContamination = FamilyContamination,
                 SensitivityForTargetedVariants = SensitivityForTargetedVariants,
                 SpecificityForTargetedVariants = SpecificityForTargetedVariants)
tab
write.table(tab, file=file.path(resdir, "Appendix2.tsv"), sep="\t", col.names=T, row.names=F, quote=F )

write.table(control,  file=file.path(resdir, "Table1.Control.results.tsv"), 
            sep="\t", col.names=T, row.names=F, quote=F)
controlQC = table(control$metric, control$status)
controlQC
write.table(controlQC, file=file.path(resdir, "Table1.Control.summary.tsv"), 
            sep="\t", col.names=T, row.names=F, quote=F)

#########################  
# Remove failed sample
#########################

metrics = c("sample_contamination_pct", 
            "sample_coverage_exceptions", 
            "sample_gc_bias", 
            "sample_non_singleton_families", 
            "sample_on_target_rate")
failedDF = autoqc %>% filter(metric %in% metrics & status == "FAIL")
write.table(failedDF, file=file.path(resdir, "failedDF.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
failedSamples =failedDF$run_sample_id

df.failed = autoqc %>% filter(run_sample_id == failedSamples)
samples2use = setdiff(unique(type$run_sample_id), failedSamples)



##  filtered data set  #here its removing the sample that failed to meet the QC metrics
snv = snv %>% filter(run_sample_id %in% samples2use)
indel = indel %>% filter(run_sample_id %in% samples2use)
cnv = cnv %>% filter(run_sample_id %in% samples2use)
fusion = fusion %>% filter(run_sample_id %in% samples2use)
nsc = nsc %>% filter(run_sample_id %in% samples2use)
autoqc = autoqc %>% filter(run_sample_id %in% samples2use)
anno = anno %>% filter(run_sample_id %in% samples2use)
type = type %>% filter(run_sample_id %in% samples2use)

#table(type$Donor, type$Study.Name)
table(type$Pool, type$Study.Name)


##############
##  FPs list
##############

#load(file.path(resdir, "work.ref.2019-01-22.RData"))

##  Map Isoform names between Ambry and GH
#table(targeted.regions.ref$Isoform[targeted.regions.ref$Gene == "NTRK1"])
#table(snv2$transcript_id[snv2$gene == "NTRK1"])
targeted.regions.ref$Isoform[targeted.regions.ref$Gene == "NTRK1"] = "NM_002529"

ref.targeted.keys = paste(targeted.regions.ref$Gene, targeted.regions.ref$Isoform, targeted.regions.ref$exon, sep=".")

##  SNVs
snv2 = snv %>% left_join(type, by=c("run_sample_id", "runid")) #joins the two tables by run_sample_id and runid columns
#snv2 are variants identified by BIP

snv.ref.df = lapply(snv.ref, function(d) d %>% filter(mut_key1 %in% snv2$mut_key1)) #snv.ref is from Ambry data set
snv.ref.df = do.call(rbind, snv.ref.df)

snv.pool.all = snv2 %>% filter(run_sample_id %in% type$run_sample_id[type$Study.Name == "pools LoB"] & #discarding all variants with these entries in the variant comment column 
                                 # filter gremline, these are never going to be somatic calls.
                                 !grepl("not reportable|ExAC common germline|GH SNP common germline", variant_comment)) #this comes from different part of the pipeline and not from GH-somatic
#these variants are not FPs because they never end up being reported to the patients.
#now what remains are the candidate variants that we need to investigate
#snv.pool.all is all the remaining variants that we need to investigate for FP
vars = unique(snv.pool.all$mut_key1)
df.snv.pool.all =  data.frame(chrom = snv.pool.all$chrom[match(vars, snv.pool.all$mut_key1)],
                              position = snv.pool.all$position[match(vars, snv.pool.all$mut_key1)],
                              gene = snv.pool.all$gene[match(vars, snv.pool.all$mut_key1)],
                              transcript_id = snv.pool.all$transcript_id[match(vars, snv.pool.all$mut_key1)],
                              exon = snv.pool.all$exon[match(vars, snv.pool.all$mut_key1)],
                              mut_nt = snv.pool.all$mut_nt[match(vars, snv.pool.all$mut_key1)],
                              mut_aa = snv.pool.all$mut_aa[match(vars, snv.pool.all$mut_key1)], 
                              mut_key = paste(snv.pool.all$gene[match(vars, snv.pool.all$mut_key1)], 
                                              snv.pool.all$mut_aa[match(vars, snv.pool.all$mut_key1)], sep="."),
                              mut_key1 = paste0("chr", snv.pool.all$chrom[match(vars, snv.pool.all$mut_key1)], 
                                                ".", snv.pool.all$position[match(vars, snv.pool.all$mut_key1)],
                                                ".", snv.pool.all$mut_nt[match(vars, snv.pool.all$mut_key1)]),
                              variant_type = "SNV", 
                              Pool = snv.pool.all$Pool[match(vars, snv.pool.all$mut_key1)], 
                              run_sample_ids = sapply(vars, function(v) paste(snv.pool.all$run_sample_id[v == snv.pool.all$mut_key1], collapse=",")),
                              percentage =  snv.pool.all$percentage[match(vars, snv.pool.all$mut_key1)],
                              mut_cnt =  round(snv.pool.all$percentage[match(vars, snv.pool.all$mut_key1)] * snv.pool.all$mol_cnt[match(vars, snv.pool.all$mut_key1)]/100),
                              stringsAsFactors=F) %>%
  mutate(nObs.pool = sapply(vars, function(v) sum(v == snv2$mut_key1[snv2$Study.Name == "pools LoB"], na.rm=T))) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(snv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.call = vars %in% snv2$mut_key1[snv2$Study.Name == "cfDNA LoB"]) %>% #tries to see of the variant is observed in cf-DNA
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(snv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.donor = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, snv2$Donor[ii], NA)})) %>%
  mutate(gDNA.call = vars %in% snv2$mut_key1[snv2$Study.Name == "gDNA LoB"]) %>%
  mutate(gDNA.maf = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, mean(snv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(gDNA.donor = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, snv2$Donor[ii], NA)})) %>%  # Check if variant on Ambry panel
  mutate(isoform = sapply(transcript_id, function(x) strsplit(x, "\\.")[[1]][1])) %>%
  mutate(var.targeted.key = paste(gene, isoform, exon, sep=".")) %>%
  mutate(Ambry.target = var.targeted.key %in%  ref.targeted.keys) %>% #replace the Ambry data by MDL variants
  arrange(as.numeric(chrom), position)

df.snv.pool.all$Ambry.maf = sapply(df.snv.pool.all$mut_key1, function(mut) {ii = which(snv.ref.df$mut_key1 == mut);
ifelse(length(ii) > 0, as.numeric(snv.ref.df$AF[ii[1]]), NA)})
df.snv.pool.all$MDL.call = df.snv.pool.all$mut_key %in% snv.MDL$mut_key
#df.snv.pool.all$MDL.maf = sapply(df.snv.pool.all$mut_key, function(mut) snv.MDL$percentage[mut == snv.MDL$mut_key])
df.snv.pool.all$MDL.maf = snv.MDL$percentage[match(df.snv.pool.all$mut_key, snv.MDL$mut_key)]/100
df.snv.pool.all$MDL.donor = sapply(df.snv.pool.all$mut_key, function(mut) 
{ii = which(snv.MDL$mut_key == mut);
out = NA;
if(length(ii) > 0) out = type.MDL$Donor[type.MDL$run_sample_id == snv.MDL$run_sample_id[ii]];
out})
#make an extra column reporting FP or TP
df.snv.pool.all
write.table(df.snv.pool.all, file=file.path(resdir, "Table2.All.SNV.tsv"), 
            sep="\t", col.names=T, row.names=F, quote=F)

snv.pool = snv2 %>% filter(run_sample_id %in% type$run_sample_id[type$Study.Name == "pools LoB"] &
                             # filter gremline
                             !grepl("not reportable|ExAC common germline|GH SNP common germline", variant_comment) &
                             # filter vars observed in Ambry
                             !(mut_key1 %in% snv.ref.df$mut_key1))

vars = unique(snv.pool$mut_key1)
df.snv.variant.pool =  data.frame(chrom = snv.pool$chrom[match(vars, snv.pool$mut_key1)],
                                  position = snv.pool$position[match(vars, snv.pool$mut_key1)],
                                  gene = snv.pool$gene[match(vars, snv.pool$mut_key1)],
                                  transcript_id = snv.pool$transcript_id[match(vars, snv.pool$mut_key1)],
                                  exon = snv.pool$exon[match(vars, snv.pool$mut_key1)],
                                  mut_nt = snv.pool$mut_nt[match(vars, snv.pool$mut_key1)],
                                  mut_aa = snv.pool$mut_aa[match(vars, snv.pool$mut_key1)], 
                                  mut_key = paste(snv.pool.all$gene[match(vars, snv.pool.all$mut_key1)], 
                                                  snv.pool.all$mut_aa[match(vars, snv.pool.all$mut_key1)], sep="."),
                                  mut_key1 = paste0("chr", snv.pool.all$chrom[match(vars, snv.pool.all$mut_key1)], 
                                                    ".", snv.pool.all$position[match(vars, snv.pool.all$mut_key1)],
                                                    ".", snv.pool.all$mut_nt[match(vars, snv.pool.all$mut_key1)]),
                                  variant_type = "SNV", 
                                  Pool = snv.pool$Pool[match(vars, snv.pool$mut_key1)], 
                                  run_sample_ids = sapply(vars, function(v) paste(snv.pool$run_sample_id[v == snv.pool$mut_key1], collapse=",")),
                                  stringsAsFactors=F) %>%
  mutate(nObs.pool = sapply(vars, function(v) sum(v == snv2$mut_key1[snv2$Study.Name == "pools LoB"], na.rm=T))) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(snv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.call = vars %in% snv2$mut_key1[snv2$Study.Name == "cfDNA LoB"]) %>%
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(snv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.donor =  sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, snv2$Donor[ii], NA)})) %>%
  mutate(gDNA.call = vars %in% snv2$mut_key1[snv2$Study.Name == "gDNA LoB"]) %>%
  mutate(gDNA.maf = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, mean(snv2$percentage[ii]/100, na.rm=T), NA)})) %>% # Check if variant on Ambry panel
  mutate(gDNA.donor = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, snv2$Donor[ii], NA)})) %>%
  mutate(isoform = sapply(transcript_id, function(x) strsplit(x, "\\.")[[1]][1])) %>%
  mutate(var.targeted.key = paste(gene, isoform, exon, sep=".")) %>%
  mutate(Ambry.target = var.targeted.key %in%  ref.targeted.keys) %>%
  filter(Ambry.target) %>%
  arrange(as.numeric(chrom), position)

df.snv.variant.pool$Ambry.maf = sapply(df.snv.variant.pool$mut_key1, function(mut) {ii = which(snv.ref.df$mut_key1 == mut);
ifelse(length(ii) > 0, as.numeric(snv.ref.df$AF[ii[1]]), NA)})
df.snv.variant.pool$MDL.call = df.snv.variant.pool$mut_key %in% snv.MDL$mut_key
#df.snv.variant.pool$MDL.maf = sapply(df.snv.variant.pool$mut_key, function(mut) snv.MDL$percentage[mut == snv.MDL$mut_key])
df.snv.variant.pool$MDL.maf = snv.MDL$percentage[match(df.snv.variant.pool$mut_key, snv.MDL$mut_key)]/100
df.snv.variant.pool$MDL.donor = sapply(df.snv.variant.pool$mut_key, function(mut) 
{ii = which(snv.MDL$mut_key == mut);
out = NA;
if(length(ii) > 0) out = type.MDL$Donor[type.MDL$run_sample_id == snv.MDL$run_sample_id[ii]];
out})

df.snv.variant.pool

##  APPENDIX 3 TABLE
df.snv.variants = df.snv.variant.pool %>% dplyr::select(variant_type, gene, chrom, position, mut_nt, mut_aa) %>%
  mutate(ClinicalllySig = "No", nObs = df.snv.variant.pool$nObs.pool) %>%
  mutate(nTotal = sapply(df.snv.variant.pool$Pool, function(p) sum(type$Pool == p & type$Study.Name == "pools LoB"))) %>%
  mutate(FractionPos = nObs/nTotal) %>%
  mutate(MAF = df.snv.variant.pool$pool.maf)

##  APPENDIX 2 TABLE
#vars = snv.pool$mut_key1[snv.pool$mut_key1 %in% vars]
#iVars = which(snv.pool$mut_key1 %in% vars)
zz = snv.pool.all %>% mutate(isoform = sapply(transcript_id, function(x) strsplit(x, "\\.")[[1]][1])) %>%
  mutate(ClinicallySig = "No") %>%  ##  MANUALLY CHECKED THIS ENTRY
  mutate(Ambry.target.key = paste(gene, isoform, exon, sep=".")) %>%
  mutate(Ambry.target = Ambry.target.key %in%  ref.targeted.keys) 
df.snv.sample =  zz %>% dplyr::select( run_sample_id,gene,chrom,position,mut_nt,mut_aa, ClinicallySig) %>%
  mutate(MAF = zz$percentage/100) %>%
  mutate(Confirmed = zz$mut_key1 %in% snv.ref.df$mut_key1 |!zz$Ambry.target) %>%
  #mutate(Confirmed = sapply(mut_key1, function(x) ifelse(any(x == snv.ref.df$mut_key1), "Yes", "No"))) %>%
  #mutate(Total_count = mol_cnt) %>%
  mutate(mol_cnt = zz$mol_cnt) %>%
  mutate(mut_cnt = round(mol_cnt * MAF))

head(df.snv.sample)


##  INDELS
indel2 = indel %>% left_join(type, by=c("run_sample_id", "runid"))

indel.ref.df = lapply(snv.ref, function(d) d %>% filter(mut_key1 %in% indel2$mut_key1))
indel.ref.df = do.call(rbind, indel.ref.df)

ids = type$run_sample_id[type$Study.Name == "pools LoB"]

##  INDELS all variants
indel.pool.all = indel2 %>% filter(run_sample_id %in% ids &
                                     # filter gremline
                                     !grepl("not reportable|ExAC common germline|GH SNP common germline", variant_comment))

vars = unique(indel.pool.all$mut_key1)

df.indel.pool.all =  data.frame(chrom = indel.pool.all$chrom[match(vars, indel.pool.all$mut_key1)],
                                position = indel.pool.all$position[match(vars, indel.pool.all$mut_key1)],
                                gene = indel.pool.all$gene[match(vars, indel.pool.all$mut_key1)],
                                transcript_id = indel.pool.all$transcript_id[match(vars, indel.pool.all$mut_key1)],
                                exon = indel.pool.all$exon[match(vars, indel.pool.all$mut_key1)],
                                mut_nt = indel.pool.all$mut_nt[match(vars, indel.pool.all$mut_key1)],
                                mut_aa = indel.pool.all$mut_aa[match(vars, indel.pool.all$mut_key1)], 
                                mut_key = paste(indel.pool.all$gene[match(vars, indel.pool.all$mut_key1)], 
                                                indel.pool.all$mut_aa[match(vars, indel.pool.all$mut_key1)], sep="."),
                                mut_key1 = paste0("chr", indel.pool.all$chrom[match(vars, indel.pool.all$mut_key1)], 
                                                  ".", indel.pool.all$position[match(vars, indel.pool.all$mut_key1)],
                                                  ".", indel.pool.all$mut_nt[match(vars, indel.pool.all$mut_key1)]),
                                variant_type = "INDEL", 
                                Pool = indel.pool.all$Pool[match(vars, indel.pool.all$mut_key1)], 
                                run_sample_ids = sapply(vars, function(v) paste(indel.pool.all$run_sample_id[v == indel.pool.all$mut_key1], collapse=",")),
                                stringsAsFactors=F) %>%
  mutate(nObs.pool = sapply(vars, function(v) sum(v == indel2$mut_key1[indel2$Study.Name == "pools LoB"], na.rm=T))) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(indel2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.call = vars %in% indel2$mut_key1[indel2$Study.Name == "cfDNA LoB"]) %>%
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(indel2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.donor = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, indel2$Donor[ii], NA)})) %>%
  mutate(gDNA.call = vars %in% indel2$mut_key1[indel2$Study.Name == "gDNA LoB"]) %>%
  mutate(gDNA.maf = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, mean(indel2$percentage[ii]/100, na.rm=T), NA)})) %>% # Check if variant on Ambry panel
  mutate(gDNA.donor = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, indel2$Donor[ii], NA)})) %>%
  mutate(isoform = sapply(transcript_id, function(x) strsplit(x, "\\.")[[1]][1])) %>%
  mutate(var.targeted.key = paste(gene, isoform, exon, sep=".")) %>%
  mutate(Ambry.target = var.targeted.key %in%  ref.targeted.keys) %>%
  arrange(as.numeric(chrom), position)

df.indel.pool.all$Ambry.maf = sapply(df.indel.pool.all$mut_key1, function(mut) {ii = which(indel.ref.df$mut_key1 == mut);
ifelse(length(ii) > 0, as.numeric(indel.ref.df$AF[ii[1]]), NA)})
df.indel.pool.all$MDL.call = df.indel.pool.all$mut_key %in% indel.MDL$mut_key
#df.indel.pool.all$MDL.maf = sapply(df.indel.pool.all$mut_key, function(mut) indel.MDL$percentage[mut == indel.MDL$mut_key])
df.indel.pool.all$MDL.maf = indel.MDL$percentage[match(df.indel.pool.all$mut_key, indel.MDL$mut_key)]/100
df.indel.pool.all$MDL.donor = sapply(df.indel.pool.all$mut_key, function(mut) 
{ii = which(indel.MDL$mut_key == mut);
out = NA;
if(length(ii) > 0) out = type.MDL$Donor[type.MDL$run_sample_id == snv.MDL$run_sample_id[ii]];
out})

df.indel.pool.all
write.table(df.indel.pool.all, file=file.path(resdir, "Table2.All.INDEL.tsv"), 
            sep="\t", col.names=T, row.names=F, quote=F)

##  Filtered INDELs
indel.pool = indel2 %>% filter(run_sample_id %in% ids &
                                 # filter gremline
                                 !grepl("not reportable|ExAC common germline|GH SNP common germline", variant_comment) &
                                 # filter vars observed in Ambry
                                 !(mut_key1 %in% indel.ref.df$mut_key1))

vars = unique(indel.pool$mut_key1)

df.indel.variant.pool =  data.frame(chrom = indel.pool$chrom[match(vars, indel.pool$mut_key1)],
                                    position = indel.pool$position[match(vars, indel.pool$mut_key1)],
                                    gene = indel.pool$gene[match(vars, indel.pool$mut_key1)],
                                    transcript_id = indel.pool$transcript_id[match(vars, indel.pool$mut_key1)],
                                    exon = indel.pool$exon[match(vars, indel.pool$mut_key1)],
                                    mut_nt = indel.pool$mut_nt[match(vars, indel.pool$mut_key1)],
                                    mut_aa = indel.pool$mut_aa[match(vars, indel.pool$mut_key1)], 
                                    mut_key = paste(indel.pool.all$gene[match(vars, indel.pool.all$mut_key1)], 
                                                    indel.pool.all$mut_aa[match(vars, indel.pool.all$mut_key1)], sep="."),
                                    mut_key1 = paste0("chr", indel.pool.all$chrom[match(vars, indel.pool.all$mut_key1)], 
                                                      ".", indel.pool.all$position[match(vars, indel.pool.all$mut_key1)],
                                                      ".", indel.pool.all$mut_nt[match(vars, indel.pool.all$mut_key1)]),
                                    variant_type = "INDEL", 
                                    Pool = indel.pool$Pool[match(vars, indel.pool$mut_key1)], 
                                    run_sample_ids = sapply(vars, function(v) paste(indel.pool$run_sample_id[v == indel.pool$mut_key1], collapse=",")),
                                    stringsAsFactors=F) %>%
  mutate(nObs.pool = sapply(vars, function(v) sum(v == indel2$mut_key1[indel2$Study.Name == "pools LoB"], na.rm=T))) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(indel2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.call = vars %in% indel2$mut_key1[indel2$Study.Name == "cfDNA LoB"]) %>%
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(indel2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.donor = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, indel2$Donor[ii], NA)})) %>%
  mutate(gDNA.call = vars %in% indel2$mut_key1[indel2$Study.Name == "gDNA LoB"]) %>%
  mutate(gDNA.maf = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, mean(indel2$percentage[ii]/100, na.rm=T), NA)})) %>% # Check if variant on Ambry panel
  mutate(gDNA.donor = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, indel2$Donor[ii], NA)})) %>%
  mutate(isoform = sapply(transcript_id, function(x) strsplit(x, "\\.")[[1]][1])) %>%
  mutate(var.targeted.key = paste(gene, isoform, exon, sep=".")) %>%
  mutate(Ambry.target = var.targeted.key %in%  ref.targeted.keys) %>%
  filter(Ambry.target)

df.indel.variant.pool$Ambry.maf = sapply(df.indel.variant.pool$mut_key1, function(mut) {ii = which(indel.ref.df$mut_key1 == mut);
ifelse(length(ii) > 0, as.numeric(indel.ref.df$AF[ii[1]]), NA)})
df.indel.variant.pool$MDL.call = df.indel.variant.pool$mut_key %in% indel.MDL$mut_key
#df.snv.pool.all$MDL.maf = sapply(df.snv.pool.all$mut_key, function(mut) snv.MDL$percentage[mut == snv.MDL$mut_key])
df.indel.variant.pool$MDL.maf = indel.MDL$percentage[match(df.indel.variant.pool$mut_key, indel.MDL$mut_key)]/100
df.indel.variant.pool$MDL.donor = sapply(df.indel.variant.pool$mut_key, function(mut) 
{ii = which(indel.MDL$mut_key == mut);
out = NA;
if(length(ii) > 0) out = type.MDL$Donor[type.MDL$run_sample_id == indel.MDL$run_sample_id[ii]];
out})

df.indel.variant.pool

##  APPENDIX 3 TABLE
df.indel.variants = df.indel.variant.pool %>% dplyr::select(variant_type, gene, chrom, position, mut_nt, mut_aa) %>%
  mutate(ClinicalllySig = "No", nObs = df.indel.variant.pool$nObs.pool) %>%  ## CHECK MANUALLY THE CASE 
  mutate(nTotal = sapply(df.indel.variant.pool$Pool, function(p) sum(type$Pool == p & type$Study.Name == "pools LoB"))) %>%
  mutate(FractionPos = nObs/nTotal) %>%
  mutate(MAF = df.indel.variant.pool$pool.maf)
df.indel.variants 

##  APPENDIX 2 TABLE
zz = indel.pool.all %>% mutate(isoform = sapply(transcript_id, function(x) strsplit(x, "\\.")[[1]][1])) %>%
  mutate(ClinicallySig = "No") %>%  ##  MANUALLY CHECKED THIS ENTRY
  mutate(Ambry.target.key = paste(gene, isoform, exon, sep=".")) %>%
  mutate(Ambry.target = Ambry.target.key %in%  ref.targeted.keys) %>%
  mutate(mol_cnt = fam_cnt) %>%
  mutate(mut_cnt = indel_fam_cnt)
df.indel.sample =  zz %>% dplyr::select( run_sample_id,gene,chrom,position,mut_nt,mut_aa) %>% #, percentage, mol_cnt, mut_key1, Ambry.target) %>%
  mutate(ClinicallySig = zz$ClinicallySig) %>%
  mutate(MAF = zz$percentage/100) %>%
  mutate(Confirmed = zz$mut_key1 %in% indel.ref.df$mut_key1 |!zz$Ambry.target) %>%
  #mutate(Confirmed = sapply(mut_key1, function(x) ifelse(any(x == indel.ref.df$mut_key1), "Yes", "No"))) %>%
  #mutate(Total_count = mol_cnt) %>%
  mutate(mol_cnt = zz$mol_cnt) %>%
  mutate(mut_cnt = zz$mut_cnt)
head(df.indel.sample)

##  CNVs

cnv2 = cnv %>% left_join(type, by=c("run_sample_id", "runid"))
cnv2$mut_key2 = paste(cnv2$Pool, cnv2$gene, "CNV", sep=".")
cnv.ref.df = cnv2

##  CNVS all variants
cnv.pool.all = cnv2 %>% filter(run_sample_id %in% type$run_sample_id[type$Study.Name == "pools LoB"] &
                                 call == 1)
cat(paste("CNVs =", nrow(cnv.pool.all), "...\n"))
vars = unique(cnv.pool.all$mut_key2)

df.cnv.pool.all =  data.frame(chrom = cnv.pool.all$chrom[match(vars, cnv.pool.all$mut_key2)],
                              position = rep("", length(vars)),
                              gene = cnv.pool.all$gene[match(vars, cnv.pool.all$mut_key2)],
                              transcript_id = rep("", length(vars)),
                              exon = rep("", length(vars)),
                              mut_nt = rep("", length(vars)),
                              mut_aa = rep("", length(vars)), 
                              mut_key = paste(cnv.pool.all$gene[match(vars, cnv.pool.all$mut_key2)], 
                                              cnv.pool.all$mut_aa[match(vars, cnv.pool.all$mut_key2)], sep="."),
                              mut_key1 = paste0(rep("chr", length(vars)), cnv.pool.all$chrom[match(vars, cnv.pool.all$mut_key2)], 
                                                rep(".", length(vars)), cnv.pool.all$position[match(vars, cnv.pool.all$mut_key2)],
                                                rep(".", length(vars)), cnv.pool.all$mut_nt[match(vars, cnv.pool.all$mut_key2)]),
                              variant_type = rep("CNV", length(vars)), 
                              Pool = cnv.pool.all$Pool[match(vars, cnv.pool.all$mut_key2)], 
                              run_sample_ids = rep("", length(vars)),
                              percentage =  cnv.pool.all$copy_number[match(vars, snv.pool.all$mut_key1)],
                              mut_cnt =  rep("", length(vars)),
                              stringsAsFactors=F) %>%
  mutate(nObs.pool = sapply(vars, function(v) sum(v == cnv2$mut_key2[cnv2$Study.Name == "pools LoB"], na.rm=T))) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == cnv2$mut_key2 & cnv2$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(cnv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.call = vars %in% cnv2$mut_key2[cnv2$Study.Name == "cfDNA LoB"]) %>%
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == cnv2$mut_key2 & cnv2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(cnv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.donor = sapply(vars, function(v) {ii = which(v == cnv2$mut_key2 & cnv2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, cnv2$Donor[ii], NA)})) %>%
  mutate(gDNA.call = vars %in% cnv2$mut_key2[cnv2$Study.Name == "gDNA LoB"]) %>%
  mutate(gDNA.maf = sapply(vars, function(v) {ii = which(v == cnv2$mut_key2 & cnv2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, mean(cnv2$percentage[ii]/100, na.rm=T), NA)})) %>% # Check if variant on Ambry panel
  mutate(gDNA.donor = sapply(vars, function(v) {ii = which(v == cnv2$mut_key2 & cnv2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, cnv2$Donor[ii], NA)})) %>%
  mutate(isoform = sapply(transcript_id, function(x) strsplit(x, "\\.")[[1]][1])) %>%
  mutate(var.targeted.key = paste(gene, isoform, exon, sep=".")) %>%
  mutate(Ambry.target = var.targeted.key %in%  ref.targeted.keys) %>%
  arrange(as.numeric(chrom), position)

df.cnv.pool.all$Ambry.maf = sapply(df.cnv.pool.all$mut_key2, function(mut) {ii = which(cnv.ref.df$mut_key2 == mut);
ifelse(length(ii) > 0, as.numeric(cnv.ref.df$AF[ii[1]]), NA)})
df.cnv.pool.all$MDL.call = df.cnv.pool.all$mut_key %in% cnv.MDL$mut_key
#df.cnv.pool.all$MDL.maf = sapply(df.cnv.pool.all$mut_key, function(mut) cnv.MDL$percentage[mut == cnv.MDL$mut_key])
df.cnv.pool.all$MDL.maf = cnv.MDL$percentage[match(df.cnv.pool.all$mut_key, cnv.MDL$mut_key)]/100
df.cnv.pool.all$MDL.donor = sapply(df.cnv.pool.all$mut_key, function(mut) 
{ii = which(cnv.MDL$mut_key == mut);
out = NA;
if(length(ii) > 0) out = type.MDL$Donor[type.MDL$run_sample_id == snv.MDL$run_sample_id[ii]];
out})

df.cnv.pool.all
write.table(df.cnv.pool.all, file=file.path(resdir, "Table2.All.CNV.tsv"), 
            sep="\t", col.names=T, row.names=F, quote=F)



##  Filtered CNVs
cnv.pool = cnv2 %>% filter(run_sample_id %in% ids )

vars = unique(cnv.pool$mut_key1)

df.cnv.variant.pool =  data.frame(chrom = cnv.pool$chrom[match(vars, cnv.pool$mut_key2)],
                                  position = rep(NA, length(vars)),
                                  gene = cnv.pool$gene[match(vars, cnv.pool$mut_key2)],
                                  transcript_id = rep(NA, length(vars)),
                                  exon = rep(NA, length(vars)),
                                  mut_nt = rep(NA, length(vars)),
                                  mut_aa = rep(NA, length(vars)),
                                  mut_key = cnv.pool$mut_key[match(vars, cnv.pool$mut_key2)], 
                                  mut_key1 = cnv.pool$mut_key2[match(vars, cnv.pool$mut_key2)], 
                                  stringsAsFactors=F) %>%
  mutate(variant_type = rep("CNV", length(vars))) %>% 
  mutate(Pool = cnv.pool$Pool[match(vars, cnv.pool$mut_key2)]) %>%
  mutate(run_sample_ids = sapply(vars, function(v) paste(cnv.pool$run_sample_id[v == cnv.pool$mut_key2], collapse=","))) %>%
  mutate(nObs.pool = sapply(vars, function(v) sum(v == cnv2$mut_key2[cnv2$Study.Name == "pools LoB"], na.rm=T))) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == cnv2$mut_key2 & cnv2$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(cnv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.call = vars %in% cnv2$mut_key2[cnv2$Study.Name == "cfDNA LoB"]) %>%
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == cnv2$mut_key2 & cnv2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(cnv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.donor = rep(NA, length(vars))) %>%
  mutate(gDNA.call = vars %in% cnv2$mut_key2[cnv2$Study.Name == "gDNA LoB"]) %>%
  mutate(gDNA.maf = sapply(vars, function(v) {ii = which(v == cnv2$mut_key2 & cnv2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, mean(cnv2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(gDNA.donor = rep(NA, length(vars))) 
df.cnv.variant.pool$isoform = df.cnv.variant.pool$var.targeted.key = rep(NA, length(vars))
df.cnv.variant.pool$Ambry.target = df.cnv.variant.pool$gene %in% targeted.regions.ref$Gene
df.cnv.variant.pool$Ambry.maf = rep(NA, nrow(df.cnv.variant.pool))
df.cnv.variant.pool$MDL.call = df.cnv.variant.pool$mut_key %in% cnv.MDL$mut_key
df.cnv.variant.pool$MDL.maf = sapply(df.cnv.pool.all$mut_key, function(mut) cnv.MDL$percentage[mut == cnv.MDL$mut_key])
df.cnv.variant.pool$MDL.maf = cnv.MDL$copy_number[match(df.cnv.variant.pool$mut_key, cnv.MDL$mut_key)]
df.cnv.variant.pool$MDL.donor = sapply(df.cnv.variant.pool$mut_key, function(mut) 
{ii = which(cnv.MDL$mut_key == mut);
out = NA;
if(length(ii) > 0) out = type.MDL$Donor[type.MDL$run_sample_id == snv.MDL$run_sample_id[ii]];
out})

df.cnv.variant.pool

##  APPENDIX 3 TABLE
df.cnv.variants = df.cnv.variant.pool %>% dplyr::select(variant_type, gene, chrom, position, mut_nt, mut_aa) %>%
  mutate(ClinicalllySig = "No", nObs = df.cnv.variant.pool$nObs.pool) %>%  ## CHECK MANUALLY THE CASE 
  mutate(nTotal = sapply(df.cnv.variant.pool$Pool, function(p) sum(type$Pool == p & type$Study.Name == "pools LoB")))
df.cnv.variants = df.cnv.variants %>% mutate(FractionPos = rep(NA, nrow(df.cnv.variants))) %>%
  mutate(MAF = df.cnv.variant.pool$pool.maf)
df.cnv.variants 

##  APPENDIX 2 TABLE
zz = cnv.pool.all %>% mutate(isoform = rep(NA, nrow(cnv.pool.all)), exon = rep(NA, nrow(cnv.pool.all))) %>%
  mutate(Ambry.target.key = paste(gene, isoform, exon, sep=".")) %>%
  mutate(Ambry.target = Ambry.target.key %in%  ref.targeted.keys) %>%
  mutate(position = rep(NA, nrow(cnv.pool.all)), 
         mut_nt = rep(NA, nrow(cnv.pool.all)),
         mut_aa = rep(NA, nrow(cnv.pool.all)), 
         percentage = rep(NA, nrow(cnv.pool.all)), 
         fam_cnt = rep(NA, nrow(cnv.pool.all)), 
         mut_key1 = rep(NA, nrow(cnv.pool.all)))
df.cnv.sample =  zz %>% dplyr::select( run_sample_id,gene,chrom,position,mut_nt,mut_aa) %>% #, percentage, fam_cnt, mut_key1, Ambry.target) %>%
  mutate(ClinicallySig = "No") %>%  ##  MANUALLY CHECKED THIS ENTRY
  mutate(MAF = zz$percentage) %>%
  mutate(Confirmed = rep(FALSE, nrow(zz))) %>%
  #mutate(Confirmed = sapply(mut_key1, function(x) ifelse(any(x == indel.ref.df$mut_key2), "Yes", "No"))) %>%
  #mutate(Total_count = mol_cnt) %>%
  mutate(mol_cnt = rep(NA, nrow(zz))) %>%
  mutate(mut_cnt = rep(NA, nrow(zz)))
head(df.cnv.sample)

##  FUSIONs
fusion2 = fusion %>% left_join(type, by=c("run_sample_id", "runid"))
fusion2$mut_key1 = paste(fusion2$Pool,  fusion2$gene_a, fusion2$pos_a, sep = ".")

fusion.pool = fusion2 %>% filter(run_sample_id %in% type$run_sample_id[type$Study.Name == "pools LoB"] &
                                   call == 1)
cat(paste("FUSIONs =", nrow(fusion.pool), "...\n"))

vars = unique(fusion.pool$mut_key1)

df.fusion.variant.pool =  data.frame(chrom = fusion.pool$chrom_a[match(vars, fusion.pool$mut_key1)],
                                     position = fusion.pool$pos_a[match(vars, fusion.pool$mut_key1)],
                                     gene = fusion.pool$gene_a[match(vars, fusion.pool$mut_key1)],
                                     stringsAsFactors=F) %>%
  mutate(transcript_id =rep(NA, length(chrom))) %>%
  mutate(exon =rep(NA, length(chrom))) %>%
  mutate(mut_nt =rep(NA, length(chrom))) %>%
  mutate(mut_aa =rep(NA, length(chrom))) %>%
  mutate(mut_key =rep(NA, length(chrom))) %>%
  mutate(mut_key1 =rep(NA, length(chrom))) %>%
  mutate(variant_type = rep("FUSION", length(chrom))) %>% 
  mutate(Pool = fusion.pool$Pool[match(vars, fusion.pool$mut_key1)]) %>%
  mutate(run_sample_ids = sapply(vars, function(v) paste(fusion.pool$run_sample_id[v == fusion.pool$mut_key1], collapse=","))) %>%
  mutate(nObs.pool = sapply(vars, function(v) sum(v == fusion2$mut_key1[fusion2$Study.Name == "pools LoB"], na.rm=T))) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == fusion2$mut_key1 & fusion2$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(fusion2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.call = vars %in% fusion2$mut_key1[fusion2$Study.Name == "cfDNA LoB"]) %>%
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == fusion2$mut_key1 & fusion2$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(fusion2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(cfDNA.donor =rep(NA, length(chrom))) %>%
  mutate(gDNA.call = vars %in% fusion2$mut_key1[fusion2$Study.Name == "gDNA LoB"]) %>%
  mutate(gDNA.maf = sapply(vars, function(v) {ii = which(v == fusion2$mut_key1 & fusion2$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, mean(fusion2$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(gDNA.donor =rep(NA, length(chrom))) %>%
  mutate(isoform =rep(NA, length(chrom))) %>%
  mutate(var.targeted.key =rep(NA, length(chrom))) %>%
  mutate(Ambry.target =rep(NA, length(chrom))) %>%
  mutate(Ambry.maf =rep(NA, length(chrom)))
df.fusion.variant.pool$MDL.call = df.fusion.variant.pool$mut_key %in% fusion.MDL$mut_key
df.fusion.variant.pool$MDL.maf = fusion.MDL$percentage[match(df.fusion.variant.pool$mut_key, fusion.MDL$mut_key)]/100
df.fusion.variant.pool$MDL.donor = sapply(df.fusion.variant.pool$mut_key, function(mut) 
{ii = which(fusion.MDL$mut_key == mut);
out = NA;
if(length(ii) > 0) out = type.MDL$Donor[type.MDL$run_sample_id == fusion.MDL$run_sample_id[ii]];
out})

df.fusion.variant.pool
write.table(df.fusion.variant.pool, file=file.path(resdir, "Table2.All.FUSION.tsv"), 
            sep="\t", col.names=T, row.names=F, quote=F)

cols2keep = c("chrom", "position", "gene", "mut_aa", "variant_type", "Pool", "run_sample_ids", "nObs.pool", "pool.maf", 
              "cfDNA.call", "cfDNA.maf", "gDNA.call", "gDNA.maf", "Ambry.target", "Ambry.maf", "MDL.call", "MDL.maf", "MDL.donor")
df.variant.pool = rbind(df.snv.variant.pool, df.indel.variant.pool, df.cnv.variant.pool, df.fusion.variant.pool)
df.variant.pool

write.table(df.variant.pool, file=file.path(resdir, "vars.pool.table.tsv"), sep="\t", col.names=T, row.names=F, quote=F)


df.sample = rbind(df.snv.sample, df.indel.sample)
write.table(df.sample, file=file.path(resdir, "APPENDIX2.sample.table.tsv"), sep="\t", col.names=T, row.names=F, quote=F)

df.FP = rbind(df.snv.variants, df.indel.variants)
write.table(df.FP, file=file.path(resdir, "APPENDIX3.FP.table.tsv"), sep="\t", col.names=T, row.names=F, quote=F)

####################################
##  List of neat samples variants
####################################

##  SNVs
snv2 = snv %>% left_join(type, by=c("run_sample_id", "runid"))

snv.neat = snv2 %>% filter(run_sample_id %in% type$run_sample_id[type$Study.Name %in% c("cfDNA LoB", "gDNA LoB")] &
                             call == 1) 
vars = unique(snv.neat$mut_key1)
df.snv.neat = data.frame(chrom = snv.neat$chrom[match(vars, snv.neat$mut_key1)],
                         position = snv.neat$position[match(vars, snv.neat$mut_key1)],
                         gene = snv.neat$gene[match(vars, snv.neat$mut_key1)],
                         mut_nt = snv.neat$mut_nt[match(vars, snv.neat$mut_key1)],
                         mut_aa = snv.neat$mut_aa[match(vars, snv.neat$mut_key1)],
                         variant_type = "SNV",
                         cfDNA.call = vars %in% snv.neat$mut_key1[snv.neat$Study.Name == "cfDNA LoB"],
                         cfDNA.maf = sapply(vars, function(v) {ii = which(v == snv.neat$mut_key1 & snv.neat$Study.Name == "cfDNA LoB"); 
                         ifelse(length(ii) > 0, mean(snv.neat$percentage[ii]/100, na.rm=T), NA)}),
                         gDNA.call = vars %in% snv.neat$mut_key1[snv.neat$Study.Name == "gDNA LoB"],
                         gDNA.maf = sapply(vars, function(v) {ii = which(v == snv.neat$mut_key1 & snv.neat$Study.Name == "gDNA LoB"); 
                         ifelse(length(ii) > 0, mean(snv.neat$percentage[ii]/100, na.rm=T), NA)}),
                         pool.call = vars %in% snv$mut_key1[snv$run_sample_id %in% type$run_sample_id[type$Study.Name == "pools LoB"]],
                         pool.maf = sapply(vars, function(v) {ii = which(v == snv2$mut_key1 & snv2$Study.Name == "pools LoB"); 
                         ifelse(length(ii) > 0, mean(snv2$percentage[ii]/100, na.rm=T), NA)}),
                         stringsAsFactors=F)

##  INDELs
indel2 = indel %>% left_join(type, by=c("run_sample_id", "runid"))

indel.neat = indel2[indel2$run_sample_id %in% type$run_sample_id[type$Study.Name %in% c("cfDNA LoB", "gDNA LoB")] &
                      indel2$call == 1,]
vars = unique(indel.neat$mut_key1)
df.indel.neat = data.frame(chrom = indel.neat$chrom[match(vars, indel.neat$mut_key1)],
                           position = indel.neat$position[match(vars, indel.neat$mut_key1)],
                           gene = indel.neat$gene[match(vars, indel.neat$mut_key1)],
                           mut_nt = indel.neat$mut_nt[match(vars, indel.neat$mut_key1)],
                           mut_aa = indel.neat$mut_aa[match(vars, indel.neat$mut_key1)],
                           variant_type = "INDEL",
                           cfDNA.call = vars %in% indel.neat$mut_key[indel.neat$Study.Name == "cfDNA LoB"],
                           cfDNA.maf = sapply(vars, function(v) {ii = which(v == indel.neat$mut_key1 & indel.neat$Study.Name == "cfDNA LoB"); 
                           ifelse(length(ii) > 0, mean(indel.neat$percentage[ii]/100, na.rm=T), NA)}),
                           gDNA.call = vars %in% indel.neat$mut_key[indel.neat$Study.Name == "gDNA LoB"],
                           gDNA.maf = sapply(vars, function(v) {ii = which(v == indel.neat$mut_key1 & indel.neat$Study.Name == "gDNA LoB"); 
                           ifelse(length(ii) > 0, mean(indel.neat$percentage[ii]/100, na.rm=T), NA)}),
                           pool.call = vars %in% indel$mut_key[indel$run_sample_id %in% type$run_sample_id[type$Study.Name == "pools LoB"]],
                           pool.maf = sapply(vars, function(v) {ii = which(v == indel2$mut_key1 & indel2$Study.Name == "pools LoB"); 
                           ifelse(length(ii) > 0, mean(indel2$percentage[ii]/100, na.rm=T), NA)}),
                           stringsAsFactors=F)

##  CNVs
#summary = make_gene_arm_summary(probes, params)
#genes2keep = summary$gene.table$gene[summary$gene.table$reportable.focal == 1 & summary$gene.table$reportable == 1]
genes2keep = c("ERBB2", "MET")
cnv2 = cnv %>% left_join(type, by=c("run_sample_id", "runid"))
cnv2$mut_key2 = paste(cnv$Pool, cnv$gene, "CNV", sep=".")

cnv.neat = cnv2[cnv2$run_sample_id %in% type$run_sample_id[type$Study.Name %in% c("cfDNA LoB")] &
                  cnv2$call == 1 & cnv2$focal_call == 1 & cnv2$gene %in% genes2keep,]
cat(paste("CNVs =", nrow(cnv.neat), "...\n"))

vars = unique(cnv.neat$mut_key2)
df.cnv.neat = data.frame(chrom = cnv.neat$chrom[match(vars, cnv.neat$mut_key2)],
                         position = cnv.neat$position[match(vars, cnv.neat$mut_key2)],
                         gene = cnv.neat$gene[match(vars, cnv.neat$mut_key2)],
                         stringsAsFactors=F) %>%  
  mutate(mut_nt = rep(NA, length(chrom))) %>%
  mutate(mut_aa = rep(NA, length(chrom))) %>%
  mutate(mut_nt = rep(NA, length(chrom))) %>%
  mutate(cfDNA.call = vars %in% cnv.neat$mut_key2[cnv.neat$Study.Name == "cfDNA LoB"]) %>%
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == cnv.neat$mut_key2 & cnv.neat$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(cnv.neat$copy_number[ii], na.rm=T), NA)})) %>%
  mutate(gDNA.call = vars %in% cnv.neat$mut_key[cnv.neat$Study.Name == "gDNA LoB"],
         cfDNA.maf = sapply(vars, function(v) {ii = which(v == cnv.neat$mut_key2 & cnv.neat$Study.Name == "gDNA LoB"); 
         ifelse(length(ii) > 0, mean(cnv.neat$copy_number[ii], na.rm=T), NA)})) %>%
  mutate(pool.call = vars %in% cnv$mut_key[cnv$Study.Name == "pools LoB"]) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == cnv$mut_key2 & cnv$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(cnv2$copy_number[ii], na.rm=T), NA)}))

# FUSIONs
fusion2 = fusion %>% left_join(type, by=c("run_sample_id", "runid"))
fusion2$mut_key2 = paste(fusion2$Pool, fusion2$gene_a, fusion2$pos_a, fusion2$mut_nt, sep="_")

fusion.neat = fusion2 %>% filter(run_sample_id %in% type$run_sample_id[type$Study.Name %in% c("cfDNA LoB", "gDNA LoB")] &
                                   call == 1) 

cat(paste("FUSIONs =", nrow(fusion.neat), "...\n"))

vars = unique(fusion.neat$mut_key)
df.fusion.neat = data.frame(chrom = fusion.neat$chrom_a[match(vars, fusion.neat$mut_key2)],
                            position = fusion.neat$pos_a[match(vars, fusion.neat$mut_key2)],
                            gene = fusion.neat$gene_a[match(vars, fusion.neat$mut_key2)],
                            stringsAsFactors=F) %>%
  mutate(mut_nt = rep(NA, length(chrom))) %>%
  mutate(mut_aa = rep(NA, length(chrom))) %>%
  mutate(cfDNA.call = vars %in% fusion.neat$mut_key2[fusion.neat$Study.Name == "cfDNA LoB"]) %>%
  mutate(cfDNA.maf = sapply(vars, function(v) {ii = which(v == fusion.neat$mut_key2 & fusion.neat$Study.Name == "cfDNA LoB"); 
  ifelse(length(ii) > 0, mean(fusion.neat$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(gDNA.call = vars %in% fusion.neat$mut_key2[fusion.neat$Study.Name == "gDNA LoB"]) %>%
  mutate(gDNA.maf = sapply(vars, function(v) {ii = which(v == fusion.neat$mut_key2 & fusion.neat$Study.Name == "gDNA LoB"); 
  ifelse(length(ii) > 0, mean(fusion.neat$percentage[ii]/100, na.rm=T), NA)})) %>%
  mutate(pool.call = vars %in% fusion$mut_key[fusion$Study.Name == "pools LoB"]) %>%
  mutate(pool.maf = sapply(vars, function(v) {ii = which(v == fusion$mut_key2 & fusion$Study.Name == "pools LoB"); 
  ifelse(length(ii) > 0, mean(fusion$percentage[ii]/100, na.rm=T), NA)}))


df.neat = rbind(df.snv.neat, df.indel.neat, df.cnv.neat, df.fusion.neat) %>% arrange(as.numeric(chrom), position)
df.neat

write.table(df.neat, file=file.path(resdir, "Table3.neat.table.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
