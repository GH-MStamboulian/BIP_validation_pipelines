##  Two methods:
##  LBP70: targeted MAFs adjusted by LBP70 observations, as in protocol
##  G360CDx: Level5 mean of observed MAFs, subsequent levels are adjustments of Level5 by dilution level

###############################
##  Table 6&7: TRUE MAF values
##  LBP70 and G360CDx methods
###############################

  #load()

  
  ##################################################
  ##  TRUE MAF LBP70 method
  ##  adjusted MAF by MDL observations (protocol)
  ##################################################
  
  ###########
  ##  SNVs
  ###########
  
  expected.maf = targeted.maf %>% filter(variant_type == "SNV")
  
  snvs.target = snv %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    filter(!is.na(mut_aa)) %>%
    mutate(mutation = paste(gene, mut_aa, sep = '.'), lbl = paste(mut_key, Type, ng_seq)) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  snvs.target = snvs.target %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq", "mut_aa")) %>%
    filter(Dilution == "Level5")
  #table(snvs.target$mut_key, snvs.target$ng_seq, snvs.target$Type)
  
 
  
  ###########
  ##  INDELs
  ###########
  
  expected.maf = targeted.maf %>% filter(variant_type == "INDEL")
  
  indels.target = indel %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    filter(!is.na(mut_aa)) %>%
    dplyr::mutate(mutation = paste(gene, mut_aa, sep = '.'), lbl = paste(mut_key, Type, ng_seq)) %>%
    dplyr::filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  indels.target = indels.target %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq", "mut_aa")) %>%
    filter(Dilution == "Level5")
  #table(indels.target$mut_key, indels.target$ng_seq, indels.target$Type)
  

  
  ##############
  ##  CNAs
  ##############
  
  expected.maf = targeted.maf %>% filter(variant_type == "CNV")
  
  cnvs.target = cnv %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    dplyr::mutate(mutation = gene, mut_key = gene, lbl = paste(gene, Type, ng_seq)) %>%
    dplyr::filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  cnvs.target = cnvs.target %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq")) %>%
    filter(Dilution == "Level5")
  cnvs.target$mut_aa = cnvs.target$gene
  #table(cnvs.target$mutation, cnvs.target$ng_seq, cnvs.target$Type)
  

  
  ##############
  ##  FUSIONs
  ##############
  
  expected.maf = targeted.maf %>% filter(variant_type == "FUSION")
  
  fusions.target = fusion %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    filter(!is.na(gene_a)) %>%
    filter(call == "1") %>%
    filter(downstream_gene == "A") %>%
    dplyr::mutate(mutation = paste("CL", gene_a, gene_b, sep = '_'), lbl = paste(Type, ng_seq, gene_a, gene_b)) %>%
    dplyr::mutate(gene = gene_a, mut_aa = paste("CL", gene_a, gene_b, sep = '_'), mut_key = paste("CL", gene_a, gene_b, sep = '_')) %>%
    dplyr::filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  fusions.target = fusions.target %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq", "mut_aa")) %>%
    filter(Dilution == "Level5")
  #table(fusions.target$mutation, fusions.target$ng_seq, fusions.target$Type)
  

  

  
 
  ########################################################
  ##  TRUE MAF = G360 CDx method
  ##  MAF = highest level MAF adjusted by dilution level
  ## targeted.maf = observed MAF
  ########################################################
  
  targeted.maf = targeted.maf.init
  expected.maf = targeted.maf %>% filter(variant_type == "SNV")
  snvs.target = snv %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    filter(!is.na(mut_aa)) %>%
    mutate(mutation = paste(gene, mut_aa, sep = '.'), lbl = paste(mut_key, Type, ng_seq)) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  snvs.target = snvs.target %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq", "mut_aa")) 
  snvs.df = snvs.target %>%
    dplyr::group_by(gene, mut_aa, mutation, ng_seq, Type, Dilution) %>%
    dplyr::summarize(targeted_maf = mean(target_maf), observed_maf = median(percentage) / 100, nObs = length(target_maf)) %>%
    mutate(Gene = gene, Pool = Type)
  snvs.df = expected.maf %>% 
    left_join(snvs.df, by=c("Gene", "mut_aa",  "mutation", "ng_seq", "Dilution")) %>% 
    dplyr::select(Gene, mut_aa, mutation, variant_type, ng_seq, Dilution, Type, targeted_maf, observed_maf, nObs)
  
  expected.maf = targeted.maf %>% filter(variant_type == "INDEL")
  indels.target = indel %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    filter(!is.na(mut_aa)) %>%
    dplyr::mutate(mutation = paste(gene, mut_aa, sep = '.'), lbl = paste(mut_key, Type, ng_seq)) %>%
    dplyr::filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  indels.target = indels.target %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq", "mut_aa"))
  indels.df = indels.target %>%
    dplyr::group_by(Gene, mut_aa, mutation, ng_seq, Type, Dilution) %>%
    dplyr::summarize(targeted_maf = mean(target_maf), observed_maf = median(percentage) / 100, nObs = length(target_maf))
  indels.df = expected.maf %>% 
    left_join(indels.df, by=c("Gene", "mut_aa",  "mutation", "ng_seq", "Dilution")) %>% 
    dplyr::select(Gene, mut_aa, mutation, variant_type, ng_seq, Dilution, Type, targeted_maf, observed_maf, nObs)
  
  expected.maf = targeted.maf %>% filter(variant_type == "CNV")
  cnvs.target = cnv %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    dplyr::mutate(mutation = gene, mut_key = gene, lbl = paste(gene, Type, ng_seq)) %>%
    dplyr::filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  cnvs.target = cnvs.target %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq"))
  cnvs.df = cnvs.target %>%
    dplyr::group_by(Gene, mut_aa, mutation, ng_seq, Type, Dilution) %>%
    dplyr::summarize(targeted_maf = mean(target_maf), observed_maf = median(copy_number), nObs = length(target_maf)) 
  cnvs.df = expected.maf %>% 
    left_join(cnvs.df, by=c("Gene", "mut_aa",  "mutation", "ng_seq", "Dilution")) %>% 
    dplyr::select(Gene, mut_aa, mutation, variant_type, ng_seq, Dilution, Type, targeted_maf, observed_maf, nObs)
  
  expected.maf = targeted.maf %>% filter(variant_type == "FUSION")
  fusions.target = fusion %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    filter(!is.na(gene_a)) %>%
    filter(call == "1") %>%
    filter(downstream_gene == "A") %>%
    dplyr::mutate(mutation = paste("CL", gene_a, gene_b, sep = '_'), lbl = paste(Type, ng_seq, gene_a, gene_b)) %>%
    dplyr::mutate(gene = gene_a, mut_aa = paste("CL", gene_a, gene_b, sep = '_'), mut_key = paste("CL", gene_a, gene_b, sep = '_')) %>%
    dplyr::filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  fusions.target = fusions.target %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq", "mut_aa"))
  fusions.df = fusions.target %>%
    dplyr::group_by(Gene, mut_aa, mutation, ng_seq, Type, Dilution) %>%
    dplyr::summarize(targeted_maf = mean(target_maf), observed_maf = median(percentage) / 100, nObs = length(target_maf)) 
  fusions.df = expected.maf %>% 
    left_join(fusions.df, by=c("Gene", "mut_aa",  "mutation", "ng_seq", "Dilution")) %>% 
    dplyr::select(Gene, mut_aa, mutation, variant_type, ng_seq, Dilution, Type, targeted_maf, observed_maf, nObs)
  
  df.v3 = rbind(snvs.df, indels.df, cnvs.df, fusions.df)
  write.table(df.v3, file=file.path(resdir, "Table7.MAF.G360CDx.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  if(0) {
    targeted.maf.new.v3 = targeted.maf.init
    targeted.maf.new.v3$target_maf = sapply(1:nrow(targeted.maf.new.v3), function(i) 
    {j = which(df.v3$mutation == targeted.maf.new.v3$mutation[i] & 
                 df.v3$ng_seq == targeted.maf.new.v3$ng_seq[i] & 
                 df.v3$Dilution == targeted.maf.new.v3$Dilution[i]);
    out = df.v3$observed_maf[j];
    if(is.na(out)) out = df.v3$targeted_maf[j]; out})
    targeted.maf.new.v3$replacedTarget = TRUE
  }
  
  # targeted MAF = G360 obs MAF for Level5 
  # and diluted for remaining levels
  
  targeted.maf.new.v3 = targeted.maf.init
  targeted.maf.new.v3$target_maf = sapply(1:nrow(targeted.maf.new.v3), function(i) 
  {k = which(targeted.maf.init$mutation == targeted.maf.init$mutation[i] &
               targeted.maf.init$ng_seq == targeted.maf.init$ng_seq[i] &
               targeted.maf.init$Dilution == "Level5");
  j = which(df.v3$mutation == targeted.maf.init$mutation[i] & 
              df.v3$ng_seq == targeted.maf.init$ng_seq[i] &
              df.v3$Dilution == "Level5");
  x = df.v3$observed_maf[j] * (targeted.maf.init$target_maf[i] / targeted.maf.init$target_maf[k]);
  if(targeted.maf.new.v3$variant_type[i] == "CNV") 
    x = 2 + (df.v3$observed_maf[j]-2) * (targeted.maf.init$target_maf[i]-2) / (targeted.maf.init$target_maf[k]-2);
  #if(!(targeted.maf.init$target_maf[k] == df.v3$targeted_maf[j])) x = NA;
  x})
  targeted.maf.new.v3$replacedTarget =  TRUE;
  
  write.table(targeted.maf.new.v3, file=file.path(resdir, "Table7.targeted.maf.G360CDx.tsv"), 
              sep="\t", col.names=T, row.names=F, quote=F)
  save(targeted.maf.new.v3, file=file.path(resdir, "Table7.targeted.maf.G360CDx.RData"))
  
  ####\#########################################
  ##  Make TABLE 7 G360CDx method 
  ##  with TRUE MAF / variant / Dilution level
  ##############################################
  
  dilutions = paste0("Level", 1:5)
  ng_seqs = c("5 ng", "30 ng")
  vars = unique(targeted.maf.new.v3$mutation)
  
  targeted.maf.new.v3.df = NULL
  for(ng_seq in ng_seqs){
    for(var in vars) {
      maf = sapply(dilutions, function(d) targeted.maf.new.v3$target_maf[targeted.maf.new.v3$mutation == var &
                                                                           targeted.maf.new.v3$ng_seq == ng_seq &
                                                                           targeted.maf.new.v3$Dilution == d])
      gene = targeted.maf.new.v3$Gene[targeted.maf.new.v3$mutation == var][1]
      mut_aa = targeted.maf.new.v3$mut_aa[targeted.maf.new.v3$mutation == var][1]
      variant_type = targeted.maf.new.v3$variant_type[targeted.maf.new.v3$mutation == var][1]
      Pool = targeted.maf.new.v3$Pool[targeted.maf.new.v3$mutation == var][1]
      replacedTarget = targeted.maf.new.v3$replacedTarget[targeted.maf.new.v3$mutation == var][1]
      targeted.maf.new.v3.df = rbind(targeted.maf.new.v3.df, 
                                     c(gene = gene, 
                                       mut_aa = mut_aa, 
                                       variant_type = variant_type, 
                                       Pool = Pool, 
                                       mutation = var,
                                       ng_seq = ng_seq,
                                       replacedTarget = replacedTarget,
                                       maf))
    }
  }
  targeted.maf.new.v3.df = as.data.frame(targeted.maf.new.v3.df, stringsAsFactors=F)
  class(targeted.maf.new.v3.df$Level1) = class(targeted.maf.new.v3.df$Level2) = 
    class(targeted.maf.new.v3.df$Level3) = class(targeted.maf.new.v3.df$Level4) = class(targeted.maf.new.v3.df$Level5) = "numeric"
  
  targeted.maf.new.v3.df
  write.table(targeted.maf.new.v3.df, file=file.path(resdir, "Table7.targeted.maf.level.G360CDx.tsv"), 
              sep="\t", col.names=T, row.names=F, quote=F)
  