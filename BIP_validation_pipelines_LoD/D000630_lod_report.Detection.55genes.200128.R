###############################
##  Table 4: Detection rates
##  Independent of TRUE MAF method
###############################

  ################################
  ##  Initialize targeted.maf with 
  ##  G360CDx method values
  ##  Used only for variants info
  ################################

  targeted.maf = targeted.maf.new.v3

  ##############################################################
  ##  Section 10.2.2.3	
  ##  Guardant360 CDx result data set derivation for 55 genes
  ##  i.e. BIP v3.5.3 method of detection: rm_reportable = 1
  ##       but CNV call %in% 2:3
  ##   Skips filtering by rm_reportable for all versions of BIP
  ###############################################################
  
  #########
  ##  SNVs
  #########
  
  
  expected.maf = targeted.maf %>% filter(variant_type == "SNV")
  snv_table.filter <- snv %>% left_join(type, by = c("run_sample_id", "runid")) %>%
  #if ("rm_reportable" %in% colnames(snv)) {
  #  snv_table.filter <- snv_table.filter %>% filter(rm_reportable == 1)
  #} else {
  #  snv_table.filter <- snv_table.filter %>% filter(call == 1)
  #}
  #snv_table.filter <- snv_table.filter %>%
    filter(!is.na(Dilution)) %>%
    filter(!is.na(mut_aa)) %>%
    mutate(mutation = paste(gene, mut_aa, sep = '.'), lbl = paste(mut_key, Type, ng_seq)) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  snv_table.filter = snv_table.filter %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq"))
  ta = table(snv_table.filter$mut_key, snv_table.filter$Dilution, snv_table.filter$ng_seq)
  
  mut_keys = rownames(ta[,,1])
  ng_seq = colnames(ta[,1,])
  dilution = colnames(ta[,,1])
  
  mut_key = rep(mut_keys, length(ng_seq))
  Input = rep(ng_seq, each=length(mut_keys))
  Pool = rep(targeted.maf$Pool[match(mut_keys, targeted.maf$mutation)], length(ng_seq))
  nReplicates = sapply(1:length(mut_key), function(i) sum(type$ng_seq == Input[i] & type$Sample_Description5 == Pool[i] & type$Dilution == "Level5"))
  
  detection.snv.filter = data.frame(mut_key = mut_key,
                                    Input = Input,
                                    Pool = Pool,
                                    variant_type = rep("SNV", length(mut_key)),
                                    nReplicates = nReplicates,
                                    rbind(ta[,,1], ta[,,2]))
  
  cat("Number of variants detected:\n")
  detection.snv.filter
  
  
  #############
  ##  INDELs
  #############
  
  expected.maf = targeted.maf %>% filter(variant_type == "INDEL")
  indel_table.filter <- indel %>% left_join(type, by = c("run_sample_id")) %>%
  #  if ("rm_reportable" %in% colnames(snv)) {
  #    indel_table.filter <- indel_table.filter %>% filter(rm_reportable == 1)
  #  } else {
  #    indel_table.filter <- indel_table.filter %>% filter(call == 1)
  #  }
  #indel_table.filter <- indel_table.filter %>%
    filter(!is.na(Dilution)) %>%
    filter(!is.na(mut_aa)) %>%
    mutate(mutation = paste(gene, mut_aa, sep = '.'), lbl = paste(Type, ng_seq, mut_key)) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  indel_table.filter = indel_table.filter %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq"))  
  
  ta = table(indel_table.filter$mut_key, indel_table.filter$Dilution, indel_table.filter$ng_seq)
  mut_keys = rownames(ta[,,1])
  ng_seq = colnames(ta[,1,])
  dilution = colnames(ta[,,1])
  
  mut_key = rep(mut_keys, length(ng_seq))
  Input = rep(ng_seq, each=length(mut_keys))
  Pool = rep(targeted.maf$Pool[match(mut_keys, targeted.maf$mutation)], length(ng_seq))
  nReplicates = sapply(1:length(mut_key), function(i) sum(type$ng_seq == Input[i] & type$Sample_Description5 == Pool[i] & type$Dilution == "Level5"))
  
  detection.indel.filter = data.frame(mut_key = mut_key,
                                      Input = Input,
                                      Pool = Pool,
                                      variant_type = rep("INDEL", length(mut_key)),
                                      nReplicates = nReplicates,
                                      rbind(ta[,,1], ta[,,2]))
  
  cat("Number of variants detected:\n")
  detection.indel.filter
  
  ##########
  ##  CNAs
  ##########
  
  expected.maf = targeted.maf %>% filter(variant_type == "CNV")
  cna_table.filter <- cnv %>% left_join(type, by = c("run_sample_id")) %>%
    #filter(call %in% 2:3) %>%
    filter(!is.na(Dilution)) %>%
    filter(!is.na(gene)) %>%
    mutate(mutation = gene) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  cna_table.filter = cna_table.filter %>% left_join(expected.maf, by = c("Dilution", "ng_seq", "mutation"))
  cna_table.filter$lbl = paste(cna_table.filter$gene, cna_table.filter$ng_seq)
  
  zz = cna_table.filter %>% 
    group_by(gene, ng_seq, Dilution) %>%
    summarise(obs_maf = median(copy_number))
  cna_table.filter = cna_table.filter %>% 
    left_join(zz, by = c("gene", "ng_seq", "Dilution"))
  cna_table.filter$Dilution = as.factor(cna_table.filter$Dilution)
  cna_table.filter$call = as.numeric(cna_table.filter$call %in% 2:3)
  iDetect = which(cna_table.filter$call == 1)
  
  ta = table(cna_table.filter$mutation[iDetect], cna_table.filter$Dilution[iDetect], cna_table.filter$ng_seq[iDetect])
  mut_keys = rownames(ta[,,1])
  ng_seq = colnames(ta[,1,])
  dilution = colnames(ta[,,1])
  
  mut_key = rep(mut_keys, length(ng_seq))
  Input = rep(ng_seq, each=length(mut_keys))
  Pool = rep(targeted.maf$Pool[match(mut_keys, targeted.maf$mutation)], length(ng_seq))
  nReplicates = sapply(1:length(mut_key), function(i) sum(type$ng_seq == Input[i] & type$Sample_Description5 == Pool[i] & type$Dilution == "Level5"))
  
  detection.cnv.filter = data.frame(mut_key = mut_key,
                                    Input = Input,
                                    Pool = Pool,
                                    variant_type = rep("CNV", length(mut_key)),
                                    nReplicates = nReplicates,
                                    rbind(ta[,,1], ta[,,2]))
  
  cat("Number of variants detected:\n")
  detection.cnv.filter
  
  #############
  ##  FUSIONs
  #############
  
  expected.maf = targeted.maf %>% filter(variant_type == "FUSION")
  fusion_table.filter <- fusion %>% left_join(type, by = c("run_sample_id")) %>%
    filter(!is.na(Dilution)) %>%
    filter(!is.na(gene_a)) %>%
    filter(call == "1") %>%
    filter(downstream_gene == "A") %>%
    mutate(mutation = paste("CL", gene_a, gene_b, sep = '_'), lbl = paste(Type, ng_seq, gene_a, gene_b)) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  fusion_table.filter = fusion_table.filter %>% left_join(expected.maf, by = c("Dilution", "ng_seq", "mutation"))
  
  ta = table(fusion_table.filter$mutation, fusion_table.filter$Dilution, fusion_table.filter$ng_seq)
  mut_keys = rownames(ta[,,1])
  ng_seq = colnames(ta[,1,])
  dilution = colnames(ta[,,1])
  
  mut_key = rep(mut_keys, length(ng_seq))
  Input = rep(ng_seq, each=length(mut_keys))
  Pool = rep(targeted.maf$Pool[match(mut_keys, targeted.maf$mutation)], length(ng_seq))
  nReplicates = sapply(1:length(mut_key), function(i) sum(type$ng_seq == Input[i] & type$Sample_Description5 == Pool[i] & type$Dilution == "Level5"))
  
  detection.fusion.filter = data.frame(mut_key = mut_key,
                                       Input = Input,
                                       Pool = Pool,
                                       variant_type = rep("FUSION", length(mut_key)),
                                       nReplicates = nReplicates,
                                       rbind(ta[,,1], ta[,,2]))
  
  cat("Number of variants detected:\n")
  detection.fusion.filter
  
  ################################
  ##  Aggregate detection counts 
  ##  individual variants
  ################################
  
  detection.filter = rbind(detection.snv.filter, detection.indel.filter, detection.cnv.filter, detection.fusion.filter)
  write.table(detection.filter, file=file.path(resdir, "Table4.detection.variants.G360CDx.55genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  ###############################
  ##  Generate derection rates
  ##  individual variants
  ###############################
  
  detection.rate.all.filter = detection.filter %>%  arrange(desc(Input))
  for(n in paste0("Level", 1:5)) detection.rate.all.filter[[n]] = round(detection.rate.all.filter[[n]]/detection.rate.all.filter$nReplicates, 4)
  detection.rate.all.filter
  write.table(detection.rate.all.filter, file=file.path(resdir, "Table4.detection.rates.variants.G360CDx.55genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  #######################################
  ##  Aggregate accross variant class
  ##  Required Outputs of detection
  #######################################
  
  CDx.mut_key = c("EGFR.L858R", "EGFR.T790M", "EGFR.p.Glu746_Ala750del")
  detection.snv.aggregated.filter = detection.snv.filter %>% group_by(Input) %>% summarise(nReplicates=sum(nReplicates), 
                                                                                           Level1 = sum(Level1), 
                                                                                           Level2 = sum(Level2), 
                                                                                           Level3 = sum(Level3),
                                                                                           Level4 = sum(Level4), 
                                                                                           Level5 = sum(Level5)) %>%
    mutate(mut_key="aggregated", Pool = NA, variant_type = "SNV")
  detection.snv.aggregated.filter = detection.snv.aggregated.filter[,colnames(detection.snv.filter)]
  detection.indel.aggregated.filter = detection.indel.filter %>% group_by(Input) %>% summarise(nReplicates=sum(nReplicates), 
                                                                                               Level1 = sum(Level1), 
                                                                                               Level2 = sum(Level2), 
                                                                                               Level3 = sum(Level3),
                                                                                               Level4 = sum(Level4), 
                                                                                               Level5 = sum(Level5)) %>%
    mutate(mut_key="aggregated", Pool = NA, variant_type = "INDEL")
  detection.indel.aggregated.filter = detection.indel.aggregated.filter[,colnames(detection.indel.filter)]
  
  detection.final.filter = rbind(detection.snv.filter %>% filter(mut_key %in% CDx.mut_key),
                                 detection.indel.filter %>% filter(mut_key %in% CDx.mut_key), 
                                 detection.snv.aggregated.filter, 
                                 detection.indel.aggregated.filter, 
                                 detection.cnv.filter,
                                 detection.fusion.filter) %>%  arrange(desc(Input))
  
  detection.final.filter
  write.table(detection.final.filter, file=file.path(resdir, "Table4.detection.G360CDx.55genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  detection.rate.final.filter = detection.final.filter
  for(n in paste0("Level", 1:5)) detection.rate.final.filter[[n]] = round(detection.rate.final.filter[[n]]/detection.rate.final.filter$nReplicates, 4)
  detection.rate.final.filter
  write.table(detection.rate.final.filter, file=file.path(resdir, "Table4.detection.rates.G360CDx.55genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  
