###############################
##  Table 4: Detection rates
##  Independent of TRUE MAF method
###############################

  #load()

  ################################
  ##  Initialize targeted.maf with 
  ##  G360CDx method values
  ##  Used only for variants info
  ################################
  
  targeted.maf = targeted.maf.new.v3
  
  ##############################################################
  ##  Section 10.2.2.2	
  ##  Guardant360 CDx result data set derivation for 74 genes
  ##  i.e. BIP v3.5.2 method of detection
  ###############################################################
  
  #########
  ##  SNVs
  #########
  
  expected.maf = targeted.maf %>% filter(variant_type == "SNV")
  snv_table <- snv %>% left_join(type, by = c("run_sample_id", "runid")) %>% 
    filter(!is.na(Dilution)) %>%
    filter(!is.na(mut_aa)) %>%
    mutate(mutation = paste(gene, mut_aa, sep = '.'), lbl = paste(mut_key, Type, ng_seq)) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  snv_table = snv_table %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq"))
  ta = table(snv_table$mut_key, snv_table$Dilution, snv_table$ng_seq)
  
  mut_keys = rownames(ta[,,1])
  ng_seq = colnames(ta[,1,])
  dilution = colnames(ta[,,1])
  
  mut_key = rep(mut_keys, length(ng_seq))
  Input = rep(ng_seq, each=length(mut_keys))
  Pool = rep(targeted.maf$Pool[match(mut_keys, targeted.maf$mutation)], length(ng_seq))
  nReplicates = sapply(1:length(mut_key), function(i) sum(type$ng_seq == Input[i] & type$Sample_Description5 == Pool[i] & type$Dilution == "Level5"))
  
  detection.snv = data.frame(mut_key = mut_key,
                             Input = Input,
                             Pool = Pool,
                             variant_type = rep("SNV", length(mut_key)),
                             nReplicates = nReplicates,
                             rbind(ta[,,1], ta[,,2]))
  
  cat("Number of variants detected:\n")
  detection.snv
  
  
  #############
  ##  INDELs
  #############
  
  expected.maf = targeted.maf %>% filter(variant_type == "INDEL")
  indel_table <- indel %>% left_join(type, by = c("run_sample_id")) %>%
    #group_by(Type, gene, mut_aa, Dilution) %>%
    filter(!is.na(Dilution)) %>%
    filter(!is.na(mut_aa)) %>%
    mutate(mutation = paste(gene, mut_aa, sep = '.'), lbl = paste(Type, ng_seq, mut_key)) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  indel_table = indel_table %>% left_join(expected.maf, by = c("mutation", "Dilution", "ng_seq"))  
  
  ta = table(indel_table$mut_key, indel_table$Dilution, indel_table$ng_seq)
  mut_keys = rownames(ta[,,1])
  ng_seq = colnames(ta[,1,])
  dilution = colnames(ta[,,1])
  
  mut_key = rep(mut_keys, length(ng_seq))
  Input = rep(ng_seq, each=length(mut_keys))
  Pool = rep(targeted.maf$Pool[match(mut_keys, targeted.maf$mutation)], length(ng_seq))
  nReplicates = sapply(1:length(mut_key), function(i) sum(type$ng_seq == Input[i] & type$Sample_Description5 == Pool[i] & type$Dilution == "Level5"))
  
  detection.indel = data.frame(mut_key = mut_key,
                               Input = Input,
                               Pool = Pool,
                               variant_type = rep("INDEL", length(mut_key)),
                               nReplicates = nReplicates,
                               rbind(ta[,,1], ta[,,2]))
  
  cat("Number of variants detected:\n")
  detection.indel
  
  ##########
  ##  CNAs
  ##########
  
  expected.maf = targeted.maf %>% filter(variant_type == "CNV")
  cna_table <- cnv %>% left_join(type, by = c("run_sample_id")) %>%
    #group_by(Type, gene, mut_aa, Dilution) %>%
    filter(!is.na(Dilution)) %>%
    filter(!is.na(gene)) %>%
    mutate(mutation = gene) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  cna_table = cna_table %>% left_join(expected.maf, by = c("Dilution", "ng_seq", "mutation"))
  cna_table$lbl = paste(cna_table$gene, cna_table$ng_seq)
  
  zz = cna_table %>% 
    group_by(gene, ng_seq, Dilution) %>%
    summarise(obs_maf = median(copy_number))
  cna_table = cna_table %>% 
    left_join(zz, by = c("gene", "ng_seq", "Dilution"))
  cna_table$Dilution = as.factor(cna_table$Dilution)
  cna_table$call = as.numeric(cna_table$call %in% 2:3)
  iDetect = which(cna_table$call == 1)
  
  ta = table(cna_table$mutation[iDetect], cna_table$Dilution[iDetect], cna_table$ng_seq[iDetect])
  mut_keys = rownames(ta[,,1])
  ng_seq = colnames(ta[,1,])
  dilution = colnames(ta[,,1])
  
  mut_key = rep(mut_keys, length(ng_seq))
  Input = rep(ng_seq, each=length(mut_keys))
  Pool = rep(targeted.maf$Pool[match(mut_keys, targeted.maf$mutation)], length(ng_seq))
  nReplicates = sapply(1:length(mut_key), function(i) sum(type$ng_seq == Input[i] & type$Sample_Description5 == Pool[i] & type$Dilution == "Level5"))
  
  detection.cnv = data.frame(mut_key = mut_key,
                             Input = Input,
                             Pool = Pool,
                             variant_type = rep("CNV", length(mut_key)),
                             nReplicates = nReplicates,
                             rbind(ta[,,1], ta[,,2]))
  
  cat("Number of variants detected:\n")
  detection.cnv
  
  #############
  ##  FUSIONs
  #############
  
  expected.maf = targeted.maf %>% filter(variant_type == "FUSION")
  fusion_table <- fusion %>% left_join(type, by = c("run_sample_id")) %>%
    #group_by(Type, gene, mut_aa, Dilution) %>%
    filter(!is.na(Dilution)) %>%
    filter(!is.na(gene_a)) %>%
    filter(call == "1") %>%
    filter(downstream_gene == "A") %>%
    mutate(mutation = paste("CL", gene_a, gene_b, sep = '_'), lbl = paste(Type, ng_seq, gene_a, gene_b)) %>%
    filter((Type == "poolA" & mutation %in% uvPoolA) | (Type == "poolB" & mutation %in% uvPoolB) ) 
  fusion_table = fusion_table %>% left_join(expected.maf, by = c("Dilution", "ng_seq", "mutation"))
  
  ta = table(fusion_table$mutation, fusion_table$Dilution, fusion_table$ng_seq)
  mut_keys = rownames(ta[,,1])
  ng_seq = colnames(ta[,1,])
  dilution = colnames(ta[,,1])
  
  mut_key = rep(mut_keys, length(ng_seq))
  Input = rep(ng_seq, each=length(mut_keys))
  Pool = rep(targeted.maf$Pool[match(mut_keys, targeted.maf$mutation)], length(ng_seq))
  nReplicates = sapply(1:length(mut_key), function(i) sum(type$ng_seq == Input[i] & type$Sample_Description5 == Pool[i] & type$Dilution == "Level5"))
  
  detection.fusion = data.frame(mut_key = mut_key,
                                Input = Input,
                                Pool = Pool,
                                variant_type = rep("FUSION", length(mut_key)),
                                nReplicates = nReplicates,
                                rbind(ta[,,1], ta[,,2]))
  
  cat("Number of variants detected:\n")
  detection.fusion
  
  ################################
  ##  Aggregate detection counts 
  ##  individual variants
  ################################
  
  detection = rbind(detection.snv, detection.indel, detection.cnv, detection.fusion)
  write.table(detection, file=file.path(resdir, "Table4.detection.variants.G360CDx.74genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  ###############################
  ##  Generate derection rates
  ##  individual variants
  ###############################
  
  detection.rate.all = detection %>%  arrange(desc(Input))
  for(n in paste0("Level", 1:5)) detection.rate.all[[n]] = round(detection.rate.all[[n]]/detection.rate.all$nReplicates, 4)
  detection.rate.all
  write.table(detection.rate.all, file=file.path(resdir, "Table4.detection.rates.variants.G360CDx.74genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  #######################################
  ##  Aggregate accross variant class
  ##  Required Outputs of detection
  #######################################
  
  CDx.mut_key = c("EGFR.L858R", "EGFR.T790M", "EGFR.p.Glu746_Ala750del")
  detection.snv.aggregated = detection.snv %>% group_by(Input) %>% summarise(nReplicates=sum(nReplicates), 
                                                                             Level1 = sum(Level1), 
                                                                             Level2 = sum(Level2), 
                                                                             Level3 = sum(Level3),
                                                                             Level4 = sum(Level4), 
                                                                             Level5 = sum(Level5)) %>%
    mutate(mut_key="aggregated", Pool = NA, variant_type = "SNV")
  detection.snv.aggregated = detection.snv.aggregated[,colnames(detection.snv)]
  detection.indel.aggregated = detection.indel %>% group_by(Input) %>% summarise(nReplicates=sum(nReplicates), 
                                                                                 Level1 = sum(Level1), 
                                                                                 Level2 = sum(Level2), 
                                                                                 Level3 = sum(Level3),
                                                                                 Level4 = sum(Level4), 
                                                                                 Level5 = sum(Level5)) %>%
    mutate(mut_key="aggregated", Pool = NA, variant_type = "INDEL")
  detection.indel.aggregated = detection.indel.aggregated[,colnames(detection.indel)]
  
  detection.final = rbind(detection.snv %>% filter(mut_key %in% CDx.mut_key),
                          detection.indel %>% filter(mut_key %in% CDx.mut_key), 
                          detection.snv.aggregated, 
                          detection.indel.aggregated, 
                          detection.cnv,
                          detection.fusion) %>%  arrange(desc(Input))
  
  detection.final
  write.table(detection.final, file=file.path(resdir, "Table4.detection.G360CDx.74genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
  detection.rate.final = detection.final 
  for(n in paste0("Level", 1:5)) detection.rate.final[[n]] = round(detection.rate.final[[n]]/detection.rate.final$nReplicates, 4)
  detection.rate.final
  write.table(detection.rate.final, file=file.path(resdir, "Table4.detection.rates.G360CDx.74genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  
    
  