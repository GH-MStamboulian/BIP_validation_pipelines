###############################
##  Table 5: LoD tables
###############################

  load(file.path(resdir, "work.RData"))
  load(file=file.path(resdir, "Table7.targeted.maf.G360CDx.RData"))  ##  targeted.maf.new.v3

  ################################################
  ##  TRUE MAF = G360CDx and 
  ##  detection == Guardant360 CDx for 55 genes
  ##  i.e. BIP v3.5.3 with rm_reportable = 1
  ##       but CNV call %in% 2:3
  ##  class levels LoD for SNVs and INDELs
  ##  skips filtering on rm_reportable for all version of BIP
  ################################################

  glm.models = list(SNV = list("5 ng" = list(), "30 ng" = list()),
                    INDEL = list("5 ng" = list(), "30 ng" = list()),
                    CNV = list("5 ng" = list(), "30 ng" = list()),
                    FUSION = list("5 ng" = list(), "30 ng" = list()))

  variants.df = list(SNV = list("5 ng" = list(), "30 ng" = list()),
                     INDEL = list("5 ng" = list(), "30 ng" = list()),
                     CNV = list("5 ng" = list(), "30 ng" = list()),
                     FUSION = list("5 ng" = list(), "30 ng" = list()))

  lod.table = NULL
  
  
  
  #########################
  ##  TRUE MAF = G360CDX
  #########################
  
  targeted.maf = targeted.maf.new.v3
  
  
  ############################################
  ##  SNV LoD: Aggregated LoD from all vars
  ############################################
  
  template = type %>% 
    filter(!is.na(Dilution)) %>%
    dplyr::select(c("run_sample_id", "Dilution", "ng_seq", "Type"))
  
  mutation2plot = unique(snv_table.filter$mutation)
  pool2plot = snv_table.filter$Type[match(mutation2plot, snv_table.filter$mutation)]
  
  for(ng in c("5 ng", "30 ng")) {
    dl = lapply(mutation2plot, function(mut) 
    {pool = pool2plot[which(mut == mutation2plot)];
    expected.maf = targeted.maf %>% filter(variant_type == "SNV" & ng_seq == ng & mutation == mut);
    du = snv_table.filter %>%
    #if ("rm_reportable" %in% colnames(du)) {
    #  du <- du %>% filter(rm_reportable == 1)
    #} else {
    #  du <- du %>% filter(call == 1)
    #}  
    #  du <- du %>%
      dplyr::filter(mutation == mut & Type == pool & ng_seq == ng) %>% 
      dplyr::select("run_sample_id", "Dilution", "ng_seq", "target_maf", "call", "mutation") %>%
      dplyr::right_join(template[template$Type == pool & template$ng_seq == ng,], 
                        by = c("run_sample_id", "Dilution", "ng_seq")) %>% 
      arrange(Dilution);
    du$mutation = mut;
    target_maf2plot = expected.maf[expected.maf$ng_seq == ng,] %>% 
      dplyr::select("ng_seq", "Dilution", "target_maf", "mutation") %>% 
      dplyr::right_join(du[,c("ng_seq", "Dilution", "call", "mutation", "Type")], 
                        by = c("ng_seq", "Dilution", "mutation"))
    iMissing = which(is.na(du$target_maf) );
    du$call = 1;
    du$call[iMissing] = 0;
    du$target_maf[iMissing] = target_maf2plot$target_maf[iMissing];
    
    #cat(paste(mut, sum(is.na(du$target_maf)), "\n"));
    du})
    dup = as.data.frame(do.call(rbind, dl))
    #table(dup$Dilution, dup$call)
    variants.df$SNV[[ng]][["aggregated"]] = dup
    
    m = glm(call ~ target_maf, data = dup,  family = quasibinomial(link = "probit"))
    #summary(m)
    glm.models$SNV[[ng]][["aggregated"]] = m
    
    newZ = data.frame(target_maf = seq(0, 0.02, 0.00001));
    #newZ = data.frame(target_maf = seq(2, 3, 0.001));
    newZ$pdet <- predict(m, newdata = newZ, type = "response")
    newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
    newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
    newZ$lowerL = newZ$pdet - newZ$se * 2
    
    
    tab = table(dup$target_maf, dup$call)
    #tab
    tab.df = data.frame(maf = as.numeric(rownames(tab)), sensitivity = tab[,2] / rowSums(tab))
    iLoD = max(which(newZ$pdet < 0.95))
    lod.probit = newZ$target_maf[iLoD]
    lod.FINAL = round(lod.probit, 3)
    dd = detection.rate.final %>% filter(mut_key == "aggregated" & variant_type == "SNV" & Input == ng)
    lvls = paste0("Level", 1:5)
    nValid = sum(dd[,lvls] >= 0.1 & dd[,lvls] <= 0.9)
    level_maf = NA
    lod.l95 = pmax(0, newZ$pdet[iLoD] - newZ$se[iLoD] * 2)
    lod.u95 = pmin(newZ$pdet[iLoD] + newZ$se[iLoD] * 2, 1)
    
    lod.table = rbind( lod.table, c(variant = "aggregated", variant_type="SNV", Input=ng, 
                                    LoD=lod.FINAL, nValid = nValid, level_maf = level_maf,
                                    lod.probit = lod.probit, detection.CI.low=lod.l95, detection.CI.upper=lod.u95))
  }
  
  ############################################
  ##  INDEL LoD: Aggregated LoD from all vars
  ############################################
  
  mutation2plot = unique(indel_table.filter$mutation)
  pool2plot = indel_table.filter$Type[match(mutation2plot, indel_table.filter$mutation)]
  lod.plots = list()
  
  
  for(ng in c("5 ng", "30 ng")) {
    dl = lapply(mutation2plot, function(mut) 
    {pool = pool2plot[which(mut == mutation2plot)];
    expected.maf = targeted.maf %>% filter(variant_type == "INDEL" & ng_seq == ng & mutation == mut);
    du = indel_table.filter %>%
    #if ("rm_reportable" %in% colnames(du)) {
    #  du <- du %>% filter(rm_reportable == 1)
    #} else {
    #  du <- du %>% filter(call == 1)
    #}  
    #du <- du %>%
      dplyr::filter(mutation == mut & Type == pool & ng_seq == ng) %>% 
      dplyr::select("run_sample_id", "Dilution", "ng_seq", "target_maf", "call", "mutation") %>%
      dplyr::right_join(template[template$Type == pool & template$ng_seq == ng,], 
                        by = c("run_sample_id", "Dilution", "ng_seq")) %>% 
      arrange(Dilution);
    du$mutation = mut;
    target_maf2plot = expected.maf[expected.maf$ng_seq == ng,] %>% 
      dplyr::select("ng_seq", "Dilution", "target_maf", "mutation") %>% 
      dplyr::right_join(du[,c("ng_seq", "Dilution", "call", "mutation", "Type")], 
                        by = c("ng_seq", "Dilution", "mutation"))
    iMissing = which(is.na(du$target_maf) );
    du$call = 1;
    du$call[iMissing] = 0;
    du$target_maf[iMissing] = target_maf2plot$target_maf[iMissing];
    
    #cat(paste(mut, sum(is.na(du$target_maf)), "\n"));
    du})
    dup = as.data.frame(do.call(rbind, dl))
    #table(dup$Dilution, dup$call)
    variants.df$INDEL[[ng]][["aggregated"]] = dup
    
    m = glm(call ~ target_maf, data = dup,  family = quasibinomial(link = "probit"))
    #summary(m)
    glm.models$INDEL[[ng]][["aggregated"]] = m
    
    newZ = data.frame(target_maf = seq(0, max(dup$target_maf, na.rm=T), 0.00001));
    #newZ = data.frame(target_maf = seq(2, 3, 0.001));
    newZ$pdet <- predict(m, newdata = newZ, type = "response")
    newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
    newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
    newZ$lowerL = newZ$pdet - newZ$se * 2
    
    
    tab = table(dup$target_maf, dup$call)
    #tab
    tab.df = data.frame(maf = as.numeric(rownames(tab)), sensitivity = tab[,2] / rowSums(tab))
    iLoD = max(which(newZ$pdet < 0.95))
    lod.probit = newZ$target_maf[iLoD]
    lod.FINAL = round(lod.probit, 3)
    dd = detection.rate.final %>% filter(mut_key == "aggregated" & variant_type == "INDEL" & Input == ng)
    lvls = paste0("Level", 1:5)
    nValid = sum(dd[,lvls] >= 0.1 & dd[,lvls] <= 0.9)
    level_maf = NA
    lod.l95 = pmax(0, newZ$pdet[iLoD] - newZ$se[iLoD] * 2)
    lod.u95 = pmin(newZ$pdet[iLoD] + newZ$se[iLoD] * 2, 1)
    
    lod.table = rbind( lod.table, c(variant = "aggregated", variant_type="INDEL", Input=ng, 
                                    LoD=lod.FINAL, nValid = nValid, level_maf = level_maf, 
                                    lod.probit = lod.probit, detection.CI.low=lod.l95, detection.CI.upper=lod.u95))
  }
  
  ###########################################
  ##  CNA LoD: individual genes
  ###########################################
  
  mutation2plot = unique(cna_table.filter$mutation)
  pools2plot = cna_table$Type[match(mutation2plot, cna_table.filter$mutation)]
  ng_seq2plot = unique(cna_table.filter$ng_seq)
  
  for(mut in mutation2plot) {
    for(ng in ng_seq2plot) {
      
      pool = pools2plot[mut == mutation2plot]
      d = cna_table.filter %>% filter(mutation == mut & ng_seq == ng)
      variants.df$CNV[[ng]][[mut]] = d
      
      m = glm(call ~ target_maf, data = d,  family = quasibinomial(link = "probit"))
      #summary(m)
      glm.models$CNV[[ng]][[mut]] = m
      
      #newZ = data.frame(titration = seq(0, 1, 0.001));
      newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = pmax(0, newZ$pdet - newZ$se * 2)
      
      iLoD = max(which(newZ$pdet < 0.95))
      lod.probit = newZ$target_maf[iLoD]
      lod.FINAL = ceiling(lod.probit * 10)/10
      level_maf = NA
      nValid = NA
      if(mut %in% detection.rate.final$mut_key) {
        dd = detection.rate.final %>% filter(mut_key == mut & Input == ng)
        lvls = paste0("Level", 1:5)
        nValid = sum(dd[,lvls] >= 0.1 & dd[,lvls] <= 0.9)
        if(nValid < 3) {
          i = min(which(dd[,lvls] >= 0.95))
          level = lvls[i]
          gene = mut
          level_maf = targeted.maf$target_maf[targeted.maf$Gene == gene & targeted.maf$variant_type == "CNV" &
                                                targeted.maf$Dilution == level & targeted.maf$ng_seq == ng]
          lod.FINAL = ceiling(level_maf * 10)/10
        }
      }
      lod.l95 = pmax(0, newZ$lowerL[iLoD])
      lod.u95 = pmin(newZ$upperL[iLoD], 1)
      
      lod.table = rbind( lod.table, c(variant = mut, variant_type="CNV", Input=ng, 
                                      LoD=lod.FINAL, nValid = nValid, level_maf = level_maf, 
                                      lod.probit = lod.probit, detection.CI.low=lod.l95, detection.CI.upper=lod.u95))
    }
  }
  
  
  ###########################################
  ##  FUSION LoD: individual genes
  ##  x-axis expected MAF
  ###########################################
  
  template = type %>% 
    filter(!is.na(Dilution)) %>%
    dplyr::select(c("run_sample_id", "Dilution", "ng_seq", "Type"))
  mutation2plot = unique(fusion_table.filter$mutation)
  dilution2plot = sort(unique(type$Dilution))
  pools2plot = fusion_table.filter$Type[match(mutation2plot, fusion_table.filter$mutation)]
  ng_seq2plot = unique(fusion_table.filter$ng_seq)
  
  for(mut in mutation2plot) {
    for(ng in ng_seq2plot) {
      expected.maf = targeted.maf %>% filter(variant_type == "FUSION" & ng_seq == ng)
      pool = pools2plot[mut == mutation2plot]
      
      totals = as.vector(table(type$Dilution[type$Type == pool & type$ng_seq == ng]))
      d = fusion_table.filter %>% filter(mutation == mut & ng_seq == ng) %>% 
        dplyr::select("run_sample_id", "Dilution", "ng_seq", "target_maf", "call", "mutation") %>%
        right_join(template[template$Type == pool & template$ng_seq == ng,], 
                   by = c("run_sample_id", "Dilution", "ng_seq")) %>% 
        arrange(Dilution)
      d$mutation = mut
      
      target_maf2plot = expected.maf[expected.maf$ng_seq == ng,] %>% 
        dplyr::select("ng_seq", "Dilution", "target_maf", "mutation") %>% 
        dplyr::right_join(d[,c("ng_seq", "Dilution", "call", "mutation", "Type")], 
                          by = c("ng_seq", "Dilution", "mutation"))
      
      d$call = 1
      iMissing = which(is.na(d$target_maf) )
      d$target_maf[iMissing] = target_maf2plot$target_maf[iMissing];
      d$call[iMissing] = 0
      variants.df$FUSION[[ng]][[mut]] = d
      
      m = glm(call ~ target_maf, data = d,  family = quasibinomial(link = "probit"))
      #summary(m)
      glm.models$FUSION[[ng]][[mut]] = m
      
      newZ = data.frame(target_maf = seq(0, max(d$target_maf, na.rm=T), 0.00001));
      #newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = pmax(0, newZ$pdet - newZ$se * 2)
      
      iLoD = max(which(newZ$pdet < 0.95))
      lod.probit = newZ$target_maf[iLoD]
      lod.FINAL = round(lod.probit, 3)
      nValid = NA
      level_maf = NA
      if(mut %in% detection.rate.final$mut_key) {
        dd = detection.rate.final %>% filter(mut_key == mut & Input == ng)
        lvls = paste0("Level", 1:5)
        nValid = sum(dd[,lvls] >= 0.1 & dd[,lvls] <= 0.9)
        if(nValid < 3) {
          i = min(which(dd[,lvls] >= 0.95))
          level = lvls[i]
          level_maf = targeted.maf$target_maf[targeted.maf$mutation == mut & targeted.maf$Dilution == level & targeted.maf$ng_seq == ng]
          lod.FINAL = round(level_maf, 3)
        }
      }
      lod.l95 = newZ$lowerL[iLoD]
      lod.u95 =newZ$upperL[iLoD]
      
      lod.table = rbind( lod.table, c(variant = mut, variant_type="FUSION", Input=ng, 
                                      LoD=lod.FINAL, nValid = nValid, level_maf = level_maf,
                                      lod.probit = lod.probit, detection.CI.low=lod.l95, detection.CI.upper=lod.u95))
    }
  }
  
  lod.table.part1 = as.data.frame(lod.table, stringsAsFactors=F)
  lod.table.part1 = lod.table.part1 %>% arrange(desc(Input))
  lod.table.part1
  
  write.table(lod.table.part1, file=file.path(resdir, "Table5.LoD.part1.G360CDx.55genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  save(lod.table.part1, file=file.path(resdir, "Table5.LoD.part1.G360CDx.55genes.RData"))
  
  
  #####################################
  ##  Individual variants
  ##  SAME conditions as above
  ##  use glm.model list from above
  #####################################
  
  ################################################
  ##  TRUE MAF = G360CDx and 
  ##  detection == Guardant360 CDx for 55 genes
  ##  i.e. previous results with G360CDx method
  ##  class levels LoD for SNVs and INDELs
  ################################################
  
  ###############################
  ##  Initialize NEW LoD TABLE
  ###############################
  
  lod.table = NULL
  
  
  #####################################
  ##  SNV LoD: individual variants LoD
  #####################################
  
  template = type %>% 
    filter(!is.na(Dilution)) %>%
    dplyr::select(c("run_sample_id", "Dilution", "ng_seq", "Type"))
  
  mutation2plot = unique(snv_table.filter$mutation)
  pool2plot = snv_table.filter$Type[match(mutation2plot, snv_table.filter$mutation)]
  
  for(mut in mutation2plot) {
    for(ng in c("5 ng", "30 ng")) {
      
      pool = pool2plot[which(mut == mutation2plot)];
      expected.maf = targeted.maf %>% filter(variant_type == "SNV" & ng_seq == ng & mutation == mut);
      du = snv_table.filter %>%
      #if ("rm_reportable" %in% colnames(du)) {
      #  du <- du %>% filter(rm_reportable == 1)
      #} else {
      #  du <- du %>% filter(call == 1)
      #}  
      # du <- du %>%
        dplyr::filter(mutation == mut & Type == pool & ng_seq == ng) %>% 
        dplyr::select("run_sample_id", "Dilution", "ng_seq", "target_maf", "call", "mutation") %>%
        dplyr::right_join(template[template$Type == pool & template$ng_seq == ng,], 
                          by = c("run_sample_id", "Dilution", "ng_seq")) %>% 
        arrange(Dilution);
      du$mutation = mut;
      target_maf2plot = expected.maf[expected.maf$ng_seq == ng,] %>% 
        dplyr::select("ng_seq", "Dilution", "target_maf", "mutation") %>% 
        dplyr::right_join(du[,c("ng_seq", "Dilution", "call", "mutation", "Type")], 
                          by = c("ng_seq", "Dilution", "mutation"))
      iMissing = which(is.na(du$target_maf) );
      du$call = 1;
      du$call[iMissing] = 0;
      du$target_maf[iMissing] = target_maf2plot$target_maf[iMissing];
      
      #cat(paste(mut, sum(is.na(du$target_maf)), "\n"));
      dup = du
      variants.df$SNV[[ng]][[mut]] = dup
      
      m = glm(call ~ target_maf, data = dup,  family = quasibinomial(link = "probit"))
      #summary(m)
      glm.models$SNV[[ng]][[mut]] = m
      
      newZ = data.frame(target_maf = seq(0, 0.02, 0.00001));
      #newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = newZ$pdet - newZ$se * 2
      
      
      tab = table(dup$target_maf, dup$call)
      #tab
      tab.df = data.frame(maf = as.numeric(rownames(tab)), sensitivity = tab[,2] / rowSums(tab))
      iLoD = max(which(newZ$pdet < 0.95))
      lod.probit = newZ$target_maf[iLoD]
      lod.FINAL = round(lod.probit, 3)
      level_maf = NA
      nValid = NA
      #if(mut %in% detection.rate.final$mut_key) {
      dd = detection.rate.all %>% filter(mut_key == mut & variant_type == "SNV" & Input == ng)
      lvls = paste0("Level", 1:5)
      nValid = sum(dd[,lvls] >= 0.1 & dd[,lvls] <= 0.9)
      ##  If there are less than 3 dilution levels with detection in [10%,90%]
      ##  use TRUE MAF of first dilution with detection >= 95% as LoD 
      if(nValid < 3) {
        i = min(which(dd[,lvls] >= 0.95))
        level = lvls[i]
        gene = strsplit(mut, "\\.")[[1]][1]
        mut_aa = paste(strsplit(mut, "\\.")[[1]][-1], collapse=".")
        level_maf = targeted.maf$target_maf[targeted.maf$Gene == gene & targeted.maf$mut_aa == mut_aa & targeted.maf$Dilution == level & targeted.maf$ng_seq == ng]
        lod.FINAL = round(level_maf, 3)
      }
      #}
      lod.l95 = pmax(0, newZ$pdet[iLoD] - newZ$se[iLoD] * 2)
      lod.u95 = pmin(newZ$pdet[iLoD] + newZ$se[iLoD] * 2, 1)
      
      lod.table = rbind( lod.table, c(variant = mut, variant_type="SNV", Input=ng, 
                                      LoD=lod.FINAL, nValid = nValid, level_maf = level_maf,
                                      lod.probit = lod.probit, detection.CI.low=lod.l95, detection.CI.upper=lod.u95))
    }
  }
  
  ############################################
  ##  INDEL LoD: individual variants LoD
  ############################################
  
  mutation2plot = unique(indel_table.filter$mutation)
  pool2plot = indel_table.filter$Type[match(mutation2plot, indel_table.filter$mutation)]
  
  for(mut in mutation2plot) {
    for(ng in c("5 ng", "30 ng")) {
      pool = pool2plot[which(mut == mutation2plot)];
      expected.maf = targeted.maf %>% filter(variant_type == "INDEL" & ng_seq == ng & mutation == mut);
      du = indel_table.filter %>%
      #if ("rm_reportable" %in% colnames(du)) {
      #  du <- du %>% filter(rm_reportable == 1)
      #} else {
      #  du <- du %>% filter(call == 1)
      #}  
      #  du <- du %>%
        dplyr::filter(mutation == mut & Type == pool & ng_seq == ng) %>% 
        dplyr::select("run_sample_id", "Dilution", "ng_seq", "target_maf", "call", "mutation") %>%
        dplyr::right_join(template[template$Type == pool & template$ng_seq == ng,], 
                          by = c("run_sample_id", "Dilution", "ng_seq")) %>% 
        arrange(Dilution);
      du$mutation = mut;
      target_maf2plot = expected.maf[expected.maf$ng_seq == ng,] %>% 
        dplyr::select("ng_seq", "Dilution", "target_maf", "mutation") %>% 
        dplyr::right_join(du[,c("ng_seq", "Dilution", "call", "mutation", "Type")], 
                          by = c("ng_seq", "Dilution", "mutation"))
      iMissing = which(is.na(du$target_maf) );
      du$call = 1;
      du$call[iMissing] = 0;
      du$target_maf[iMissing] = target_maf2plot$target_maf[iMissing];
      
      #cat(paste(mut, sum(is.na(du$target_maf)), "\n"));
      du
      dup = du
      variants.df$INDEL[[ng]][[mut]] = dup
      
      m = glm(call ~ target_maf, data = dup,  family = quasibinomial(link = "probit"))
      #summary(m)
      glm.models$INDEL[[ng]][[mut]] = m
      
      newZ = data.frame(target_maf = seq(0, max(dup$target_maf, na.rm=T), 0.00001));
      #newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = newZ$pdet - newZ$se * 2
      
      
      tab = table(dup$target_maf, dup$call)
      #tab
      tab.df = data.frame(maf = as.numeric(rownames(tab)), sensitivity = tab[,2] / rowSums(tab))
      iLoD = max(which(newZ$pdet < 0.95))
      lod.probit = newZ$target_maf[iLoD]
      lod.FINAL = round(lod.probit, 3)
      nValid = NA
      level_maf = NA
      #if(mut %in% detection.rate.final$mut_key) {
      dd = detection.rate.all %>% filter(mut_key == mut & variant_type == "INDEL" & Input == ng)
      lvls = paste0("Level", 1:5)
      nValid = sum(dd[,lvls] >= 0.1 & dd[,lvls] <= 0.9)
      ##  If there are less than 3 dilution levels with detection in [10%,90%]
      ##  use TRUE MAF of first dilution with detection >= 95% as LoD 
      if(nValid < 3) {
        i = min(which(dd[,lvls] >= 0.95))
        level = lvls[i]
        gene = strsplit(mut, "\\.")[[1]][1]
        mut_aa = paste(strsplit(mut, "\\.")[[1]][-1], collapse=".")
        level_maf = targeted.maf$target_maf[targeted.maf$Gene == gene & targeted.maf$mut_aa == mut_aa & targeted.maf$Dilution == level & targeted.maf$ng_seq == ng]
        lod.FINAL = round(level_maf, 3)
      }
      #}
      lod.l95 = pmax(0, newZ$pdet[iLoD] - newZ$se[iLoD] * 2)
      lod.u95 = pmin(newZ$pdet[iLoD] + newZ$se[iLoD] * 2, 1)
      
      lod.table = rbind( lod.table, c(variant = mut, variant_type="INDEL", Input=ng, 
                                      LoD=lod.FINAL, nValid = nValid, level_maf = level_maf,
                                      lod.probit = lod.probit, detection.CI.low=lod.l95, detection.CI.upper=lod.u95))
      
    }
  }
  
  ###########################################
  ##  CNA LoD: individual genes LoD
  ###########################################
  
  mutation2plot = unique(cna_table.filter$mutation)
  pools2plot = cna_table.filter$Type[match(mutation2plot, cna_table.filter$mutation)]
  ng_seq2plot = unique(cna_table.filter$ng_seq)
  
  for(mut in mutation2plot) {
    for(ng in ng_seq2plot) {
      
      pool = pools2plot[mut == mutation2plot]
      d = cna_table.filter %>% filter(mutation == mut & ng_seq == ng)
      m = glm(call ~ target_maf, data = d,  family = quasibinomial(link = "probit"))
      #summary(m)
      
      #newZ = data.frame(titration = seq(0, 1, 0.001));
      newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = pmax(0, newZ$pdet - newZ$se * 2)
      
      iLoD = max(which(newZ$pdet < 0.95))
      lod.probit = newZ$target_maf[iLoD]
      lod.FINAL = ceiling(lod.probit * 10)/10
      nValid = NA
      level_maf = NA
      if(mut %in% detection.rate.final$mut_key) {
        dd = detection.rate.final %>% filter(mut_key == mut & Input == ng)
        lvls = paste0("Level", 1:5)
        nValid = sum(dd[,lvls] >= 0.1 & dd[,lvls] <= 0.9)
        
        ##  If there are less than 3 dilution levels with detection in [10%,90%]
        ##  use TRUE MAF of first dilution with detection >= 95% as LoD 
        if(nValid < 3) {
          i = min(which(dd[,lvls] >= 0.95))
          level = lvls[i]
          gene = mut
          level_maf = targeted.maf$target_maf[targeted.maf$Gene == gene & targeted.maf$variant_type == "CNV" & targeted.maf$Dilution == level & targeted.maf$ng_seq == ng]
          lod.FINAL = ceiling(level_maf * 10)/10
        }
      }
      lod.l95 = pmax(0, newZ$lowerL[iLoD])
      lod.u95 = pmin(newZ$upperL[iLoD], 1)
      
      lod.table = rbind( lod.table, c(variant = mut, variant_type="CNV", Input=ng, 
                                      LoD=lod.FINAL, nValid = nValid, level_maf = level_maf,
                                      lod.probit = lod.probit, detection.CI.low=lod.l95, detection.CI.upper=lod.u95))
    }
  }
  
  
  ###########################################
  ##  FUSION LoD: individual genes LoD
  ###########################################
  
  template = type %>% 
    filter(!is.na(Dilution)) %>%
    dplyr::select(c("run_sample_id", "Dilution", "ng_seq", "Type"))
  mutation2plot = unique(fusion_table.filter$mutation)
  dilution2plot = sort(unique(type$Dilution))
  pools2plot = fusion_table.filter$Type[match(mutation2plot, fusion_table.filter$mutation)]
  ng_seq2plot = unique(fusion_table.filter$ng_seq)
  
  for(mut in mutation2plot) {
    for(ng in ng_seq2plot) {
      expected.maf = targeted.maf %>% filter(variant_type == "FUSION" & ng_seq == ng)
      pool = pools2plot[mut == mutation2plot]
      
      totals = as.vector(table(type$Dilution[type$Type == pool & type$ng_seq == ng]))
      d = fusion_table %>% filter(mutation == mut & ng_seq == ng) %>% 
        dplyr::select("run_sample_id", "Dilution", "ng_seq", "target_maf", "call", "mutation") %>%
        right_join(template[template$Type == pool & template$ng_seq == ng,], 
                   by = c("run_sample_id", "Dilution", "ng_seq")) %>% 
        arrange(Dilution)
      d$mutation = mut
      
      target_maf2plot = expected.maf[expected.maf$ng_seq == ng,] %>% 
        dplyr::select("ng_seq", "Dilution", "target_maf", "mutation") %>% 
        dplyr::right_join(d[,c("ng_seq", "Dilution", "call", "mutation", "Type")], 
                          by = c("ng_seq", "Dilution", "mutation"))
      
      d$call = 1
      iMissing = which(is.na(d$target_maf) )
      d$target_maf[iMissing] = target_maf2plot$target_maf[iMissing];
      d$call[iMissing] = 0
      
      m = glm(call ~ target_maf, data = d,  family = quasibinomial(link = "probit"))
      summary(m)
      
      newZ = data.frame(target_maf = seq(0, max(d$target_maf, na.rm=T), 0.00001));
      #newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = newZ$pdet - newZ$se * 2
      
      iLoD = max(which(newZ$pdet < 0.95))
      lod.probit = newZ$target_maf[iLoD]
      lod.FINAL = round(lod.probit, 3)
      nValid = NA
      level_maf = NA
      if(mut %in% detection.rate.final$mut_key) {
        dd = detection.rate.final %>% filter(mut_key == mut & Input == ng)
        lvls = paste0("Level", 1:5)
        nValid = sum(dd[,lvls] >= 0.1 & dd[,lvls] <= 0.9)
        
        ##  If there are less than 3 dilution levels with detection in [10%,90%]
        ##  use TRUE MAF of first dilution with detection >= 95% as LoD 
        if(nValid < 3) {
          i = min(which(dd[,lvls] >= 0.95))
          level = lvls[i]
          level_maf = targeted.maf$target_maf[targeted.maf$mutation == mut & targeted.maf$Dilution == level & targeted.maf$ng_seq == ng]
          lod.FINAL = round(level_maf, 3)
        }
      }
      lod.l95 = pmax(0, newZ$lowerL[iLoD])
      lod.u95 = pmin(newZ$upperL[iLoD], 1)
      
      lod.table = rbind( lod.table, c(variant = mut, variant_type="FUSION", Input=ng, 
                                      LoD=lod.FINAL, nValid = nValid, level_maf = level_maf,
                                      lod.probit = lod.probit, detection.CI.low=lod.l95, detection.CI.upper=lod.u95))
    }
  }
  
  lod.table.part2 = as.data.frame(lod.table, stringsAsFactors=F)
  lod.table.part2 = lod.table.part2 %>% arrange(desc(Input))
  #lod.table.part2
  write.table(lod.table.part2, file=file.path(resdir, "Table5.LoD.part2.G360CDx.55genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  save(lod.table.part2, file=file.path(resdir, "Table5.LoD.part2.G360CDx.55genes.RData"))
  
  ######################
  ##  TABLE 5 Report
  ######################
  
  lod.table.final = rbind(lod.table.part2 %>% filter(variant %in% CDx.mut_key),
                          lod.table.part1) %>% arrange(desc(Input))
  #lod.table.final$LoD = round(as.numeric(lod.table.final$LoD), 3)
  #lod.table.final$LoD[ lod.table.final$variant_type != "CNV"] = 100*lod.table.final$LoD[ lod.table.final$variant_type != "CNV"]
  lod.table.final
  write.table(lod.table.final, file=file.path(resdir, "Table5.LoD.G360CDx.55genes.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  save(lod.table.final, file=file.path(resdir, "Table5.LoD.G360CDx.55genes.RData"))
  
  save(glm.models, variants.df, file=file.path(resdir, "glm.models.LoD.G360CDx.55genes.200128.RData"))
  
  ###########################
  ##  FIGURE 2: LoD plots
  ###########################
  
  ## SNV LoD
  
  #mutation2plot = c("aggregated", intersect(names(glm.models$SNV[['5 ng']]), CDx.mut_key))
  mutation2plot = names(glm.models$SNV[['5 ng']])
  ng_seq2plot = c("5 ng", "30 ng")
  
  pdf(file=file.path(resdir, "FIGURE2.LoD_PLOTS_SNV.G360CDx.55genes.v1.pdf"), height=7, width=14)
  op = par(mfrow=c(1,2))
  for(mut in mutation2plot) {
    for(ng in ng_seq2plot) {
      
      dup = variants.df$SNV[[ng]][[mut]]
      tab = table(dup$target_maf, dup$call)
      #tab
      
      m = glm.models$SNV[[ng]][[mut]]
      newZ = data.frame(target_maf = seq(0, 0.04, 0.00001));
      #newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = newZ$pdet - newZ$se * 2
      
      tab.df = data.frame(maf = as.numeric(rownames(tab)), sensitivity = tab[,2] / rowSums(tab))
      lod = newZ$target_maf[max(which(newZ$pdet < 0.95))]
      
      xlim = c(0, max(dup %>% filter(ng_seq == ng) %>% dplyr::select(target_maf)))
      if(ng == "5 ng") xlim[2] = pmax(xlim[2], 0.04)
      if(ng == "30 ng") xlim[2] = pmax(xlim[2], 0.012)
      xlim = 100*xlim
      
      plot(100*newZ$target_maf, newZ$pdet, pch = ".", cex = 4, col = 4, font = 2, ylim = c(0, 1), xlim=xlim,
           xlab = "Targeted MAF (%)", ylab = "Probability of detection", main = paste("LoD SNV", mut, ",", ng, "samples")); 
      grid()
      points(100*newZ$target_maf, newZ$pdet - newZ$se * 2, pch = ".", cex = 2, col = 2); grid()
      points(100*newZ$target_maf, pmin(1.0, newZ$pdet + newZ$se * 2), pch = ".", cex = 2, col = 2); grid()
      
      points(100*as.numeric(rownames(tab)), tab[,2] / rowSums(tab), pch = 19, cex=2, col = 3)
      points(100*as.numeric(rownames(tab)), tab[,2] / rowSums(tab), pch = 19, cex=2, col = 3)
      abline(h = c(0.1, 0.8, 0.9, 0.95), col = "gray", lwd = 3, lty = 2)
      abline(v = 100*lod, lwd=2, col=3, lty=2)
      text(100*lod, 0.2, paste("LoD =", round(100*lod, 2), "%"), cex=1.5, pos=4)
    }
  }
  par(op)
  dev.off()
  
  ## INDEL LoD
  
  #mutation2plot = c("aggregated", intersect(names(glm.models$INDEL[['5 ng']]), CDx.mut_key))
  mutation2plot = names(glm.models$INDEL[['5 ng']])
  ng_seq2plot = c("5 ng", "30 ng")
  
  pdf(file=file.path(resdir, "FIGURE2.LoD_PLOTS_INDEL.G360CDx.55genes.v1.pdf"), height=7, width=14)
  op = par(mfrow=c(1,2))
  for(mut in mutation2plot) {
    for(ng in ng_seq2plot) {
      
      dup = variants.df$INDEL[[ng]][[mut]]
      tab = table(dup$target_maf, dup$call)
      #tab
      
      m = glm.models$INDEL[[ng]][[mut]]
      newZ = data.frame(target_maf = seq(0, 0.04, 0.00001));
      #newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = newZ$pdet - newZ$se * 2
      
      tab.df = data.frame(maf = as.numeric(rownames(tab)), sensitivity = tab[,2] / rowSums(tab))
      lod = newZ$target_maf[max(which(newZ$pdet < 0.95))]
      
      xlim = c(0, max(dup %>% filter(ng_seq == ng) %>% dplyr::select(target_maf)))
      if(ng == "5 ng") xlim[2] = pmax(xlim[2], 0.04)
      if(ng == "30 ng") xlim[2] = pmax(xlim[2], 0.012)
      xlim = 100 * xlim
      
      plot(100 * newZ$target_maf, newZ$pdet, pch = ".", cex = 4, col = 4, font = 2, ylim = c(0, 1),xlim = xlim,
           xlab = "Targeted MAF (%)", ylab = "Probability of detection", main = paste("LoD INDEL", mut, ",", ng, "samples")); 
      grid()
      points(100 * newZ$target_maf, newZ$pdet - newZ$se * 2, pch = ".", cex = 2, col = 2); grid()
      points(100 * newZ$target_maf, pmin(1.0, newZ$pdet + newZ$se * 2), pch = ".", cex = 2, col = 2); grid()
      
      points(100 * as.numeric(rownames(tab)), tab[,2] / rowSums(tab), pch = 19, cex=2, col = 3)
      points(100 * as.numeric(rownames(tab)), tab[,2] / rowSums(tab), pch = 19, cex=2, col = 3)
      abline(h = c(0.1, 0.8, 0.9, 0.95), col = "gray", lwd = 3, lty = 2)
      abline(v = 100 * lod, lwd=2, col=3, lty=2)
      text(100 * lod, 0.2, paste("LoD =", round(100*lod, 2), "%"), cex=1.5, pos=4)
    }
  }
  par(op)
  dev.off()
  
  ## CNA LoD
  
  mutation2plot = names(glm.models$CNV[['5 ng']])
  ng_seq2plot = c("5 ng", "30 ng")
  
  pdf(file=file.path(resdir, "FIGURE2.LoD_PLOTS_CNV.G360CDx.55genes.v1.pdf"), height=7, width=14)
  op = par(mfrow=c(1,2))
  for(mut in mutation2plot) {
    for(ng in ng_seq2plot) {
      
      dup = variants.df$CNV[[ng]][[mut]]
      tab = table(dup$target_maf, dup$call)
      #tab
      
      m = glm.models$CNV[[ng]][[mut]]
      newZ = data.frame(target_maf = seq(2, 3, 0.0001));
      #newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = newZ$pdet - newZ$se * 2
      
      tab.df = data.frame(maf = as.numeric(rownames(tab)), sensitivity = tab[,2] / rowSums(tab))
      lod = newZ$target_maf[max(which(newZ$pdet < 0.95))]
      
      plot(newZ$target_maf, newZ$pdet, pch = ".", cex = 4, col = 4, font = 2, ylim = c(0, 1),
           xlab = "Targeted MAF", ylab = "Probability of detection", main = paste("LoD CNV", mut, ",", ng, "samples")); 
      grid()
      points(newZ$target_maf, pmax(0, newZ$pdet - 2*newZ$se), pch = ".", cex = 2, col = 2); grid()
      points(newZ$target_maf, pmin(1.0, newZ$pdet + 2*newZ$se), pch = ".", cex = 2, col = 2); grid()
      
      points(rownames(tab), tab[,2] / rowSums(tab), pch = 19, cex=2, col = 3)
      abline(h = c(0.1, 0.8, 0.9, 0.95), col = "gray", lwd = 3, lty = 2)
      abline(v = lod, lwd=2, col=3, lty=2)
      text(lod, 0.2, paste("LoD =", round(lod, 2)), cex=1.5, pos=4)
    }
  }
  par(op)
  dev.off()
  
  ## FUSION LoD
  
  mutation2plot = names(glm.models$FUSION[['5 ng']])
  ng_seq2plot = c("5 ng", "30 ng")
  
  pdf(file=file.path(resdir, "FIGURE2.LoD_PLOTS_FUSION.G360CDx.55genes.v1.pdf"), height=7, width=14)
  op = par(mfrow=c(1,2))
  for(mut in mutation2plot) {
    for(ng in ng_seq2plot) {
      
      dup = variants.df$FUSION[[ng]][[mut]]
      tab = table(dup$target_maf, dup$call)
      #tab
      
      m = glm.models$FUSION[[ng]][[mut]]
      newZ = data.frame(target_maf = seq(0, 0.04, 0.00001));
      #newZ = data.frame(target_maf = seq(2, 3, 0.001));
      newZ$pdet <- predict(m, newdata = newZ, type = "response")
      newZ$se <- predict(m, newdata = newZ, type = "response", se.fit=TRUE)$se.fit
      newZ$upperL = pmin(1.0, newZ$pdet + newZ$se * 2)
      newZ$lowerL = newZ$pdet - newZ$se * 2
      
      tab.df = data.frame(maf = as.numeric(rownames(tab)), sensitivity = tab[,2] / rowSums(tab))
      lod = newZ$target_maf[max(which(newZ$pdet < 0.95))]
      
      xlim = c(0, max(dup %>% filter(ng_seq == ng) %>% dplyr::select(target_maf)))
      if(ng == "5 ng") xlim[2] = pmax(xlim[2], 0.04)
      if(ng == "30 ng") xlim[2] = pmax(xlim[2], 0.012)
      xlim = 100 * xlim
      
      plot(100 * newZ$target_maf, newZ$pdet, pch = ".", cex = 4, col = 4, font = 2, ylim = c(0, 1), xlim=xlim,
           xlab = "Targeted MAF (%)", ylab = "Probability of detection", main = paste("LoD FUSION", mut, ",", ng, "samples")); 
      grid()
      points(100 * newZ$target_maf, newZ$pdet - newZ$se * 2, pch = ".", cex = 2, col = 2); grid()
      points(100 * newZ$target_maf, pmin(1.0, newZ$pdet + newZ$se * 2), pch = ".", cex = 2, col = 2); grid()
      
      points(100 * as.numeric(rownames(tab)), tab[,2] / rowSums(tab), pch = 19, cex=2, col = 3)
      points(100 * as.numeric(rownames(tab)), tab[,2] / rowSums(tab), pch = 19, cex=2, col = 3)
      abline(h = c(0.1, 0.8, 0.9, 0.95), col = "gray", lwd = 3, lty = 2)
      abline(v = 100 * lod, lwd=2, col=3, lty=2)
      text(100 * lod, 0.2, paste("LoD =", round(100*lod, 2), "%"), cex=1.5, pos=4)
    }
  }
  par(op)
  dev.off()
  
  
  ##########################
  ##  LoD table all genes
  ##########################
  
  #load(file=file.path(resdir, "glm.models.LoD.G360CDx.55genes.200128.RData"))
  
  lvls = paste0("Level", 1:5)
  nTotals = c(14, 20)
  names(nTotals) = c("30 ng", "5 ng")
  
  ##  SNVs
  tab.snv = read.table(file=file.path(resdir, "Table5.LoD.part2.G360CDx.55genes.tsv"), sep="\t", header=T, stringsAsFactors = F) %>%
    filter(variant_type == "SNV") %>%
    mutate(mut_key = variant, lod.maf = lod.probit,   lod.FINAL = LoD) %>%  
    dplyr::select("mut_key", "Input", "lod.probit", "nValid", "level_maf", "lod.maf", "lod.FINAL")
  tab.snv.agg = read.table(file=file.path(resdir, "Table5.LoD.part1.G360CDx.55genes.tsv"), sep="\t", header=T, stringsAsFactors = F) %>%
    filter(variant_type == "SNV") %>%
    mutate(mut_key = variant, lod.maf = lod.probit,   lod.FINAL = LoD) %>%  
    dplyr::select("mut_key", "Input", "lod.probit", "nValid", "level_maf", "lod.maf", "lod.FINAL")
  tab.snv = rbind(tab.snv, tab.snv.agg) %>%
    arrange(Input)
  
  ##  INDELS
  tab.indel = read.table(file=file.path(resdir, "Table5.LoD.part2.G360CDx.55genes.tsv"), sep="\t", header=T, stringsAsFactors = F) %>%
    filter(variant_type == "INDEL") %>%
    mutate(mut_key = variant, lod.maf = lod.probit,   lod.FINAL = LoD) %>%  
    dplyr::select("mut_key", "Input", "lod.probit", "nValid", "level_maf", "lod.maf", "lod.FINAL")
  tab.indel.agg = read.table(file=file.path(resdir, "Table5.LoD.part1.G360CDx.55genes.tsv"), sep="\t", header=T, stringsAsFactors = F) %>%
    filter(variant_type == "INDEL") %>%
    mutate(mut_key = variant, lod.maf = lod.probit,   lod.FINAL = LoD) %>%  
    dplyr::select("mut_key", "Input", "lod.probit", "nValid", "level_maf", "lod.maf", "lod.FINAL")
  tab.indel = rbind(tab.indel, tab.indel.agg)%>%
    arrange(Input)
  
  ##  CNVs
  tab.cnv = read.table(file=file.path(resdir, "Table5.LoD.G360CDx.55genes.tsv"), sep="\t", header=T, stringsAsFactors=F) %>% filter(variant_type == "CNV") %>%
    mutate(mut_key = variant, lod.maf = lod.probit,   lod.FINAL = LoD) %>%  
    dplyr::select("mut_key", "Input", "lod.probit", "nValid", "level_maf", "lod.maf", "lod.FINAL")
  
  ##  FUSIONS
  tab.fusion = read.table(file=file.path(resdir, "Table5.LoD.G360CDx.55genes.tsv"), sep="\t", header=T, stringsAsFactors=F) %>% filter(variant_type == "FUSION") %>%
    mutate(mut_key = variant, lod.maf = lod.probit,   lod.FINAL = LoD) %>%  
    dplyr::select("mut_key", "Input", "lod.probit", "nValid", "level_maf", "lod.maf", "lod.FINAL")
  
  tab.all = rbind(tab.snv, tab.indel, tab.cnv, tab.fusion)
  write.table(tab.all, file=file.path(resdir, "Table5.LoD.VARIANTSL.G360CDx.55genes.tsv"), sep="\t", col.names=T, row.names=F, quote = F)
  
  cat("DONE LoD G360CDx method for 55 genes..\n")
