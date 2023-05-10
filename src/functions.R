options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(tidyverse)
library(tidybulk)
library(tidyHeatmap)
library(data.table)
library(survminer)
library(survival)
library(cowplot)
library(scales)
library(ggsci)
library(foreach)
library(reshape2)
library(Hmisc)
library(gridExtra)
library(formattable)
library(psych)
library(clusterProfiler)
library(ggpmisc)
library(ggrepel)
library(Seurat)
library(patchwork)
library(illuminaHumanv4.db)
library(preprocessCore)
library(grid)
library(ggsci)
library("RColorBrewer")
library(ggeasy)
library(singscore)
library(GSEABase)

#url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human"
#download.file(url, destfile = "kegg_2021.gmt")
kegg <- read.gmt("kegg_2021.gmt")

#url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Human"
#download.file(url, destfile = "wp.gmt")
wp <- read.gmt("wp.gmt")

#url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2021"
#download.file(url, destfile = "gobp.gmt")
gobp <- read.gmt("gobp.gmt")

#url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019"
#download.file(url, destfile = "biop.gmt")
biop <- read.gmt("biop.gmt")

#url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016"
#download.file(url, destfile = "react.gmt")
react <- read.gmt("react.gmt")


NKT = c('ReNK', 'IL2NK', 'SPANK', "T Helper", "Naive CD8 T", "GD T", 'CD4 Tcm', 
        'CD4 Tem', 'CD8 Tcm',  'CD8 Tem', "Treg")

#*Aggregate samples and do TMM normalization
TCGA_transcript <- function(cancer){
  as_tibble(
    data.table::setDT(
      map_dfr(as.list(list.files(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/"))), 
              function(i){
                data.table::fread(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/", i)) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}

TCGA_transcript_tumor <- function(cancer){
  as_tibble(
    data.table::setDT(
      map_dfr(as.list(list.files(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/"))), 
              function(i){
                data.table::fread(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/", i)) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/tumor_gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}


TCGA_transcript_normal <- function(cancer){
  as_tibble(
    data.table::setDT(
      map_dfr(as.list(list.files(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/"))), 
              function(i){
                data.table::fread(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/", i)) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/normal_gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}


TCGA_transcript_tumor_raw <- function(cancer){
map_dfr(as.list(list.files(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/"))), 
              function(i){
                data.table::fread(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/Transcript/", i)) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-", cancer, "/tumor_gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)

}

#Combine clinical data
clinical_combine <- function(x, cancer) {
  x %>% 
    inner_join(read.csv(paste0("D:/PhD Projects/Rprojects/fastTCGA/data/Isoforms/TCGA-",cancer, "/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}



#Run CIBERSORT
CIBERSORT <- function(x, c1, c2, c3, ref, meth, act) {
  as_tibble(setDT(x)[, lapply(.SD, sum), by = c(c1, c2), .SDcols = c3]) %>%
    deconvolve_cellularity(
      .sample = !!sym(c1),
      .transcript = !!sym(c2),
      .abundance = !!sym(c3),
      reference = ref,
      method = meth,
      action = act
    )
}

########### Analysis

pw_diff_selection <- function(data, genepair){
  x <- data %>%
    filter(grepl(genepair, comp2)) %>%
    filter(grepl(genepair, comp1)) %>%
    filter(!grepl("HLLH", strata)) %>%
    mutate(gene_combine = genepair, 
           comp1 = substring(comp1, nchar(comp1)-2),
           comp2 = substring(comp2, nchar(comp2)-2)) %>%
    dplyr::select(comp1, comp2,P, adjP, gene_combine) %>%
    na.omit() %>%
    unique()
  rownames(x) <- NULL
  x
}


ggsurvp_cat_item <- function(x, p, nrow, palette, xlab, ylab, title){
  ggsurvplot(
    survfit(
      Surv(total_living_days, vital_status) ~ item,
      data = x
    ),
    data = x,
    facet.by = c("cat"),
    scales = "free_x",
    conf.int = T,
    risk.table = F,
    conf.int.alpha = 0.15,
    pval = p,
    nrow = nrow,
    legend.title = " ",
    short.panel.labs = T,
    palette = palette
  ) + theme_bw() + 
    theme(text=element_text(size=30, family="sans"),
          legend.position = "bottom") +
    labs(x = xlab, y = ylab, title = title) +
    guides(linetype = FALSE)
}


ggsurvp_cat_item_inf <- function(x, p, nrow, palette, xlab, ylab, title){
  ggsurvplot(
    survfit(
      Surv(total_living_days, vital_status) ~ item,
      data = x
    ),
    data = x,
    facet.by = c("cat", "infection"),
    scales = "free_x",
    conf.int = T,
    risk.table = F,
    conf.int.alpha = 0.15,
    pval = p,
    nrow = nrow,
    legend.title = " ",
    short.panel.labs = T,
    palette = palette
  ) + theme_bw() + 
    theme(text=element_text(size=30, family="sans"),
          legend.position = "bottom") +
    labs(x = xlab, y = ylab, title = title) +
    guides(linetype = FALSE)
}

split_median <- function(x, cat, item){
  x %>%
    group_by(!!as.name(cat)) %>%
    mutate(item = factor(Hmisc::cut2(!!as.name(item), g = 2), labels = c(1:nlevels(Hmisc::cut2(!!as.name(item), g = 2))))) %>%
    ungroup() %>%
    mutate(item = gsub("1", "L", item)) %>%
    mutate(item = gsub("2", "H", item))
}



HPVpos_DEG_analysis <- function(inf_gf_cyto_horm, foldername, gf_km_box_cor_width, cyto_km_box_cor_width, horm_km_box_cor_width){
#*correlation heatmap----
    x <- map_dfr(list("All HNSC patients", "HPV-free","HPV-infected"), function(c){
    x <- cell %>%
      inner_join(gene %>%
                   filter(symbol %in% c(unique(inf_gf_cyto_horm$symbol))) %>%
                   dplyr::select(sample, symbol, raw_count_scaled) %>%
                   pivot_wider(names_from = "symbol", values_from = "raw_count_scaled")) %>%
      inner_join(hpv_inf) %>%
      mutate(infection = "All HNSC patients") %>%
      rbind(cell %>%
              inner_join(gene %>%
                           filter(symbol %in% c(unique(inf_gf_cyto_horm$symbol)))%>%
                           dplyr::select(sample, symbol, raw_count_scaled) %>%
                           pivot_wider(names_from = "symbol", values_from = "raw_count_scaled")) %>%
              inner_join(hpv_inf)) %>%
      filter(infection == c) %>%
      dplyr::select(-infection, -hpv)
    
    a <- reshape::melt(cor(x %>% dplyr::select(-sample))[-c(1:24), c(1:24)])
    a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample)))
    a %>% inner_join(as.data.frame(a.p$P[-c(1:24), c(1:24)]) %>% 
                       mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:24), c(1:24)]))) %>% 
                       gather(X2, p.value, -X1)) %>%
      mutate(infection = c)
    
    
  }) %>%
    inner_join(gf_cyto_horm, by = c("X1" = "gene")) %>%
    mutate(infection = fct_relevel(infection, "All HNSC patients", "HPV-free","HPV-infected")) %>%
      mutate(type = gsub("growth_factors_&_receptors", "GF(R)", type)) %>%
      mutate(type = gsub("cytokin/chemokine_&_receptors", "Cyto/Chem(R)", type)) %>%
      mutate(type = gsub("hormones_&_receptors", "Horm(R)", type)) %>%
      mutate(type = fct_relevel(type, "GF(R)", "Cyto/Chem(R)","Horm(R)"))
  
  x$X1 <- factor(x$X1, levels=c(unique(inf_gf_cyto_horm$symbol)), ordered=TRUE)
  x$X2 <- factor(x$X2, levels=colnames(as.data.frame(ref)), ordered=TRUE)
  
  
  p_cor <- as_tibble(x ) %>%
    tidybulk::rename("Cell Type" = "X1", "GFs" = "X2") %>%
    ggplot(aes(x = `Cell Type`, y = GFs, fill = value)) +
    geom_tile(stat = "identity") +
    facet_wrap(infection~type, nrow = 1, scales = "free_x") +
    theme_classic() +
    theme(text=element_text(size=16, family="sans"),
          panel.spacing.x = unit(0,"line"),
          axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
          legend.position = "bottom") +
    labs(x = "", y = "", title = "") + 
    scale_fill_gradient2(low = "#006699", high = "#b30000", mid = "white",
                         midpoint = 0, 
                         name = "Correlation",
                         limits=range(-0.4, 0.4),
                         breaks = seq(-0.4, 0.4, 0.4))
  gt = ggplotGrob(p_cor)
  N <- x %>% group_by(type) %>% 
    summarise(count = length(unique(X1))) %>% 
    `[[`(2)
  
  # Get the column index in the gt layout corresponding to the panels.
  panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
  
  # Replace the default panel widths with relative heights.
  gt$widths[panelI] <- unit(N, "null")
  
  # Add extra width between panels (assuming two panels)
  gt$widths[panelI[1] + 1] = unit(1, "cm")
  ggsave(plot = gt, paste0(foldername, " CorCell_GF_Cyto_Horm.pdf"),
         device = "pdf", height = 10, width = 25)
  
  #####* KM plot/Box plot/Viral-Load Dot plot----
  p_multi <- map2(list("growth_factors_&_receptors", "cytokin/chemokine_&_receptors", "hormones_&_receptors"),
                  list("GF(R)", "Cyto/Chem(R)", "Horm(R)"),
                  function(type0, typetitle){
                    type = typetitle
                    map(as.list(unique((inf_gf_cyto_horm %>%
                                          filter(type == type0))$symbol)),
                        function(tar_gene, type0){
                          x <- gene %>% filter(symbol == tar_gene) %>%
                            inner_join(hpv_inf) %>%
                            mutate(infection = "All HNSC patients") %>%
                            mutate(cat = symbol, item = raw_count_scaled) %>%
                            split_median("cat", "item") %>%
                            mutate(item = fct_relevel(item, "L", "H")) %>%
                            rbind(gene %>% filter(symbol == tar_gene) %>%
                                    inner_join(hpv_inf) %>%
                                    filter(infection == "HPV-infected") %>%
                                    mutate(cat = symbol, item = raw_count_scaled) %>%
                                    split_median("cat", "item") %>%
                                    mutate(item = fct_relevel(item, "L", "H"))) %>%
                            rbind(gene %>% filter(symbol == tar_gene) %>%
                                    inner_join(hpv_inf) %>%
                                    filter(infection == "HPV-free") %>%
                                    mutate(cat = symbol, item = raw_count_scaled) %>%
                                    split_median("cat", "item") %>%
                                    mutate(item = fct_relevel(item, "L", "H"))) %>%
                            clinical_combine(cancer) %>%
                            mutate(cat = infection) %>%
                            mutate(cat = fct_relevel(cat, "All HNSC patients", "HPV-free", "HPV-infected")) %>%
                            as.data.frame()
                          
                          
                          tb <- map_dfr(as.list(c("All HNSC patients", "HPV-free","HPV-infected")), function(c){
                            y <- x %>% filter(cat == c)
                            data.frame(P = c(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,p.adjust.method = "none",
                                                                             data = y)[3])[1,1]),
                                       adj.P = c(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,p.adjust.method = "BH",
                                                                                 data = y)[3])[1,1])) %>% 
                              mutate_if(is.numeric, round, digits = 3) %>%
                              mutate(infection = c)
                          }) %>%
                            right_join(data.frame(infection = c("All HNSC patients", "HPV-free","HPV-infected")))
                          
                          tbs <- lapply(split(tb, tb$infection), function(x) { x["infection"] <- NULL; x })
                          df <- tibble(cat = levels(as.factor(c("All HNSC patients", "HPV-free","HPV-infected"))), 
                                       tbl = tbs)
                          df$cat <- factor(df$cat, levels=c(c("All HNSC patients", "HPV-free","HPV-infected"), ordered=TRUE))
                          
                          p1 <- ggsurvp_cat_item(x, F, 3, "jco", "", "", paste0(tar_gene, " (", type, ")")) + geom_table_npc(data = df, aes(npcx = 0.7, npcy = 0.9, label = tbl),
                                                                                                                             hjust = 0, vjust = 1,
                                                                                                                             table.theme = ttheme_gtlight) 
                          #### Box and scatter plots
                          x <- gene %>%
                            filter(symbol == tar_gene) %>%
                            inner_join(hpv_inf) 
                          
                          #All patients: tar_gene expression and infecteion
                          p1a <- x %>% 
                            mutate(viral = "",
                                   infection = fct_relevel(infection, "HPV-free", "HPV-infected")) %>%
                            ggplot(aes(x=viral, y=log10(raw_count_scaled+1), fill=infection)) + 
                            geom_boxplot() +
                            ylim(0, 5) +
                            stat_compare_means(aes(group = infection), label = "p.signif") +
                            scale_fill_manual(values = c("#4dbbd5ff", "#e64b35ff")) +
                            theme_bw() + 
                            theme(text=element_text(size=16, family="sans"),
                                  legend.position = "bottom",
                                  axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
                            labs(x = "All HNSC patients", y = paste0(tar_gene, " (TMM normalized Log10)"), fill = "HPV infection", title = "")
                          
                          #Only infected patients: tar_gene expression and viral load
                          p1b <- x %>%
                            filter(infection == "HPV-infected") %>%
                            ggplot(aes(y = log10(raw_count_scaled+1), x = log10(hpv +1))) + 
                            geom_point() +
                            stat_smooth(method = "lm", col = "red") +
                            stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.8) +
                            labs(x = "HPV viral load (Log10) \n HPV infected HNSC patients", y = paste0(tar_gene, " (TMM normalized Log10)")) +
                            theme_bw() +
                            ylim(0, 5) +
                            theme(text=element_text(size=16, family="sans"),
                                  legend.position = "bottom") 
                          
                          plot_grid(p1, p1a, p1b, ncol = 1, rel_heights = c(3,1.2,1))
                        })     
                  })
  
  p_gf <- plot_grid(plotlist = p_multi[[1]],   nrow = 1)
  p_cyto <- plot_grid(plotlist = p_multi[[2]], nrow = 1)
  p_horm <- plot_grid(plotlist = p_multi[[3]], nrow = 1)
  
  ggsave(plot = p_gf,
           paste0(foldername, " GF&VL.pdf"), device = "pdf", height = 19, width = gf_km_box_cor_width, limitsize = F)
  ggsave(plot = p_cyto,
         paste0(foldername, " Cyto_Chem&VL.pdf"), device = "pdf", height = 19, width = cyto_km_box_cor_width, limitsize = F)
  ggsave(plot = p_horm,
         paste0(foldername, " Horm&VL.pdf"), device = "pdf", height = 19, width = horm_km_box_cor_width, limitsize = F)
  
}




Figure_1 <- function(cancer, gene, cell, hpv_inf, clinical){
  ######* Panel A: Plot survival infected vs non-infected ######
  x <- hpv_inf %>%
    inner_join(clinical) %>%
    mutate(cat = cancer, item = infection) %>%  
    mutate(item = ifelse(item == "uninfected", "HPV-", "HPV+")) %>%
    mutate(item = fct_relevel(item, "HPV-", "HPV+")) %>%
    as.data.frame()
  
  tb <- data.frame(P = c(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,p.adjust.method = "none",
                                                         data = x)[3])[1,1]),
                   adj.P = c(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,p.adjust.method = "BH",
                                                             data = x)[3])[1,1])) %>% 
    mutate_if(is.numeric, round, digits = 3) %>%
    mutate(cancer = cancer)
  tbs <- lapply(split(tb, tb$cancer), function(x) { x["cancer"] <- NULL; x })
  df <- tibble(cat = cancer, 
               tbl = tbs)
  
  p1a <- ggsurvp_cat_item(x, F, 1, "jco", "", "", "") + 
    geom_table_npc(data = df, aes(npcx = 0.75, npcy = 0.8, label = tbl),
                   hjust = 0, vjust = 1,
                   table.theme = ttheme_gtlight)
  
  ######* Panel B: Box plot Immune cell fractions #####
  x <- cell %>%
    pivot_longer(-sample, names_to = "cat", values_to = "item") %>%
    mutate(cat = ifelse(cat %in% c("Epi", "Endo", "Fibro"), cat, "Immune cell")) %>%
    inner_join(hpv_inf) %>%
    group_by(sample, cat) %>%
    mutate(item = sum(item)) %>%
    ungroup() %>%
    unique() %>%
    filter(cat == "Immune cell") %>%
    mutate(infection = ifelse(infection == "uninfected", "HPV-", "HPV+")) %>%
    mutate(infection = fct_relevel(infection, "HPV-", "HPV+"))
  
  p1b01 <- x %>% ggplot(aes(x=cat, y=logit(item), fill=infection)) + 
    geom_boxplot() +
    stat_compare_means(aes(group = infection), label = "p.signif") +
    scale_fill_manual(values = c("#4dbbd5ff", "#e64b35ff")) +
    theme_bw() + 
    theme(text=element_text(size=20, family="sans"),
          legend.position = "bottom",
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    labs(x = "", y = "Immune Infiltration Logit", fill = "", title = "")
  
  
  p1b02 <- gene %>%
    filter(symbol %in% c("PTPRC")) %>%
    inner_join(hpv_inf) %>%
    mutate(infection = ifelse(infection == "uninfected", "HPV-", "HPV+")) %>%
    mutate(infection = fct_relevel(infection, "HPV-", "HPV+")) %>%
    ggplot(aes(x = symbol, y = z_norm, fill = infection)) +
    geom_boxplot()+
    stat_compare_means(aes(group = infection), label = "p.signif") +
    scale_fill_manual(values = c("#4dbbd5ff", "#e64b35ff")) +
    theme_bw() + 
    theme(text=element_text(size=20, family="sans"),
          legend.position = "bottom",
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    labs(x = "", y = " TMM normalized Log10", fill = "", title = "")
  
  
  p1b <- plot_grid(p1b01, p1b02, nrow = 1) 
  
  
  #######* Panel C: Plot survival 4 groups survival/whole immune infiltration
  x <- cell %>%
    pivot_longer(-sample, names_to = "cat", values_to = "item") %>%
    mutate(cat = ifelse(cat %in% c("Epi", "Endo", "Fibro"), cat, "Immune cell")) %>%
    inner_join(hpv_inf) %>%
    group_by(sample, cat) %>%
    mutate(item = sum(item)) %>%
    ungroup() %>%
    unique() %>%
    filter(cat == "Immune cell") %>%
    split_median("infection", "item") %>%
    unite("item", c("infection", "item"), sep = "/", remove = T) %>%
    inner_join(clinical) %>%
    mutate(item = fct_relevel(item, "uninfected/L", "uninfected/H", "infected/L", "infected/H" ),
           cat = "HPV/Immune Infiltration") %>%
    as.data.frame()
  
  tb <- as.data.frame(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,p.adjust.method = "none",
                                                      data = x)[3]) %>%
                        mutate(comb1 = rownames(.)) %>%
                        pivot_longer(-comb1, names_to = "comb2", values_to = "P") %>%
                        inner_join(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,p.adjust.method = "BH",
                                                                   data = x)[3]) %>%
                                     mutate(comb1 = rownames(.)) %>%
                                     pivot_longer(-comb1, names_to = "comb2", values_to = "adj.P")) %>% 
                        rbind(data.frame(comb1 = "infected/H",
                                         comb2 = "other",
                                         P = as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,p.adjust.method = "none",
                                                                             data = x %>%
                                                                               mutate(item = ifelse(item == "infected/H", "infected/H", "other")))[3])[1,1],
                                         adj.P = as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,p.adjust.method = "none",
                                                                                 data = x %>%
                                                                                   mutate(item = ifelse(item == "infected/H", "infected/H", "other")))[3])[1,1])) %>%
                        mutate_all(funs(str_replace(., "p.value.", " "))) %>% 
                        mutate(comb1 = str_replace(comb1, "\\.", "/"),
                               comb2 = str_replace(comb2, "\\.", "/"),
                               adj.P = as.numeric(adj.P),
                               P = as.numeric(P)) %>%
                        mutate_if(is.numeric, round, digits = 3) %>%
                        na.omit() %>%
                        mutate(cat = "HPV/Immune Infiltration"))[-c(3,4),]
  tb <- tb %>% 
    mutate_all(funs(str_replace(., "uninfected", "HPV-"))) %>%
    mutate_all(funs(str_replace(., "infected", "HPV+")))
  
  tbs <- lapply(split(tb, tb$cat), function(x) { x["cat"] <- NULL; x })
  df <- tibble(cat = "HPV/Immune Infiltration", 
               tbl = tbs)
  x <-  x %>%
    mutate(item = gsub("uninfected", "HPV-", item)) %>%
    mutate(item = gsub("infected", "HPV+", item)) %>%
    mutate(item = fct_relevel(item, "HPV-/L", "HPV-/H", "HPV+/L", "HPV+/H" ))
  
  
  
  p1c <- ggsurvp_cat_item(x, F, 1, "jco", "", "Survival probability", "") + 
    geom_table_npc(data = df, aes(npcx = 0.5, npcy = 0.8, label = tbl),
                   hjust = 0, vjust = 1,
                   table.theme = ttheme_gtlight) 
  
  
  plot_grid(p1a, p1b, p1c, nrow = 1, rel_widths = c(1,1.5,1), align = "hv")
  
}



Figure_2 <- function(cancer, gene, cell, hpv_inf, clinical){
  ######* Panel A: Box plot of all cells HPV- vs HPV+ ######
  x <- cell %>%
    pivot_longer(-sample, names_to = "cat", values_to = "item") %>%
    inner_join(hpv_inf) %>%
    unique() %>%
    mutate(infection = ifelse(infection == "uninfected", "HPV-", "HPV+")) %>%
    mutate(type = ifelse(cat %in% c("Macro M1","Macro M2", "Mono","Neutro", "Eosino", "Mast","iDC", "mDC"), "Myeloid", "other")) %>%
    mutate(type = ifelse(cat %in% c("Endo", "Epi", "Fibro"), "Stromal", type)) %>%
    mutate(type = ifelse(!type %in% c("Myeloid", "Stromal"), "Lymphoid", type)) %>%
    mutate(cat = fct_relevel(cat, cellorder),
           infection = fct_relevel(infection, "HPV-", "HPV+"),
           type = fct_relevel(type, "Stromal", "Myeloid", "Lymphoid")) 
  
  
  p2a <- x %>% ggplot(aes(x=cat, y=logit(item), fill=infection)) + 
    geom_boxplot() +
    facet_wrap(.~type, scales = "free") +
    stat_compare_means(aes(group = infection), label = "p.signif") +
    scale_fill_manual(values = c("#4dbbd5ff", "#e64b35ff")) +
    theme_bw() + 
    ylim(-12.5, 0) +
    theme(text=element_text(size=20, family="sans"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
    labs(x = "", y = "Immune Infiltration Logit", fill = "", title = "", tag = "A")
  
  gt = ggplotGrob(p2a)
  N <- x %>% group_by(type) %>% 
    summarise(count = length(unique(cat))) %>% 
    `[[`(2)
  
  # Get the column index in the gt layout corresponding to the panels.
  panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
  
  # Replace the default panel widths with relative heights.
  gt$widths[panelI] <- unit(N, "null")
  
  # Add extra width between panels (assuming two panels)
  gt$widths[panelI[1] + 1] = unit(1, "cm")
  p2a <- gt
  
  ######* Panel B: KM plot 4 groups HPV/Mem B&SPANK ######
  x <- cell %>%
    pivot_longer(-sample, names_to = "cat", values_to = "item") %>%
    filter(cat %in% c("Memory B", "SPANK")) %>%
    inner_join(hpv_inf) %>%
    group_by(infection, cat) %>%
    mutate(item = factor(Hmisc::cut2(item, g = 2), labels = c("L", "H"))) %>%
    ungroup() %>%
    unite("item", c("infection", "item"), sep = "/", remove = T) %>%
    inner_join(clinical)
  
  y <- x %>% mutate(item = ifelse(item == "infected/H", "infected/H", "other"))
  
  tb <- as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ cat + item,
                                        p.adjust.method = "none",
                                        data = x)[3]) %>%
    mutate(cat = rownames(.)) %>%
    pivot_longer(-cat, names_to = "comb2", values_to = "P") %>%
    na.omit() %>%
    separate("cat", c("cat", "comb1"), sep = ",") %>%
    mutate(cat = gsub(".*=", "", cat),
           comb1 = gsub(".*=", "", comb1),
           comb2 = gsub(".*cat.", "", comb2)) %>%
    separate("comb2", c("cat2", "comb2"), sep = "..item.") %>%
    mutate(cat2 = gsub("\\.", " ", cat2)) %>%
    filter(cat == cat2) %>%
    mutate(comb2 = gsub("\\.\\.", "", comb2)) %>%
    mutate(comb2 = gsub("\\.", "/", comb2)) %>%
    unite("test", c("comb1", "comb2"), sep = "/", remove = F) %>%
    filter(!test %in% c("uninfected/H/infected/L", "uninfected/L/infected/H")) %>%
    inner_join(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ cat + item,
                                               p.adjust.method = "BH",
                                               data = x)[3]) %>%
                 mutate(cat = rownames(.)) %>%
                 pivot_longer(-cat, names_to = "comb2", values_to = "adj.P") %>%
                 na.omit() %>%
                 separate("cat", c("cat", "comb1"), sep = ",") %>%
                 mutate(cat = gsub(".*=", "", cat),
                        comb1 = gsub(".*=", "", comb1),
                        comb2 = gsub(".*cat.", "", comb2)) %>%
                 separate("comb2", c("cat2", "comb2"), sep = "..item.") %>%
                 mutate(cat2 = gsub("\\.", " ", cat2)) %>%
                 filter(cat == cat2) %>%
                 mutate(comb2 = gsub("\\.\\.", "", comb2)) %>%
                 mutate(comb2 = gsub("\\.", "/", comb2)) %>%
                 unite("test", c("comb1", "comb2"), sep = "/", remove = F) %>%
                 filter(!test %in% c("uninfected/H/infected/L", "uninfected/L/infected/H")))  %>%
    dplyr::select(-test, -cat2) %>%
    mutate_if(is.numeric, round, digits = 3) %>%
    rbind(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ cat + item,
                                          p.adjust.method = "none",
                                          data = y)[3]) %>%
            mutate(cat = rownames(.)) %>%
            pivot_longer(-cat, names_to = "comb2", values_to = "P") %>%
            na.omit() %>%
            separate("cat", c("cat", "comb1"), sep = ",") %>%
            mutate(cat = gsub(".*=", "", cat),
                   comb1 = gsub(".*=", "", comb1),
                   comb2 = gsub(".*cat.", "", comb2)) %>%
            separate("comb2", c("cat2", "comb2"), sep = "..item.") %>%
            mutate(cat2 = gsub("\\.", " ", cat2)) %>%
            filter(cat == cat2) %>%
            mutate(comb2 = gsub("\\.\\.", "", comb2)) %>%
            mutate(comb2 = gsub("\\.", "/", comb2)) %>%
            unite("test", c("comb1", "comb2"), sep = "/", remove = F) %>%
            inner_join(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ cat + item,
                                                       p.adjust.method = "BH",
                                                       data = y)[3]) %>%
                         mutate(cat = rownames(.)) %>%
                         pivot_longer(-cat, names_to = "comb2", values_to = "adj.P") %>%
                         na.omit() %>%
                         separate("cat", c("cat", "comb1"), sep = ",") %>%
                         mutate(cat = gsub(".*=", "", cat),
                                comb1 = gsub(".*=", "", comb1),
                                comb2 = gsub(".*cat.", "", comb2)) %>%
                         separate("comb2", c("cat2", "comb2"), sep = "..item.") %>%
                         mutate(cat2 = gsub("\\.", " ", cat2)) %>%
                         filter(cat == cat2) %>%
                         mutate(comb2 = gsub("\\.\\.", "", comb2)) %>%
                         mutate(comb2 = gsub("\\.", "/", comb2)) %>%
                         unite("test", c("comb1", "comb2"), sep = "/", remove = F) ) %>%
            dplyr::select(-test, -cat2) %>%
            mutate_if(is.numeric, round, digits = 3)) %>%
    mutate(`Strata 1` = comb2, `Strata 2` = comb1) %>% 
    mutate_all(funs(str_replace(., "uninfected", "HPV-"))) %>%
    mutate_all(funs(str_replace(., "infected", "HPV+"))) %>%
    dplyr::select(cat, `Strata 1`, `Strata 2`, P, adj.P)
  
  tbs <- lapply(split(tb, tb$cat), function(x) { x["cat"] <- NULL; x })
  df <- tibble(cat = levels(as.factor(c("Memory B", "SPANK"))), 
               tbl = tbs)
  x <-  x %>%
    mutate(item = gsub("uninfected", "HPV-", item)) %>%
    mutate(item = gsub("infected", "HPV+", item)) %>%
    mutate(item = fct_relevel(item, "HPV-/L", "HPV-/H", "HPV+/L", "HPV+/H" )) %>%
    as.data.frame()
  
  p2b <- ggsurvp_cat_item(x, F, 1, "jco", "", "Survival probability", "") + 
    geom_table_npc(data = df, aes(npcx = 0.5, npcy = 0.8, label = tbl),
                   hjust = 0, vjust = 1,
                   table.theme = ttheme_gtlight) +
    labs(tag = "B")
  
  plot_grid(p2a, p2b, ncol = 1)
  
}

strata_p <- function(x, method0, p_name0){

  as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ cat + item,
                                  p.adjust.method = method0,
                                  data = x)[3]) %>%
    mutate(cat = rownames(.)) %>%
    pivot_longer(-cat, names_to = "comb2", values_to = p_name0) %>%
    na.omit() %>%
    separate("cat", c("cat", "comb1"), sep = ",") %>%
    mutate(cat = gsub(".*=", "", cat),
           comb1 = gsub(".*=", "", comb1),
           comb2 = gsub(".*cat.", "", comb2)) %>%
    mutate(comb2 = gsub("\\.", "/", comb2),
           cat1 = gsub(" ", "/", cat)) %>%
    separate("comb2", c("cat2", "comb2"), sep = "//item/") %>%   ### keep cat2 in complicated analysis for double check
    unite("test", c("comb1", "comb2"), sep = "/", remove = F) %>%
    filter(!test %in% c("H/L/L/H", "L/H/H/L", "L/L/H/H", "H/H/L/L")) %>%
    filter(cat1 == cat2) %>% 
    dplyr::select(cat, comb1, comb2, !!as.name(p_name0))
}

paste_char2 <- function(char1, char2, sep0){
  c(sapply(char1, function(x) {paste(x,char2, sep = sep0)}))
}


