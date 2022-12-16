source("src/functions.R")
ref <- readRDS("data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]

#distinguish hpv_status and state
x <- all_factor %>%
  dplyr::select(sample, hpv_status_ish, hpv_status_p16, state, PVLog) %>%
  mutate(hpv_status_TCGA = ifelse(hpv_status_ish == "Positive"|hpv_status_p16 == "Positive", "Positive", NA)) %>%
  mutate(hpv_status_TCGA = ifelse(hpv_status_ish == "Negative"|hpv_status_p16 == "Negative", "Negative", hpv_status_TCGA)) %>%
  mutate(diff = ifelse(state == "infected"&hpv_status_TCGA != "Positive", "different", NA)) %>%
  mutate(diff = ifelse(state == "uninfected"&hpv_status_TCGA != "Negative", "different", diff))

x %>% filter(diff == "different") %>% write.csv("output/hpv_status_TCGA&paper.csv", row.names = F)
################### A table of PC correlated factors
x <- factor_numerical %>%
  dplyr::select(sample, gender, hpv_status_ish, hpv_status_p16, state, PVLog, epithelial, 
                neutrophil, dendritic_myeloid_mature, perineural_invasion)

a <- Hmisc::rcorr(as.matrix(x %>% dplyr::select(-sample)))
a <- as.data.frame(a$r) %>% 
  mutate(Var1 = rownames(a$r)) %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "corr") %>%
  ggplot(aes(x = Var1, y = Var2, fill = corr)) +
  geom_tile(stat = "identity") +
  theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "") + 
  scale_fill_gradient2(low = "#00284d", high = "#b30000", mid = "white",
                       midpoint = 0, 
                       name = "Correlation")

#+ scale_fill_distiller(palette = "RdBu", limits = c(-1, 1))

ggsave(plot = a, "output/PC_correlated_factor.pdf", device = "pdf", height = 10, width = 10)

a <- Hmisc::rcorr(as.matrix(x %>% dplyr::select(-sample)))
a <- as_tibble(as.data.frame(a$r) %>% 
                 mutate(Var1 = rownames(a$r)) %>%
                 pivot_longer(-Var1, names_to = "Var2", values_to = "corr")) %>%
  tidyHeatmap::heatmap(
    Var1,
    Var2,
    corr,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(1,-1, length.out = 6),
      
      RColorBrewer::brewer.pal(6, "RdBu")
    )
    
  ) #%>% tidyHeatmap::save_pdf("output/linear.pdf", height = 30, width = 30)

################### A table including the tissue_source_site of patients and their infection status as well as viral load (PVLog). 
x <- HNSC_PV %>% 
  inner_join(read_csv("data/clinical_hnsc_numerical.csv")) %>%
  inner_join(read_csv("data/clinical_hnsc_binary.csv")) %>%
  inner_join(read_csv("data/clinical_hnsc_categorical.csv")) %>%
  dplyr::select(sample, state, tissue_source_site, PVLog) %>%
  unique()

x %>% group_by(tissue_source_site) %>% summarise(n = n()) %>% view()

x %>%  
  group_by(state, tissue_source_site) %>%
  summarise(n = n()) %>%
  mutate(perc = n/sum(n)) %>%
  mutate(perc_lab = paste0(round(perc, 3)*100, "%")) %>%
  ggplot(aes(x = state, y = perc, fill = tissue_source_site, label = perc_lab)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(size = 1, position = position_stack(vjust = 0.5)) +
  theme_bw()


x %>%  
  ggplot(aes(x = state, y = PVLog, fill = tissue_source_site)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_()


#################BOXplot
x %>% 
  ggplot(aes(x = tissue_source_site, y = PVLog, fill = tissue_source_site)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.6, aes( fill = tissue_source_site)) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Tissue Source Site (Laboratories)", y = "Log Transformed Viral Load", title = "TCGA-HNSC", fill = "TSS") + 
  theme(text=element_text(size=16, family="sans")) 

ggsave("output/TSS_boxplot.pdf", height = 10, width = 21, device = "pdf")

##################### All factor correlation
factor_numerical <- HNSC_PV %>%
  dplyr::select(sample, PVLog, state, celltype, logit_fraction) %>%
  pivot_wider(names_from = "celltype", values_from = "logit_fraction") %>% 
  inner_join(read_csv("data/clinical_hnsc_numerical.csv")) %>%
  inner_join(read_csv("data/clinical_hnsc_binary.csv")) %>%
  mutate(state = ifelse(state == "uninfected", 0, 1))

factor_categorical <- HNSC_PV %>%
  group_by(sample, state) %>%
  summarise() %>%
  ungroup %>%
  inner_join(read_csv("data/clinical_hnsc_categorical.csv")) %>%
  dplyr::select(-state)

x <- factor_categorical %>%
  gather(cat, item, -sample)

factor_cat <-   foreach(col = colnames(factor_categorical)[-c(1, 11, 12, 22)], .combine = inner_join) %do% {
  number <- x %>% filter(cat == col) %>%
    group_by(item) %>%
    summarise() %>%
    mutate(num = ifelse(is.na(item) == T, NA, row_number()))
  
  colnames(number)[1] <- col
  
  y <- factor_categorical %>% inner_join(number) %>%
    dplyr::select(sample, num)
  colnames(y)[2] <- col
  
  y
} %>%
  inner_join(factor_categorical[, c(1,11,12,22)])


############## PCA analysis and factors
hnsc_pca <- HNSC %>%
  as_tibble() %>%
  reduce_dimensions(.element = sample, .feature = symbol, .abundance = raw_count_scaled, method="PCA", .dims = 10, action = "get") %>%
  dplyr::select(-TMM, -multiplier)

HNSC %>%
  as_tibble() %>%
  reduce_dimensions(method="PCA", .dims = 20) 

p <- hnsc_pca %>% 
  inner_join(all_factor %>%
               dplyr::select(sample, tumor_grade)) %>%
  GGally::ggpairs(columns = 2:11, ggplot2::aes(colour=tumor_grade))
ggsave(plot = p, "output/tumor_grade.pdf", device = "pdf", height = 30, width = 30)

p <- hnsc_pca %>% 
  inner_join(all_factor %>%
               dplyr::select(sample, PVLog)) %>%
  mutate(PVLog = factor(Hmisc::cut2(PVLog, g = 2),labels = c(1:nlevels(Hmisc::cut2(PVLog, g = 2))))) %>%
  GGally::ggpairs(columns = 2:11, ggplot2::aes(colour=PVLog))
ggsave(plot = p, "output/PCA/PVLog.pdf", device = "pdf", height = 30, width = 30)

read_csv("data/PCs.csv") %>%
  mutate(x = "PC") %>%
  ggplot(aes(x = x, y = Variance, fill = PC, label = paste0("PC",PC, ":", round(Variance, 3)*100, "%"))) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Spectral") +
  coord_polar(theta="y") +
  theme_bw()  + 
  theme(text=element_text(size=16, family="sans"))
ggsave("output/PC_variance.pdf", device = "pdf", height = 15, width = 17)


linearpca <- Hmisc::rcorr(as.matrix(hnsc_pca %>% 
                                      inner_join(factor_numerical) %>% dplyr::select(-sample)))

x <- as.data.frame(linearpca$r[-c(1:10),1:10]) %>% 
  mutate(Var1 = rownames(linearpca$r[-c(1:10),1:10])) %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "corr")

x$Var2 <- factor(x$Var2, levels=c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), ordered=TRUE)

linearpca <- x %>% ggplot(aes(x = Var1, y = Var2, fill = corr)) +
  geom_tile(stat = "identity") +
  scale_fill_distiller(palette = "RdBu", name = "Correlation") +
  theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "") 

ggsave(plot = linearpca, "output/linearpca.pdf", device = "pdf", height = 10, width = 15)


x <- foreach(i = c(2:11), .combine = bind_rows) %do% {
  foreach(j = c(2:27), .combine = bind_rows) %do% {
    data.frame(Var1 = c(colnames(hnsc_pca)[i]), 
               Var2= c(colnames(factor_categorical)[j]), 
               aov_cohens_f = cohens_f(aov(as.numeric(as.data.frame(hnsc_pca[,i])[,1])~as.character(as.data.frame(factor_categorical[,j])[,1])))[1,2])
  }
}  

x$Var1 <- factor(x$Var1, levels=c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), ordered=TRUE)

anovapca <- x %>%
  ggplot(aes(x = Var2, y = Var1, fill = aov_cohens_f)) +
  geom_tile(stat = "identity") +
  scale_fill_distiller(palette = "RdBu", name = "Cohen's f", limits = c(min = 0, max = 0.8)) +
  theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "") 

ggsave(plot = anovapca, "output/anovapca.pdf", device = "pdf", height = 10, width = 10)


x <- HNSC %>% inner_join(hnsc_pca) %>%
  dplyr::select(sample, symbol, raw_count_scaled, PC7) %>%
  pivot_wider(names_from = "symbol", values_from = "raw_count_scaled")



PC7genes <- map_dfr(as.list(3:34448), function(i){
  data.frame(symbol = c(colnames(x[,i])), PC7cor = c(cor(x[,2], x[,i])[1,1]))
}) %>%
  mutate(cor_value = ifelse(PC7cor >= 0, PC7cor, -PC7cor)) %>%
  arrange(cor_value, desc(cor_value)) %>%
  na.omit() %>%
  as_tibble() %>%
  describe_transcript(.transcript = symbol)

PC7genes %>% write.csv("PC7relatedgenes.csv", row.names = F)

####################   ANOVA analysis
anova <- foreach(i = c(2:8, 10:53), .combine = bind_rows) %do% {
  foreach(j = c(2:10, 12:27), .combine = bind_rows) %do% {
    data.frame(Var1 = c(colnames(factor_numerical)[i]), 
               Var2= c(colnames(factor_categorical)[j]), 
               aov_cohens_f = cohens_f(aov(as.numeric(as.data.frame(factor_numerical[,i])[,1])~as.character(as.data.frame(factor_categorical[,j])[,1])))[1,2])
  }
}  %>%
  ggplot(aes(x = Var1, y = Var2, fill = aov_cohens_f)) +
  geom_tile(stat = "identity") +
  scale_fill_distiller(palette = "RdBu", name = "Cohen's f", limits = c(min = 0, max = 2)) +
  theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "") 

ggsave(plot = anova, "output/anova.pdf", device = "pdf", height = 15, width = 20)


p_anova <- as_tibble(anova) %>%
  tidyHeatmap::heatmap(
    Var1,
    Var2,
    aov_cohens_f,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE
    
  ) # %>% tidyHeatmap::save_pdf("output/anova.pdf", height = 30, width = 20)


##############################  Linear correlation
linear <- Hmisc::rcorr(as.matrix(factor_numerical %>% dplyr::select(-sample)))
linear <- as.data.frame(linear$r) %>% 
  mutate(Var1 = rownames(linear$r)) %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "corr") %>%
  ggplot(aes(x = Var1, y = Var2, fill = corr)) +
  geom_tile(stat = "identity") +
  scale_fill_distiller(palette = "RdBu", name = "Correlation") +
  theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "") 

ggsave(plot = linear, "output/linear.pdf", device = "pdf", height = 20, width = 20)


p_linear <- as_tibble(linear) %>%
  tidyHeatmap::heatmap(
    Var1,
    Var2,
    corr,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(2,-2, length.out = 11),
      
      RColorBrewer::brewer.pal(11, "RdBu")
    )
    
  ) %>% tidyHeatmap::save_pdf("output/linear.pdf", height = 30, width = 30)


##################### Memory B cells/SPANK and GF/GFR genes
x <- HNSC %>%
  inner_join(read_csv("data/Human growth factor gene list.csv")) %>%
  dplyr::select(sample, symbol, raw_count_scaled) %>%
  unique() %>%
  pivot_wider(names_from = "symbol", values_from = "raw_count_scaled") %>%
  full_join(HNSC_cell %>% dplyr::select(sample, b_memory, nk_primed_IL2_PDGFD)) 

a <- Hmisc::rcorr(as.matrix(x %>% dplyr::select(-sample)))

as.data.frame(a$r[-c(1:222), c(1:222)]) %>% 
  mutate(Var1 = rownames(a$r[-c(1:222), c(1:222)])) %>%
  mutate(Var1 = ifelse(Var1 == "b_memory", "Memory B", "SPANK")) %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "corr") %>%
  ggplot(aes(x = Var1, y = Var2, fill = corr)) +
  geom_tile(stat = "identity") +
  scale_fill_distiller(palette = "RdBu", name = "Correlation") +
  theme_bw() +
  theme(text=element_text(size=16, family="sans"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "") 

ggsave("output/b_spank_GFs.pdf", device = "pdf", height = 40, width = 5)

##################### Heatmap: Cells and stage
x <- as.data.frame(HNSC_cell %>%
                     dplyr::select(sample, b_memory, contains("t_"), nk_primed_IL2_PDGFD, macrophage_M1) %>%
                     dplyr::select(-mast_cell) %>%
                     inner_join(HNSC_PV %>% 
                                  filter(state == "infected") %>%
                                  dplyr::select(sample, PVLog) %>%
                                  inner_join(read_csv("data/clinical_hnsc.csv") %>%
                                               dplyr::select(sample, stage) %>%
                                               mutate(stage = gsub("Stage IV.", "Stage IV", stage)) %>%
                                               filter(grepl("Stage", stage)) %>%
                                               mutate(stage = ifelse(stage == "Stage I"|stage == "Stage II", "Develop(I/II)", "Advanced(III/Iv)"))) %>%
                                  filter(stage == "Develop(I/II)") %>%
                                  unique() %>%
                                  spread(stage, PVLog))
) 

a <- reshape::melt(cor(x %>% dplyr::select(-sample))[-c(1:11), c(1:11)])
a <- a %>% 
  tidybulk::rename("Correlation" = "value") %>%
  mutate(Var2 = rownames(as.data.frame(a))) %>% 
  inner_join(read_csv("data/cellname.csv")) %>%
  mutate(Stage = "Developing (I/II)")

as_tibble(a ) %>%
  tidyHeatmap::heatmap(
    Stage,
    `Cell Type`,
    Correlation,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(2,-2, length.out = 11),
      
      RColorBrewer::brewer.pal(11, "RdBu")
      
    )
    
  ) %>%
  tidyHeatmap::save_pdf("output/develop_stage_HNSC_cell_corr.pdf", height = 2, width = 10)


x <- as.data.frame(HNSC_cell %>%
                     dplyr::select(sample, b_memory, contains("t_"), nk_primed_IL2_PDGFD, macrophage_M1) %>%
                     dplyr::select(-mast_cell) %>%
                     inner_join(HNSC_PV %>% 
                                  filter(state == "infected") %>%
                                  dplyr::select(sample, PVLog) %>%
                                  inner_join(read_csv("data/clinical_hnsc.csv") %>%
                                               dplyr::select(sample, stage) %>%
                                               mutate(stage = gsub("Stage IV.", "Stage IV", stage)) %>%
                                               filter(grepl("Stage", stage)) %>%
                                               mutate(stage = ifelse(stage == "Stage I"|stage == "Stage II", "Develop(I/II)", "Advanced(III/Iv)"))) %>%
                                  filter(stage == "Advanced(III/Iv)") %>%
                                  unique() %>%
                                  spread(stage, PVLog))
) 

a <- reshape::melt(cor(x %>% dplyr::select(-sample))[-c(1:11), c(1:11)])
a <- a %>% 
  tidybulk::rename("Correlation" = "value") %>%
  mutate(Var2 = rownames(as.data.frame(a))) %>% 
  inner_join(read_csv("data/cellname.csv")) %>%
  mutate(Stage = "Advanced(III/Iv)")

as_tibble(a ) %>%
  tidyHeatmap::heatmap(
    Stage,
    `Cell Type`,
    Correlation,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(2,-2, length.out = 11),
      
      RColorBrewer::brewer.pal(11, "RdBu")
      
    )
    
  ) %>%
  tidyHeatmap::save_pdf("output/advanced_stage_HNSC_cell_corr.pdf", height = 2, width = 10)


##################### Heatmap: Genes and stage
x <- as.data.frame(foreach(i = c("IL17REL", "FGF3", "IL1R2", "NGF", "FGF19", "CSF2"), .combine = bind_rows) %do% {
  HNSC %>% filter(symbol == i) %>%
    dplyr::select(sample, symbol, raw_count_scaled) 
} %>% spread(symbol, raw_count_scaled) %>%
  inner_join(HNSC_PV %>% 
               filter(state == "infected") %>%
               dplyr::select(sample, PVLog) %>%
               inner_join(read_csv("data/clinical_hnsc.csv") %>%
                            dplyr::select(sample, stage) %>%
                            mutate(stage = gsub("Stage IV.", "Stage IV", stage)) %>%
                            filter(grepl("Stage", stage)) %>%
                            mutate(stage = ifelse(stage == "Stage I"|stage == "Stage II", "Develop(I/II)", "Advanced(III/Iv)"))) %>%
               filter(stage == "Develop(I/II)") %>%
               unique() %>%
               spread(stage, PVLog)))


a <- reshape::melt(cor(x %>% dplyr::select(-sample))[-c(1:6), c(1:6)])
a <- a %>% 
  tidybulk::rename("Correlation" = "value") %>%
  mutate(Symbol = rownames(as.data.frame(a))) %>% 
  mutate(Stage = "Developing (I/II)")

as_tibble(a ) %>%
  tidyHeatmap::heatmap(
    Stage,
    Symbol,
    Correlation,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(2,-2, length.out = 11),
      
      RColorBrewer::brewer.pal(11, "RdBu")
      
    )
    
  ) %>%
  tidyHeatmap::save_pdf("output/develop_stage_HNSC_gene_corr.pdf", height = 2, width = 10)


x <- as.data.frame(foreach(i = c("IL17REL", "FGF3", "IL1R2", "NGF", "FGF19", "CSF2"), .combine = bind_rows) %do% {
  HNSC %>% filter(symbol == i) %>%
    dplyr::select(sample, symbol, raw_count_scaled) 
} %>% spread(symbol, raw_count_scaled) %>%
  inner_join(HNSC_PV %>% 
               filter(state == "infected") %>%
               dplyr::select(sample, PVLog) %>%
               inner_join(read_csv("data/clinical_hnsc.csv") %>%
                            dplyr::select(sample, stage) %>%
                            mutate(stage = gsub("Stage IV.", "Stage IV", stage)) %>%
                            filter(grepl("Stage", stage)) %>%
                            mutate(stage = ifelse(stage == "Stage I"|stage == "Stage II", "Develop(I/II)", "Advanced(III/Iv)"))) %>%
               filter(stage == "Advanced(III/Iv)") %>%
               unique() %>%
               spread(stage, PVLog)))

a <- reshape::melt(cor(x %>% dplyr::select(-sample))[-c(1:6), c(1:6)])
a <- a %>% 
  tidybulk::rename("Correlation" = "value") %>%
  mutate(Symbol = rownames(as.data.frame(a))) %>% 
  mutate(Stage = "Advanced(III/Iv)")

as_tibble(a ) %>%
  tidyHeatmap::heatmap(
    Stage,
    Symbol,
    Correlation,
    cluster_rows = FALSE,
    
    cluster_columns = FALSE,
    
    palette_value = circlize::colorRamp2(
      
      seq(2,-2, length.out = 11),
      
      RColorBrewer::brewer.pal(11, "RdBu")
      
    )
    
  ) %>%
  tidyHeatmap::save_pdf("output/advanced_stage_HNSC_gene_corr.pdf", height = 2, width = 10)

