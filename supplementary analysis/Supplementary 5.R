source("src/functions.R")
ref = readRDS("D:/PhD Projects/Rprojects/fastTCGA/data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]
colnames(ref) <- (data.frame(Var2 = colnames(ref)) %>%
                    inner_join(read_csv("data/cellname.csv")))$`Cell Type`
cancer = "HNSC"

gene <- TCGA_transcript_tumor(cancer) 
gene_norm <- TCGA_transcript_normal(cancer) 
cell <-  CIBERSORT(gene, "sample", "symbol", "raw_count", ref, "cibersort", "get") 

hpv_inf <- read_csv("data/HNSC_HPV.csv") %>%
  right_join(read_csv("data/pid.csv")) %>%
  mutate(hpv = ifelse(is.na(`Max aligns Human Papillomavirus`) == T, 0, `Max aligns Human Papillomavirus`)) %>%
  mutate(infection = ifelse(hpv == 0, "uninfected", "infected")) %>%
  dplyr::select(sample, infection, hpv)

gf_cyto_horm <- read_csv("data/consolidated_list_of_gf_hormones_&_cytokines_with_annotation.csv")
gf_cyto_horm_ht <- gene %>%
  filter(symbol %in% gf_cyto_horm$gene) %>%
  group_by(symbol) %>%
  summarise(mean_count = mean(raw_count_scaled)) %>%
  filter(mean_count >= 10) 
gf_cyto_horm <- gf_cyto_horm %>% filter(gene %in% gf_cyto_horm_ht$symbol)


infect_diff <- gene %>%
  inner_join(hpv_inf) %>%
  as_tibble() %>%
  tidybulk(sample, symbol, raw_count) %>%
  test_differential_abundance(~0 + infection, 
                              .contrasts = c("infectioninfected - infectionuninfected"),
                              action = "get") %>%
  rename_with(function(x){gsub("_.*", "",x)}) %>%
  mutate(change = ifelse(logFC >=1.5 & FDR <= 0.05, "Up", "None")) %>%
  mutate(change = ifelse(logFC <=-1.5 & FDR <= 0.05, "Down", change))

up_inf <- infect_diff %>% filter(change == "Up") %>%
  arrange(desc(logFC))
up_inf_gf_cyto_horm <- up_inf %>% inner_join(gf_cyto_horm, by = c("symbol" = "gene")) %>%
  filter(!symbol %in% c("IL19", "IL1R2", "CHGB", "LRP2",
                        "GHRH", "CALCA", "AR")) 

down_inf <-infect_diff %>% filter(change == "Down") %>%
  arrange(desc(logFC))
down_inf_gf_cyto_horm <- down_inf %>% inner_join(gf_cyto_horm, by = c("symbol" = "gene")) %>%
  filter(!symbol %in% c("MTNR1A", "CGB3", "TRH", "GHSR", "UCN3", 
                        "FGF20", "FGF4")) 

hpv_inf <- read_csv("data/HNSC_HPV.csv") %>%
  right_join(read_csv("data/pid.csv")) %>%
  mutate(hpv = ifelse(is.na(`Max aligns Human Papillomavirus`) == T, 0, `Max aligns Human Papillomavirus`)) %>%
  mutate(infection = ifelse(hpv == 0, "HPV-free", "HPV-infected")) %>%
  dplyr::select(sample, infection, hpv)


#### Correlation heatmap
x <- map_dfr(list("All HNSC patients"), function(c){
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