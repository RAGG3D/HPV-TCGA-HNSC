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

HPVpos_DEG_analysis(up_inf_gf_cyto_horm, "paper figures/output/Figure 4", 16, 16, 16)






#########* Combination of MemB and SPANK in infected patients #####
x <- cell %>%
  dplyr::select(sample, `Memory B`, SPANK) %>%
  mutate(`Memory B` = factor(Hmisc::cut2(`Memory B`, g = 2), labels = c("L", "H"))) %>%
  mutate(SPANK = factor(Hmisc::cut2(SPANK, g = 2), labels = c("L", "H"))) %>%
  clinical_combine(cancer) %>%
  unite("item", c("Memory B", "SPANK"), sep = "/", remove = T) %>%
  mutate(cat = "MemB/SPANK") %>%
  inner_join(hpv_inf) %>%
  filter(infection == "infected") %>%
  as.data.frame() %>%
  mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H")) 

tb <- map_dfr(as.list(c("MemB/SPANK")), function(t){
  y <- x %>% filter(cat == t)
  as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                  p.adjust.method = "none",
                                  data = y)[3]) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname, names_to = "comb2", values_to = "P") %>%
    na.omit() %>%
    mutate(comb2 = gsub("\\.", "/", comb2)) %>%
    mutate(comb2 = gsub(".*value/", "", comb2)) %>%
    mutate_if(is.numeric, round, digits = 3) %>%
    tidybulk::rename("comb1" = "rowname") %>%
    mutate(cat = t,
           test = paste(comb1, comb2, sep = "/")) %>%
    filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
    dplyr::select(cat, comb1, comb2, P) %>%
    inner_join(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                               p.adjust.method = "BH",
                                               data = y)[3]) %>%
                 rownames_to_column() %>%
                 pivot_longer(-rowname, names_to = "comb2", values_to = "P.adj") %>%
                 na.omit() %>%
                 mutate(comb2 = gsub("\\.", "/", comb2)) %>%
                 mutate(comb2 = gsub(".*value/", "", comb2)) %>%
                 mutate_if(is.numeric, round, digits = 3) %>%
                 tidybulk::rename("comb1" = "rowname") %>%
                 mutate(cat = t,
                        test = paste(comb1, comb2, sep = "/")) %>%
                 filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
                 dplyr::select(cat, comb1, comb2, P.adj))
})


tbs <- lapply(split(tb, tb$cat), function(x) { x["cat"] <- NULL; x })
df <- tibble(cat = levels(as.factor(c("MemB/SPANK"))), 
             tbl = tbs)

p7 <- ggsurvp_cat_item(x, F, 1, "jco", "", "", "HNSC HPV infected Tumor") + geom_table_npc(data = df, aes(npcx = 0.58, npcy = 1, label = tbl),
                                                                                           hjust = 0, vjust = 1,
                                                                                           table.theme = ttheme_gtlight) 
ggsave(plot = p7,
       "paper figures/output/Figure 7.pdf", device = "pdf",
       height = 6, width = 6)

#########* Combination of IL17R/MemB SPANK in infected patients ######
x <- map_dfr(as.list(c("IL17RB", "IL17REL")), function(il17r){
  cell %>%
    dplyr::select(sample, `Memory B`, SPANK) %>%
    mutate(`Memory B` = factor(Hmisc::cut2(`Memory B`, g = 2), labels = c("L", "H"))) %>%
    mutate(SPANK = factor(Hmisc::cut2(SPANK, g = 2), labels = c("L", "H"))) %>%
    pivot_longer(-sample, names_to = "celltype", values_to = "fraction") %>%
    inner_join(gene %>%
                 filter(symbol %in% c(il17r)) %>%
                 split_median("symbol", "raw_count_scaled")) %>%
    clinical_combine(cancer) %>%
    unite("item", c("fraction", "item"), sep = "/", remove = T) %>%
    unite("cat", c("celltype", "symbol"), sep = "/", remove = T) %>%
    inner_join(hpv_inf) %>%
    filter(infection == "infected") %>%
    as.data.frame() %>%
    mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H"))
})


tb <- map_dfr(as.list(unique(x$cat)), function(t){
  y <- x %>% filter(cat == t)
  as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                  p.adjust.method = "none",
                                  data = y)[3]) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname, names_to = "comb2", values_to = "P") %>%
    na.omit() %>%
    mutate(comb2 = gsub("\\.", "/", comb2)) %>%
    mutate(comb2 = gsub(".*value/", "", comb2)) %>%
    mutate_if(is.numeric, round, digits = 3) %>%
    tidybulk::rename("comb1" = "rowname") %>%
    mutate(cat = t,
           test = paste(comb1, comb2, sep = "/")) %>%
    filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
    dplyr::select(cat, comb1, comb2, P) %>%
    inner_join(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                               p.adjust.method = "BH",
                                               data = y)[3]) %>%
                 rownames_to_column() %>%
                 pivot_longer(-rowname, names_to = "comb2", values_to = "P.adj") %>%
                 na.omit() %>%
                 mutate(comb2 = gsub("\\.", "/", comb2)) %>%
                 mutate(comb2 = gsub(".*value/", "", comb2)) %>%
                 mutate_if(is.numeric, round, digits = 3) %>%
                 tidybulk::rename("comb1" = "rowname") %>%
                 mutate(cat = t,
                        test = paste(comb1, comb2, sep = "/")) %>%
                 filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
                 dplyr::select(cat, comb1, comb2, P.adj))
})


tbs <- lapply(split(tb, tb$cat), function(x) { x["cat"] <- NULL; x })
df <- tibble(cat = levels(as.factor(unique(x$cat))), 
             tbl = tbs)

p8 <- ggsurvp_cat_item(x, F, 2, "jco", "", "", "HNSC HPV infected Tumor") + geom_table_npc(data = df, aes(npcx = 0.6, npcy = 1, label = tbl),
                                                                                           hjust = 0, vjust = 1,
                                                                                           table.theme = ttheme_gtlight) 
ggsave(plot = p8,
       "paper figures/output/Figure 8.pdf", device = "pdf",
       height = 12, width = 12)

#########* Combination of T/MemB SPANK in infected patients ######
x <- map_dfr(as.list(c("IL17RB", "IL17REL")), function(il17r){
  cell %>%
    dplyr::select(sample, `Memory B`, SPANK) %>%
    mutate(`Memory B` = factor(Hmisc::cut2(`Memory B`, g = 2), labels = c("L", "H"))) %>%
    mutate(SPANK = factor(Hmisc::cut2(SPANK, g = 2), labels = c("L", "H"))) %>%
    pivot_longer(-sample, names_to = "celltype", values_to = "fraction") %>%
    inner_join(gene %>%
                 filter(symbol %in% c(il17r)) %>%
                 split_median("symbol", "raw_count_scaled")) %>%
    clinical_combine(cancer) %>%
    unite("item", c("fraction", "item"), sep = "/", remove = T) %>%
    unite("cat", c("celltype", "symbol"), sep = "/", remove = T) %>%
    inner_join(hpv_inf) %>%
    filter(infection == "infected") %>%
    as.data.frame() %>%
    mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H"))
})


tb <- map_dfr(as.list(unique(x$cat)), function(t){
  y <- x %>% filter(cat == t)
  as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                  p.adjust.method = "none",
                                  data = y)[3]) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname, names_to = "comb2", values_to = "P") %>%
    na.omit() %>%
    mutate(comb2 = gsub("\\.", "/", comb2)) %>%
    mutate(comb2 = gsub(".*value/", "", comb2)) %>%
    mutate_if(is.numeric, round, digits = 3) %>%
    tidybulk::rename("comb1" = "rowname") %>%
    mutate(cat = t,
           test = paste(comb1, comb2, sep = "/")) %>%
    filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
    dplyr::select(cat, comb1, comb2, P) %>%
    inner_join(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                               p.adjust.method = "BH",
                                               data = y)[3]) %>%
                 rownames_to_column() %>%
                 pivot_longer(-rowname, names_to = "comb2", values_to = "P.adj") %>%
                 na.omit() %>%
                 mutate(comb2 = gsub("\\.", "/", comb2)) %>%
                 mutate(comb2 = gsub(".*value/", "", comb2)) %>%
                 mutate_if(is.numeric, round, digits = 3) %>%
                 tidybulk::rename("comb1" = "rowname") %>%
                 mutate(cat = t,
                        test = paste(comb1, comb2, sep = "/")) %>%
                 filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
                 dplyr::select(cat, comb1, comb2, P.adj))
})


tbs <- lapply(split(tb, tb$cat), function(x) { x["cat"] <- NULL; x })
df <- tibble(cat = levels(as.factor(unique(x$cat))), 
             tbl = tbs)

p8 <- ggsurvp_cat_item(x, F, 2, "jco", "", "", "HNSC HPV infected Tumor") + geom_table_npc(data = df, aes(npcx = 0.6, npcy = 1, label = tbl),
                                                                                           hjust = 0, vjust = 1,
                                                                                           table.theme = ttheme_gtlight) 
ggsave(plot = p8,
       "paper figures/output/Figure 8.pdf", device = "pdf",
       height = 12, width = 12)


######* Combination of IL-17RB/IL-17B and IL-17RB/IL-17A, IL17REL/IL-17 *#######
x <- map_dfr(as.list(c("IL17A", "IL17B")), function(il17r){
  gene %>%
    filter(symbol %in% c("IL17RB")) %>%
    split_median("symbol", "raw_count_scaled") %>%
    tidybulk::rename("IL17RB" = "item") %>%
    dplyr::select(sample, IL17RB) %>%
    inner_join(gene %>%
                 filter(symbol %in% c(il17r)) %>%
                 split_median("symbol", "raw_count_scaled") %>%
                 tidybulk::rename(il17r = "item")) %>%
    unite("item", c("IL17RB", il17r), sep = "/", remove = T) %>%
    mutate(cat = paste0("IL17RB/", il17r)) %>%
    inner_join(hpv_inf) %>%
    filter(infection == "infected") 
}) %>%
  dplyr::select(sample, cat, item) %>%
  rbind(map_dfr(as.list(c("IL17A", "IL17C")), function(il17r){
    gene %>%
      filter(symbol %in% c("IL17REL")) %>%
      split_median("symbol", "raw_count_scaled") %>%
      tidybulk::rename("IL17REL" = "item") %>%
      dplyr::select(sample, IL17REL) %>%
      inner_join(gene %>%
                   filter(symbol %in% c(il17r)) %>%
                   split_median("symbol", "raw_count_scaled") %>%
                   tidybulk::rename(il17r = "item")) %>%
      unite("item", c("IL17REL", il17r), sep = "/", remove = T) %>%
      mutate(cat = paste0("IL17REL/", il17r)) %>%
      inner_join(hpv_inf) %>%
      filter(infection == "infected") %>%
      dplyr::select(sample, cat, item)
  })) %>%
  clinical_combine(cancer) %>%
  as.data.frame() %>%
  mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H"))


tb <- map_dfr(as.list(unique(x$cat)), function(t){
  y <- x %>% filter(cat == t)
  as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                  p.adjust.method = "none",
                                  data = y)[3]) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname, names_to = "comb2", values_to = "P") %>%
    na.omit() %>%
    mutate(comb2 = gsub("\\.", "/", comb2)) %>%
    mutate(comb2 = gsub(".*value/", "", comb2)) %>%
    mutate_if(is.numeric, round, digits = 3) %>%
    tidybulk::rename("comb1" = "rowname") %>%
    mutate(cat = t,
           test = paste(comb1, comb2, sep = "/")) %>%
    filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
    dplyr::select(cat, comb1, comb2, P) %>%
    inner_join(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                               p.adjust.method = "BH",
                                               data = y)[3]) %>%
                 rownames_to_column() %>%
                 pivot_longer(-rowname, names_to = "comb2", values_to = "P.adj") %>%
                 na.omit() %>%
                 mutate(comb2 = gsub("\\.", "/", comb2)) %>%
                 mutate(comb2 = gsub(".*value/", "", comb2)) %>%
                 mutate_if(is.numeric, round, digits = 3) %>%
                 tidybulk::rename("comb1" = "rowname") %>%
                 mutate(cat = t,
                        test = paste(comb1, comb2, sep = "/")) %>%
                 filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
                 dplyr::select(cat, comb1, comb2, P.adj))
})


tbs <- lapply(split(tb, tb$cat), function(x) { x["cat"] <- NULL; x })
df <- tibble(cat = levels(as.factor(unique(x$cat))), 
             tbl = tbs)

pa <- ggsurvp_cat_item(x, F, 2, "jco", "", "", "HNSC HPV infected Tumor") + geom_table_npc(data = df, aes(npcx = 0.6, npcy = 1, label = tbl),
                                                                                           hjust = 0, vjust = 1,
                                                                                           table.theme = ttheme_gtlight) 

x <- map_dfr(as.list(c("IL17A", "IL17B")), function(il17r){
  gene %>%
    filter(symbol %in% c("IL17RB")) %>%
    split_median("symbol", "raw_count_scaled") %>%
    tidybulk::rename("IL17RB" = "item") %>%
    dplyr::select(sample, IL17RB) %>%
    inner_join(gene %>%
                 filter(symbol %in% c(il17r)) %>%
                 split_median("symbol", "raw_count_scaled") %>%
                 tidybulk::rename(il17r = "item")) %>%
    unite("item", c("IL17RB", il17r), sep = "/", remove = T) %>%
    mutate(cat = paste0("IL17RB/", il17r)) %>%
    inner_join(hpv_inf) %>%
    filter(infection == "uninfected") 
}) %>%
  dplyr::select(sample, cat, item) %>%
  rbind(map_dfr(as.list(c("IL17A", "IL17C")), function(il17r){
    gene %>%
      filter(symbol %in% c("IL17REL")) %>%
      split_median("symbol", "raw_count_scaled") %>%
      tidybulk::rename("IL17REL" = "item") %>%
      dplyr::select(sample, IL17REL) %>%
      inner_join(gene %>%
                   filter(symbol %in% c(il17r)) %>%
                   split_median("symbol", "raw_count_scaled") %>%
                   tidybulk::rename(il17r = "item")) %>%
      unite("item", c("IL17REL", il17r), sep = "/", remove = T) %>%
      mutate(cat = paste0("IL17REL/", il17r)) %>%
      inner_join(hpv_inf) %>%
      filter(infection == "uninfected") %>%
      dplyr::select(sample, cat, item)
  })) %>%
  clinical_combine(cancer) %>%
  as.data.frame() %>%
  mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H"))


tb <- map_dfr(as.list(unique(x$cat)), function(t){
  y <- x %>% filter(cat == t)
  as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                  p.adjust.method = "none",
                                  data = y)[3]) %>%
    rownames_to_column() %>%
    pivot_longer(-rowname, names_to = "comb2", values_to = "P") %>%
    na.omit() %>%
    mutate(comb2 = gsub("\\.", "/", comb2)) %>%
    mutate(comb2 = gsub(".*value/", "", comb2)) %>%
    mutate_if(is.numeric, round, digits = 3) %>%
    tidybulk::rename("comb1" = "rowname") %>%
    mutate(cat = t,
           test = paste(comb1, comb2, sep = "/")) %>%
    filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
    dplyr::select(cat, comb1, comb2, P) %>%
    inner_join(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ item,
                                               p.adjust.method = "BH",
                                               data = y)[3]) %>%
                 rownames_to_column() %>%
                 pivot_longer(-rowname, names_to = "comb2", values_to = "P.adj") %>%
                 na.omit() %>%
                 mutate(comb2 = gsub("\\.", "/", comb2)) %>%
                 mutate(comb2 = gsub(".*value/", "", comb2)) %>%
                 mutate_if(is.numeric, round, digits = 3) %>%
                 tidybulk::rename("comb1" = "rowname") %>%
                 mutate(cat = t,
                        test = paste(comb1, comb2, sep = "/")) %>%
                 filter(!test %in% c("H/L/L/H", "H/H/L/L", "L/L/H/H", "L/H/H/L")) %>%
                 dplyr::select(cat, comb1, comb2, P.adj))
})


tbs <- lapply(split(tb, tb$cat), function(x) { x["cat"] <- NULL; x })
df <- tibble(cat = levels(as.factor(unique(x$cat))), 
             tbl = tbs)

pb <- ggsurvp_cat_item(x, F, 2, "jco", "", "", "HNSC HPV infected Tumor") + geom_table_npc(data = df, aes(npcx = 0.6, npcy = 1, label = tbl),
                                                                                           hjust = 0, vjust = 1,
                                                                                           table.theme = ttheme_gtlight) 

ggsave(plot = plot_grid(pa, pb, nrow = 1),
       "paper figures/output/suplementary 7c.pdf", device = "pdf",
       height = 14, width = 24) 
