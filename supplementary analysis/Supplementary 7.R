source("src/functions.R")
ref = readRDS("D:/PhD Projects/Rprojects/fastTCGA/data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]
colnames(ref) <- (data.frame(Var2 = colnames(ref)) %>%
                    inner_join(read_csv("data/cellname.csv")))$`Cell Type`
cancer = "HNSC"

gene <- TCGA_transcript_tumor(cancer) 
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

##### KM combination of IL17REL/IL17RB and immune infiltration
x <- hpv_inf %>%
  filter(infection == "infected") %>%
  mutate(infection = "HPV load") %>%
  unique() %>%
  mutate(hpv = factor(Hmisc::cut2(hpv, g = 2), labels = c("L", "H"))) %>%
  inner_join(gene %>%
               filter(symbol %in% up_inf_gf_cyto_horm$symbol) %>%
               dplyr::select(sample, symbol, raw_count_scaled) %>%
               split_median("symbol", "raw_count_scaled") %>%
               tidybulk::rename("IL17R" = "item")) %>%
  clinical_combine(cancer) %>%
  unite("cat", c("infection", "symbol"), sep = "/", remove = T) %>%
  unite("item", c("hpv", "IL17R"), sep = "/", remove = T) %>%
  as.data.frame() %>%
  mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H")) 

tb <- map_dfr(as.list(paste("HPV load", up_inf_gf_cyto_horm$symbol, sep = "/")), function(t){
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
df <- tibble(cat = levels(as.factor(paste("HPV load", up_inf_gf_cyto_horm$symbol, sep = "/"))), 
             tbl = tbs)

supp2 <- ggsurvp_cat_item(x, F, 2, "jco", "", "", "HNSC HPV infected Tumor") + geom_table_npc(data = df, aes(npcx = 0.5, npcy = 1, label = tbl),
                                                                                              hjust = 0, vjust = 1,
                                                                                              table.theme = ttheme_gtlight) 
ggsave(plot = supp2,
       "paper figures/supplementary/Supplementary 2.pdf", device = "pdf",
       height = 10, width = 24)
