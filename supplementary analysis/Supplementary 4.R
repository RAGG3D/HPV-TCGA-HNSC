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