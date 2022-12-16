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

##### KM combination of IL17REL/IL17RB and HPV viral loads
x <- hpv_inf %>%
  filter(infection == "infected") %>%
  mutate(infection = "HPV load") %>%
  unique() %>%
  mutate(hpv = factor(Hmisc::cut2(hpv, g = 2), labels = c("L", "H"))) %>%
  inner_join(gene %>%
               filter(symbol %in% c("IL17REL", "IL17RB")) %>%
               dplyr::select(sample, symbol, raw_count_scaled) %>%
               split_median("symbol", "raw_count_scaled") %>%
               tidybulk::rename("IL17R" = "item")) %>%
  clinical_combine(cancer) %>%
  unite("cat", c("infection", "symbol"), sep = "/", remove = T) %>%
  unite("item", c("hpv", "IL17R"), sep = "/", remove = T) %>%
  as.data.frame() %>%
  mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H")) 

tb <- map_dfr(as.list(c("HPV load/IL17REL", "HPV load/IL17RB")), function(t){
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
df <- tibble(cat = levels(as.factor(c("HPV load/IL17REL", "HPV load/IL17RB"))), 
             tbl = tbs)

p5 <- ggsurvp_cat_item(x, F, 1, "jco", "", "", "HNSC HPV infected Tumor") + geom_table_npc(data = df, aes(npcx = 0.5, npcy = 1, label = tbl),
                                                                                           hjust = 0, vjust = 1,
                                                                                           table.theme = ttheme_gtlight) 
ggsave(plot = p5,
       "paper figures/output/Figure 5.pdf", device = "pdf",
       height = 6.5, width = 10)
