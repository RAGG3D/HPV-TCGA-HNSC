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

c("IL17A", "IL17B", "IL17C", "IL17D", "IL17F", "IL25",
  "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE", "IL17REL")

x <- hpv_inf %>%
  filter(infection == "infected") %>%
  mutate(infection = "HPV load") %>%
  unique() %>%
  mutate(hpv = factor(Hmisc::cut2(hpv, g = 2), labels = c("L", "H"))) %>%
  inner_join(gene %>%
               filter(symbol %in% c("IL17A", "IL17B", "IL17C",
                                    "IL17RA", "IL17RB", "IL17RC", "IL17RE", "IL17REL")) %>%
               dplyr::select(sample, symbol, raw_count_scaled) %>%
               split_median("symbol", "raw_count_scaled") %>%
               tidybulk::rename("IL17R" = "item")) %>%
  clinical_combine(cancer) %>%
  unite("cat", c("infection", "symbol"), sep = "/", remove = T) %>%
  unite("item", c("hpv", "IL17R"), sep = "/", remove = T) %>%
  as.data.frame() %>%
  mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H")) 

tb <- map_dfr(as.list(paste("HPV load", c("IL17A", "IL17B", "IL17C",
                                     "IL17RA", "IL17RB", "IL17RC", "IL17RE", "IL17REL"), sep = "/")), function(t){
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
df <- tibble(cat = levels(as.factor(paste("HPV load", c("IL17A", "IL17B", "IL17C",
                                                        "IL17RA", "IL17RB", "IL17RC", "IL17RE", "IL17REL"), sep = "/"))), 
             tbl = tbs)

p <- ggsurvp_cat_item(x, F, 2, "jco", "", "", "HNSC HPV infected Tumor") + geom_table_npc(data = df, aes(npcx = 0.5, npcy = 1, label = tbl),
                                                                                           hjust = 0, vjust = 1,
                                                                                           table.theme = ttheme_gtlight) +
  labs(tag = "A")

pb01 <- gene %>% 
  filter(symbol %in% c("IL17A")) %>%
  mutate(IL17A = raw_count_scaled) %>%
  cbind(gene %>%
          filter(symbol == "IL17RB") %>%
          mutate(IL17RB = raw_count_scaled))  %>%
  dplyr::select(IL17A, IL17RB) %>%
  ggplot(aes(y = log10(IL17A+1), x = log10(IL17RB +1))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.8)  +
  labs(x = "log10(TMM norm IL17RB)", y = "log10(TMM norm IL17A)", tag = "") +
  theme_bw() +
  theme(text=element_text(size=30, family="sans"),
        legend.position = "bottom") 


pb02 <- gene %>% 
  filter(symbol %in% c("IL17A")) %>%
  mutate(IL17A = raw_count_scaled) %>%
  cbind(gene %>%
               filter(symbol == "IL17REL") %>%
               mutate(IL17REL = raw_count_scaled))  %>%
  dplyr::select(IL17A, IL17REL) %>%
  ggplot(aes(y = log10(IL17A+1), x = log10(IL17REL +1))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.8)  +
  labs(x = "log10(TMM norm IL17REL)", y = "log10(TMM norm IL17A)", tag = "") +
  theme_bw() +
  theme(text=element_text(size=30, family="sans"),
        legend.position = "bottom") 

pb03 <- gene %>% 
  filter(symbol %in% c("IL17RB")) %>%
  mutate(IL17RB = raw_count_scaled) %>%
  cbind(gene %>%
          filter(symbol == "IL17REL") %>%
          mutate(IL17REL = raw_count_scaled))  %>%
  dplyr::select(IL17RB, IL17REL) %>%
  ggplot(aes(y = log10(IL17RB+1), x = log10(IL17REL +1))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.8)  +
  labs(x = "log10(TMM norm IL17REL)", y = "log10(TMM norm IL17RB)", tag = "") +
  theme_bw() +
  theme(text=element_text(size=30, family="sans"),
        legend.position = "bottom") 

pb <- plot_grid(pb01, pb02, pb03, nrow = 1)

ggsave(plot = plot_grid(p, pb, ncol = 1, rel_heights = c(2, 1)),
       "paper figures/supplementary 7.pdf", device = "pdf", height = 17, width = 18)
