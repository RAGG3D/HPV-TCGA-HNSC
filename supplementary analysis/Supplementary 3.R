source("src/functions.R")
ref = readRDS("D:/PhD Projects/Rprojects/fastTCGA/data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]
colnames(ref) <- (data.frame(Var2 = colnames(ref)) %>%
                    inner_join(read_csv("data/cellname.csv")))$`Cell Type`

cellorder <- c("Endo", "Epi", "Fibro", "Macro M1", "Macro M2", "Mono", "Neutro",
               "Eosino", "Mast", "iDC", "mDC", "Naive B", "Memory B", "ReNK",       
               "IL2NK", "SPANK", "CD4 Tcm", "CD4 Tem", "Naive CD8 T", "CD8 Tcm", "CD8 Tem",  "GD T",
               "T Helper", "Treg" )
cancer = "HNSC"

gene <- TCGA_transcript_tumor(cancer) 
cell <-  CIBERSORT(gene, "sample", "symbol", "raw_count", ref, "cibersort", "get") 

hpv_inf <- read_csv("data/HNSC_HPV.csv") %>%
  right_join(read_csv("data/pid.csv")) %>%
  mutate(hpv = ifelse(is.na(`Max aligns Human Papillomavirus`) == T, 0, `Max aligns Human Papillomavirus`)) %>%
  mutate(infection = ifelse(hpv == 0, "uninfected", "infected")) %>%
  dplyr::select(sample, infection, hpv)


x <- cell %>%
  pivot_longer(-sample, names_to = "cat", values_to = "item") %>%
  inner_join(hpv_inf) %>%
  group_by(infection, cat) %>%
  mutate(item = factor(Hmisc::cut2(item, g = 2), labels = c(1:nlevels(Hmisc::cut2(item, g = 2))))) %>%
  ungroup() %>%
  mutate(item = gsub("1", "L", item), 
         item = gsub("2", "H", item)) %>%
  unite("item", c("infection", "item"), sep = "/", remove = T) %>%
  clinical_combine(cancer) 

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
df <- tibble(cat = levels(as.factor(cellorder)), 
             tbl = tbs)
x <-  x %>%
  mutate(item = gsub("uninfected", "HPV-", item)) %>%
  mutate(item = gsub("infected", "HPV+", item)) %>%
  mutate(item = fct_relevel(item, "HPV-/L", "HPV-/H", "HPV+/L", "HPV+/H" )) %>%
  as.data.frame()

supp3 <- ggsurvp_cat_item(x, F, 4, "jco", "", "Survival probability", "") + 
  geom_table_npc(data = df, aes(npcx = 0.5, npcy = 0.9, label = tbl),
                 hjust = 0, vjust = 1,
                 table.theme = ttheme_gtlight) 


ggsave(plot = supp3, 
       "paper figures/supplementary/Supplementary 3.pdf", device = "pdf", 
       height = 20, width = 30)
