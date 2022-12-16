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

######* Panel A: Plot survival infected vs non-infected ######
x <- hpv_inf %>%
  clinical_combine(cancer) %>%
  mutate(cat = cancer, item = infection) %>%  
  mutate(item = ifelse(item == "uninfected", "HPV-", "HPV+")) %>%
  mutate(item = fct_relevel(item, "HPV-", "HPV+")) %>%
  as.data.frame()

tb <- data.frame(P = c(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ cat + item,p.adjust.method = "none",
                                                       data = x)[3])[1,1]),
                 adj.P = c(as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ cat + item,p.adjust.method = "BH",
                                                           data = x)[3])[1,1])) %>% 
  mutate_if(is.numeric, round, digits = 3) %>%
  mutate(cancer = cancer)
tbs <- lapply(split(tb, tb$cancer), function(x) { x["cancer"] <- NULL; x })
df <- tibble(cat = cancer, 
             tbl = tbs)

p1a <- ggsurvp_cat_item(x, F, 1, "jco", "", "", "") + 
  geom_table_npc(data = df, aes(npcx = 0.75, npcy = 0.8, label = tbl),
                 hjust = 0, vjust = 1,
                 table.theme = ttheme_gtlight) +
  labs(tag = "A")
  
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
  labs(x = "", y = "Immune Infiltration Logit", fill = "", title = "", tag = "B")


p1b02 <- gene %>%
  filter(symbol %in% c("PTPRC")) %>%
  inner_join(hpv_inf) %>%
  mutate(infection = ifelse(infection == "uninfected", "HPV-", "HPV+")) %>%
  mutate(infection = fct_relevel(infection, "HPV-", "HPV+")) %>%
  ggplot(aes(x = symbol, y = log10(raw_count_scaled + 1), fill = infection)) +
  geom_boxplot()+
  stat_compare_means(aes(group = infection), label = "p.signif") +
  scale_fill_manual(values = c("#4dbbd5ff", "#e64b35ff")) +
  theme_bw() + 
  scale_y_log10() +
  theme(text=element_text(size=20, family="sans"),
        legend.position = "bottom",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
  labs(x = "", y = " TMM normalized Log10", fill = "", title = "", tag = "")


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
  clinical_combine(cancer) %>%
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
                 table.theme = ttheme_gtlight) +
  labs(tag = "C")

###############
#### Box and scatter plots
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

#Only infected patients: tar_gene expression and viral load
p1e <- x %>%
  filter(infection == "HPV+") %>%
  ggplot(aes(y = logit(item), x = log10(hpv +1))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.8) +
  labs(x = "HPV viral load (Log10) \n HPV infected HNSC patients", y = "logit(Immune TS)", tag = "E") +
  theme_bw() +
  theme(text=element_text(size=30, family="sans"),
        legend.position = "bottom") 


########
#### Box and scatter plots
x <- gene %>%
  filter(symbol == "PTPRC") %>%
  inner_join(hpv_inf) 

#Only infected patients: tar_gene expression and viral load
p1f <- x %>%
  filter(infection == "infected") %>%
  ggplot(aes(y = log10(raw_count_scaled+1), x = log10(hpv +1))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.8) +
  labs(x = "HPV viral load (Log10) \n HPV infected HNSC patients", y = paste0("", " (TMM normalized Log10)"),
       tag = " ") +
  theme_bw() +
  theme(text=element_text(size=30, family="sans"),
        legend.position = "bottom") 


ggsave(plot = plot_grid(p1a, p1b, p1c, nrow = 1, rel_widths = c(1,1.5,1), align = "hv"),
       "paper figures/output/Figure 1.pdf", height = 7, width = 22)


ggsave(plot = plot_grid(p1a, plot_grid(p1e, p1f, nrow = 1) , p1c, nrow = 1, rel_widths = c(1,1.6,1), align = "hv"),
       "paper figures/output/Figure 1.pdf", height = 6.5, width = 22)



