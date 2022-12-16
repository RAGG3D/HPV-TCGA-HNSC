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


###### IL-17A, IL-17B, IL-17C, IL-17E
x <- gene %>%
  filter(symbol %in% c("IL17A", "IL17C", "IL17D", "IL17B", "IL25", "IL17F",
                       "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE")) %>%
  inner_join(hpv_inf) %>%
  filter(infection == "infected") %>%
  split_median("symbol", "raw_count_scaled")%>%
  mutate(cat = symbol) %>%
  clinical_combine(cancer) %>%
  as.data.frame() %>%
  mutate(item = fct_relevel(item, "L", "H"))

p01 <- ggsurvp_cat_item(x, T, 3, "jco", "", "Survival probability", "") 

ggsave(plot = p01, "paper figures/output/IL17 KM.pdf", device = "pdf",
       height = 15, width = 15)



########
p02 <- gene %>%
  filter(symbol %in% c("IL17A", "IL17C", "IL17D", "IL17B", "IL25", "IL17F",
                       "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE")) %>%
  inner_join(hpv_inf) %>%
  mutate(infection = ifelse(infection == "uninfected", "HPV-", "HPV+")) %>%
  mutate(infection = fct_relevel(infection, "HPV-", "HPV+")) %>%
  ggplot(aes(x=infection, y=log10(raw_count_scaled + 1), fill=infection)) + 
  geom_boxplot() +
  stat_compare_means(aes(group = infection), label = "p.signif") +
  scale_fill_manual(values = c("#4dbbd5ff", "#e64b35ff")) +
  facet_wrap(.~symbol) +
  theme_bw() + 
  theme(text=element_text(size=30, family="sans"),
        legend.position = "bottom",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
  labs(x = "", y = "IL-17 family expression log10", fill = "", title = "")
ggsave(plot = p02, "paper figures/output/IL17 box.pdf", device = "pdf",
       height = 15, width = 16)
########
p03 <- gene %>%
  filter(symbol %in% c("IL17A", "IL17C", "IL17D", "IL17B", "IL25", "IL17F",
                       "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE")) %>%
  inner_join(hpv_inf) %>%
  filter(infection == "infected") %>%
  ggplot(aes(y = log10(raw_count_scaled + 1), x = log10(hpv +1))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.8) +
  facet_wrap(.~symbol, ncol = 3) +
  labs(x = "HPV viral load (Log10) \n HPV infected HNSC patients", y = "IL-17 family expression log10", tag = "C") +
  theme_bw() +
  theme(text=element_text(size=30, family="sans"),
        legend.position = "bottom") 
ggsave(plot = p03, "paper figures/output/IL17 point.pdf", device = "pdf",
       height = 18, width = 15)




