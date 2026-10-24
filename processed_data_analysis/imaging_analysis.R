library(readxl)
library(reshape2)
library(emmeans)
library(lme4)
library(lmerTest)
library(ggplotify)
library(ggpubr)
library(sjPlot)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 90, hjust = 1)))


#IF data SAFB
SAFB_meta <- read_excel("/home/shared_folder/imaging_data/SAFB_meta.xlsx")

SAFB_WT <- melt(read_excel("/home/shared_folder/imaging_data/SAFB_WT.xlsx") %>% dplyr::select(-Surfaces))
SAFB_WT = na.omit(SAFB_WT)
colnames(SAFB_WT) = c("cell", "value")
SAFB_meta_WT = SAFB_meta %>% dplyr::filter(condition == "WT")
SAFB_WT = SAFB_WT %>% dplyr::left_join(SAFB_meta_WT, by = "cell")

SAFB_RBM20KO <- melt(read_excel("/home/shared_folder/imaging_data/SAFB_RBM20KO.xlsx") %>% dplyr::select(-Surfaces))
SAFB_RBM20KO = na.omit(SAFB_RBM20KO)
colnames(SAFB_RBM20KO) = c("cell", "value")
SAFB_meta_RBM20KO = SAFB_meta %>% dplyr::filter(condition == "RBM20KO")
SAFB_RBM20KO = SAFB_RBM20KO %>% dplyr::left_join(SAFB_meta_RBM20KO, by = "cell")

SAFB_pTTNKO <- melt(read_excel("/home/shared_folder/imaging_data/SAFB_pTTNKO.xlsx") %>% dplyr::select(-Surfaces))
SAFB_pTTNKO = na.omit(SAFB_pTTNKO)
colnames(SAFB_pTTNKO) = c("cell", "value")
SAFB_meta_pTTNKO = SAFB_meta %>% dplyr::filter(condition == "pTTNKO")
SAFB_pTTNKO = SAFB_pTTNKO %>% dplyr::left_join(SAFB_meta_pTTNKO, by = "cell")

SAFB_doubleKO <- melt(read_excel("/home/shared_folder/imaging_data/SAFB_doubleKO.xlsx") %>% dplyr::select(-Surfaces))
SAFB_doubleKO = na.omit(SAFB_doubleKO)
colnames(SAFB_doubleKO) = c("cell", "value")
SAFB_meta_doubleKO = SAFB_meta %>% dplyr::filter(condition == "doubleKO")
SAFB_doubleKO = SAFB_doubleKO %>% dplyr::left_join(SAFB_meta_doubleKO, by = "cell")

SAFB_all = rbind(SAFB_WT, SAFB_RBM20KO, SAFB_pTTNKO, SAFB_doubleKO)

n_foci_per_cell = c()
for(i in SAFB_all$value){
  filt = SAFB_all %>% dplyr::filter(value >= i)
  sum = filt %>% dplyr::mutate(n = 1) %>% dplyr::group_by(cell) %>% dplyr::summarise(foci = sum(n))
  n_foci_per_cell = append(n_foci_per_cell, sum(sum$foci)/length(unique(SAFB_all$cell)))
}

SAFB_all$n_foci_per_cell = n_foci_per_cell
n_foci_SAFB = SAFB_all %>% dplyr::mutate(n = 1) %>% dplyr::group_by(condition, diff, cell) %>% dplyr::summarise(n_foci = sum(n))
n_foci_SAFB$diff = factor(n_foci_SAFB$diff, levels = c("Diff1", "Diff2"))
n_foci_SAFB$condition = factor(n_foci_SAFB$condition, levels = c("WT", "RBM20KO", "pTTNKO", "doubleKO"))
model_n_foci_SAFB <- lmer(n_foci ~ condition + (1 | diff), data = n_foci_SAFB)
emm <- emmeans(model_n_foci_SAFB, pairwise ~ condition)
emm_df =  as.data.frame(emm[["contrasts"]])
g1 = c()
g2 = c()
for(i in emm_df$contrast){
  g1 = append(g1, strsplit(i, " - ")[[1]][1])
  g2 = append(g2, strsplit(i, " - ")[[1]][2])
}
pval_df = data.frame(group1 = g1, group2 = g2, p = round(emm_df$p.value, digits = 5), y.position = seq(from = max(n_foci_SAFB$n_foci) + 0.5, to = max(n_foci_SAFB$n_foci) + 11.5, by = 2))
pval_df = pval_df %>% dplyr::mutate(p = ifelse(p < 0.00001, "< 0.00001 ****",
                                               ifelse(p <= 0.0001, paste0(format(p, scientific = FALSE), " ***"),
                                                      ifelse(p <= 0.001, paste0(format(p, scientific = FALSE), " **"),
                                                             ifelse(p <= 0.05, paste0(format(p, scientific = FALSE), " *"),
                                                                    paste0(format(p, scientific = FALSE), " ns"))))))
ggplot(n_foci_SAFB, aes(x = condition, y = n_foci)) + geom_violin() + geom_boxplot(fill = "#DEA0FD", alpha = 0.5) +
  stat_pvalue_manual(data = pval_df, size = 3, coord.flip = TRUE) + xlab("") + ylab("number of foci per cell") + coord_flip()
ggsave("/home/shared_folder/output_plots/SAFB_nfoci_boxplot.pdf")

#IF data QKI
QKI_meta <- read_excel("/home/shared_folder/imaging_data/QKI_meta.xlsx")

QKI_WT <- melt(read_excel("/home/shared_folder/imaging_data/QKI_WT.xlsx") %>% dplyr::select(-Surfaces))
QKI_WT = na.omit(QKI_WT)
colnames(QKI_WT) = c("cell", "value")
QKI_meta_WT = QKI_meta %>% dplyr::filter(condition == "WT")
QKI_WT = QKI_WT %>% dplyr::left_join(QKI_meta_WT, by = "cell")

QKI_RBM20KO <- melt(read_excel("/home/shared_folder/imaging_data/QKI_RBM20KO.xlsx") %>% dplyr::select(-Surfaces))
QKI_RBM20KO = na.omit(QKI_RBM20KO)
colnames(QKI_RBM20KO) = c("cell", "value")
QKI_meta_RBM20KO = QKI_meta %>% dplyr::filter(condition == "RBM20KO")
QKI_RBM20KO = QKI_RBM20KO %>% dplyr::left_join(QKI_meta_RBM20KO, by = "cell")

QKI_pTTNKO <- melt(read_excel("/home/shared_folder/imaging_data/QKI_pTTNKO.xlsx") %>% dplyr::select(-Surfaces))
QKI_pTTNKO = na.omit(QKI_pTTNKO)
colnames(QKI_pTTNKO) = c("cell", "value")
QKI_meta_pTTNKO = QKI_meta %>% dplyr::filter(condition == "pTTNKO")
QKI_pTTNKO = QKI_pTTNKO %>% dplyr::left_join(QKI_meta_pTTNKO, by = "cell")

QKI_doubleKO <- melt(read_excel("/home/shared_folder/imaging_data/QKI_doubleKO.xlsx") %>% dplyr::select(-Surfaces))
QKI_doubleKO = na.omit(QKI_doubleKO)
colnames(QKI_doubleKO) = c("cell", "value")
QKI_meta_doubleKO = QKI_meta %>% dplyr::filter(condition == "doubleKO")
QKI_doubleKO = QKI_doubleKO %>% dplyr::left_join(QKI_meta_doubleKO, by = "cell")

QKI_all = rbind(QKI_WT, QKI_RBM20KO, QKI_pTTNKO, QKI_doubleKO)

n_foci_per_cell = c()
for(i in QKI_all$value){
  filt = QKI_all %>% dplyr::filter(value >= i)
  sum = filt %>% dplyr::mutate(n = 1) %>% dplyr::group_by(cell) %>% dplyr::summarise(foci = sum(n))
  n_foci_per_cell = append(n_foci_per_cell, sum(sum$foci)/length(unique(QKI_all$cell)))
}

QKI_all$n_foci_per_cell = n_foci_per_cell
n_foci_QKI = QKI_all %>% dplyr::mutate(n = 1) %>% dplyr::group_by(condition, diff, cell) %>% dplyr::summarise(n_foci = sum(n))
n_foci_QKI$diff = factor(n_foci_QKI$diff, levels = c("Diff1", "Diff2"))
n_foci_QKI$condition = factor(n_foci_QKI$condition, levels = c("WT", "RBM20KO", "pTTNKO", "doubleKO"))
model_n_foci_QKI <- lmer(n_foci ~ condition + (1 | diff), data = n_foci_QKI)
emm <- emmeans(model_n_foci_QKI, pairwise ~ condition)
emm_df =  as.data.frame(emm[["contrasts"]])
g1 = c()
g2 = c()
for(i in emm_df$contrast){
  g1 = append(g1, strsplit(i, " - ")[[1]][1])
  g2 = append(g2, strsplit(i, " - ")[[1]][2])
}
pval_df = data.frame(group1 = g1, group2 = g2, p = round(emm_df$p.value, digits = 5), y.position = seq(from = max(n_foci_QKI$n_foci) + 0.5, to = max(n_foci_QKI$n_foci) + 11.5, by = 2))
pval_df = pval_df %>% dplyr::mutate(p = ifelse(p < 0.00001, "< 0.00001 ****",
                                               ifelse(p <= 0.0001, paste0(format(p, scientific = FALSE), " ***"),
                                                      ifelse(p <= 0.001, paste0(format(p, scientific = FALSE), " **"),
                                                             ifelse(p <= 0.05, paste0(format(p, scientific = FALSE), " *"),
                                                                    paste0(format(p, scientific = FALSE), " ns"))))))
ggplot(n_foci_QKI, aes(x = condition, y = n_foci)) + geom_violin() + geom_boxplot(fill = "#009E73", alpha = 0.5) +
  stat_pvalue_manual(data = pval_df, size = 3, coord.flip = TRUE) + xlab("") + ylab("number of foci per cell") + coord_flip()
ggsave("/home/shared_folder/output_plots/QKI_nfoci_boxplot.pdf")
