library(readxl)
library(reshape2)
library(emmeans)
library(lme4)
library(lmerTest)
library(ggplotify)
library(ggpubr)
library(sjPlot)
library(clusterProfiler)
library(ensembldb)
library("org.Hs.eg.db", character.only = TRUE)
library(viridis)
library(VennDiagram)
library(fgsea)
library(pheatmap)
library(limma)
library(topGO)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 90, hjust = 1)))

##proteomics data GSEA - figure 6 panel C
proteomics_pulldown <- read_excel("/home/shared_folder/pulldown_data/proteomics_data/proteomics_pulldown.xlsx", skip = 2)[,1:19]
proteomics_pulldown$Log2FC_KO_WT = as.numeric(proteomics_pulldown$Log2FC_KO_WT)
proteomics_pulldown$T_significance[is.na(proteomics_pulldown$T_significance)] = ""

genes = -proteomics_pulldown$T
names(genes) = proteomics_pulldown$Gene
genes = sort(genes, decreasing = TRUE)
gse = gseGO(geneList = genes, 
            ont = "BP", 
            keyType = "SYMBOL",
            pvalueCutoff = 1, 
            OrgDb = org.Hs.eg.db, 
            pAdjustMethod = "BH",
            by = "fgsea")
res_gsea_BP = gse@result
BP_terms_WT = c("DNA repair", "chromatin organization")

res_gsea_BP$n = 1:nrow(res_gsea_BP)
for(i in BP_terms_WT){
  res_i = res_gsea_BP %>% dplyr::filter(Description == i)
  n = res_i$n
  gseaplot(gse, by ="runningScore", color.line =  "#3283FE", title = paste0(gse$Description[n], " p = ", signif(gse$p.adjust[n]), digits = 3), geneSetID = n) +
    theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7)) + theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
  ggsave(paste0("/home/shared_folder/output_plots/runningScore_", gse$Description[n], ".pdf"))
  core_genes = res_i$core_enrichment
  core_genes = strsplit(core_genes, "/")[[1]]
}

#core genes in pulldown - figure 6 panel D
res_i = res_gsea_BP %>% dplyr::filter(Description == i)
core_genes = res_i$core_enrichment
core_genes = strsplit(core_genes, "/")[[1]]
ggplot(proteomics_pulldown %>% dplyr::filter(Gene %in% core_genes) %>% dplyr::mutate(SAFB = ifelse(Gene == "SAFB", "SAFB", "no"), sign = ifelse(qvalue < 0.05, Gene, "")), aes(x = -Log2FC_KO_WT, y = -log10(qvalue))) +
  geom_point(aes(fill = SAFB, color = SAFB), size = 5, shape = 21) +  scale_fill_manual(values = c("grey", "magenta")) + scale_color_manual(values = c("grey", "black")) +
  ggtitle(i) + theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = sign), max.overlaps = 5,	size = 3, hjust = 2) + xlab("log2 fold change WT/RBM20KO") 
ggsave("/home/shared_folder/output_plots/volcano_plot_SAFB.pdf")

#IF data SAFB - figure 7 panel A
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

#IF data QKI - figure 7 panel B
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


#panel G
ELAVL1_dge <- read.csv("/home/shared_folder/gene_expression_data/DGE_KD_data/KD_data/ELAVL1.csv", row.names=1) %>%
  dplyr::select(c(gene_ID, gene_name, logFC))
colnames(ELAVL1_dge) = c("gene_ID", "gene_name", "log2FC_ELAVL1_KD")
QKI_dge <- read.csv("/home/shared_folder/gene_expression_data/DGE_KD_data/KD_data/QKI.csv", row.names=1) %>%
  dplyr::select(c(gene_ID, gene_name, logFC))
colnames(QKI_dge) = c("gene_ID", "gene_name", "log2FC_QKI_KD")
SAFB_dge <- read.csv("/home/shared_folder/gene_expression_data/DGE_KD_data/KD_data/SAFB.csv", row.names=1) %>%
  dplyr::select(c(gene_ID, gene_name, logFC))
colnames(SAFB_dge) = c("gene_ID", "gene_name", "log2FC_SAFB_KD")

DGE_all = ELAVL1_dge %>% dplyr::left_join(QKI_dge, by = c("gene_name", "gene_ID"))  %>% dplyr::left_join(SAFB_dge, by = c("gene_name", "gene_ID")) %>%
  dplyr::inner_join(QKI_KO_RNAseq, by = c("gene_name", "gene_ID"))

down_ELAVL1 = (ELAVL1_dge %>% dplyr::filter(logFC < -2 & adj.P.Val < 0.05))$gene_ID
down_QKI = (QKI_dge %>% dplyr::filter(logFC < -2 & adj.P.Val < 0.05))$gene_ID
down_SAFB = (SAFB_dge %>% dplyr::filter(logFC < -2 & adj.P.Val < 0.05))$gene_ID

up_ELAVL1 = (SE.ELAVL1_KD %>% dplyr::filter(delta_PSI_ELAVL1KD > 0.1 & FDR_ELAVL1KD < 0.05))$variant_ID
up_QKI = (SE.QKI_KD %>% dplyr::filter(delta_PSI_QKIKD > 0.1 & FDR_QKIKD < 0.05))$variant_ID
up_SAFB = (SE.SAFB_KD %>% dplyr::filter(delta_PSI_SAFBKD > 0.1 & FDR_SAFBKD < 0.05))$variant_ID

WTvsRBM20.KO <- read.csv("/home/shared_folder/gene_expression_data/DGE_KD_data/DGE_RUES2_data/WTvsRBM20.KO.csv", row.names=1)
WTvsTTN.PROM.KO <- read.csv("/home/shared_folder/gene_expression_data/DGE_KD_data/DGE_RUES2_data/WTvsTTN.PROM.KO.csv", row.names=1)
WTvsTTN.KO <- read.csv("/home/shared_folder/gene_expression_data/DGE_KD_data/DGE_RUES2_data/WTvsTTN.KO.csv", row.names=1)

factory_genes_RUES2 = WTvsRBM20.KO %>% dplyr::select(c("gene_name", "gene_ID", "logFC", "adj.P.Val")) %>% dplyr::mutate(logFC = -logFC)
colnames(factory_genes_RUES2) = c("gene_name", "gene_ID", "logFC_RBM20KO", "adj.P.Val_RBM20KO")
ttnpromKO_genes = WTvsTTN.PROM.KO %>% dplyr::select(c("gene_name", "gene_ID", "logFC", "adj.P.Val")) %>% dplyr::mutate(logFC = -logFC)
colnames(ttnpromKO_genes) = c("gene_name", "gene_ID", "logFC_TTNpromKO", "adj.P.Val_TTNpromKO")
factory_genes_RUES2 = factory_genes_RUES2 %>% dplyr::left_join(ttnpromKO_genes, by = c("gene_name", "gene_ID"))
ttnKO_genes = WTvsTTN.KO %>% dplyr::select(c("gene_name", "gene_ID", "logFC", "adj.P.Val")) %>% dplyr::mutate(logFC = -logFC)
colnames(ttnKO_genes) = c("gene_name", "gene_ID", "logFC_TTNKO", "adj.P.Val_TTNKO")
factory_genes_RUES2 = factory_genes_RUES2 %>% dplyr::left_join(ttnKO_genes, by = c("gene_name", "gene_ID"))

factory_genes_RUES2 = factory_genes_RUES2 %>% dplyr::mutate(ELAVL1_KD = ifelse(gene_ID %in% (ELAVL1_dge %>% dplyr::filter(col_to_plot == "down_SCR"))$gene_ID, "down_ELAVL1_KD", ifelse(gene_ID %in% (ELAVL1_dge %>% dplyr::filter(col_to_plot == "up_ELAVL1_KD"))$gene_ID, "up_ELAVL1_KD", ""))) %>%
  dplyr::mutate(QKI_KD = ifelse(gene_ID %in% (QKI_dge %>% dplyr::filter(col_to_plot == "down_SCR"))$gene_ID, "down_QKI_KD", ifelse(gene_ID %in% (QKI_dge %>% dplyr::filter(col_to_plot == "up_QKI_KD"))$gene_ID, "up_QKI_KD", ""))) %>%
  dplyr::mutate(SAFB_KD = ifelse(gene_ID %in% (SAFB_dge %>% dplyr::filter(col_to_plot == "down_SCR"))$gene_ID, "down_SAFB_KD", ifelse(gene_ID %in% (SAFB_dge %>% dplyr::filter(col_to_plot == "up_SAFB_KD"))$gene_ID, "up_SAFB_KD", "")))

#factory KD anticorrelation - figure 7 panel C
ggplot() +
  geom_point(aes(x = factory_genes_RUES2$logFC_RBM20KO, y = factory_genes_RUES2$logFC_TTNpromKO), size = 2, color = "grey", alpha = 0.5) +
  stat_density_2d(aes(x = (factory_genes_RUES2 %>% dplyr::filter(ELAVL1_KD == "down_ELAVL1_KD"))$logFC_RBM20KO, y = (factory_genes_RUES2 %>% dplyr::filter(ELAVL1_KD == "down_ELAVL1_KD"))$logFC_TTNpromKO), color = "#6db6ff") +
  stat_density_2d(aes(x = (factory_genes_RUES2 %>% dplyr::filter(ELAVL1_KD == "up_ELAVL1_KD"))$logFC_RBM20KO, y = (factory_genes_RUES2 %>% dplyr::filter(ELAVL1_KD == "up_ELAVL1_KD"))$logFC_TTNpromKO), color = "#ff6db6") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_x_continuous(limits = c(-5, 5)) + scale_y_continuous(limits = c(-5, 5)) + xlab("RBM20KO vs WT logFC") + ylab("TTNpromKO vs WT logFC")
ggsave("/home/shared_folder/output_plots/ELAVL1_dge_factory_density.pdf")

ggplot() +
  geom_point(aes(x = factory_genes_RUES2$logFC_RBM20KO, y = factory_genes_RUES2$logFC_TTNpromKO), size = 2, color = "grey", alpha = 0.5) +
  stat_density_2d(aes(x = (factory_genes_RUES2 %>% dplyr::filter(QKI_KD == "down_QKI_KD"))$logFC_RBM20KO, y = (factory_genes_RUES2 %>% dplyr::filter(QKI_KD == "down_QKI_KD"))$logFC_TTNpromKO), color = "#6db6ff") +
  stat_density_2d(aes(x = (factory_genes_RUES2 %>% dplyr::filter(QKI_KD == "up_QKI_KD"))$logFC_RBM20KO, y = (factory_genes_RUES2 %>% dplyr::filter(QKI_KD == "up_QKI_KD"))$logFC_TTNpromKO), color = "#ff6db6") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_x_continuous(limits = c(-5, 5)) + scale_y_continuous(limits = c(-5, 5)) + xlab("RBM20KO vs WT logFC") + ylab("TTNpromKO vs WT logFC")
ggsave("/home/shared_folder/output_plots/QKI_dge_factory_density.pdf")

ggplot() +
  geom_point(aes(x = factory_genes_RUES2$logFC_RBM20KO, y = factory_genes_RUES2$logFC_TTNpromKO), size = 2, color = "grey", alpha = 0.5) +
  stat_density_2d(aes(x = (factory_genes_RUES2 %>% dplyr::filter(SAFB_KD == "down_SAFB_KD"))$logFC_RBM20KO, y = (factory_genes_RUES2 %>% dplyr::filter(SAFB_KD == "down_SAFB_KD"))$logFC_TTNpromKO), color = "#6db6ff") +
  stat_density_2d(aes(x = (factory_genes_RUES2 %>% dplyr::filter(SAFB_KD == "up_SAFB_KD"))$logFC_RBM20KO, y = (factory_genes_RUES2 %>% dplyr::filter(SAFB_KD == "up_SAFB_KD"))$logFC_TTNpromKO), color = "#ff6db6") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_x_continuous(limits = c(-5, 5)) + scale_y_continuous(limits = c(-5, 5)) + xlab("RBM20KO vs WT logFC") + ylab("TTNpromKO vs WT logFC")
ggsave("/home/shared_folder/output_plots/SAFB_dge_factory_density.pdf")

#TIDs GSEA - figure 7 panel D
act_genes <- read.csv("/home/shared_folder/chromatin_data/ChromHMM_data/act_genes.csv", row.names=1)
rep_genes <- read.csv("/home/shared_folder/chromatin_data/ChromHMM_data/DGE_RUES2_data/rep_genes.csv", row.names=1)
biv_genes <- read.csv("/home/shared_folder/chromatin_data/ChromHMM_data/DGE_RUES2_data/biv_genes.csv", row.names=1)

pathway = list(tids = c(act_genes$geneID, rep_genes$geneID, biv_genes$geneID))

rank_ELAVL1 = ELAVL1_dge$logFC
names(rank_ELAVL1) = ELAVL1_dge$gene_ID
fgsea_ELAVL1 <- fgsea(pathways = pathway,
                      stats = rank_ELAVL1,
                      scoreType = 'std',
                      minSize = 0,
                      maxSize = Inf)
p = plotEnrichment(pathway2[["tids"]], rank_ELAVL1) +
  labs(title = paste0("TIDs enclosed genes ranking enrichment ELAVL1 KD vs SCR\nNES = ", signif(fgsea_ELAVL1$NES[1], digits = 3), "\np-val = ", signif(fgsea_ELAVL1$pval[1], digits = 3)))
p[["layers"]][[1]][["aes_params"]]$colour = "#FEAF16"
as.ggplot(print(p))
ggsave("/home/shared_folder/output_plots/ELAVL1_GSEA_tids_all.pdf", width = 8, height = 5)

rank_QKI = QKI_dge$logFC
names(rank_QKI) = QKI_dge$gene_ID
fgsea_QKI <- fgsea(pathways = pathway,
                   stats = rank_QKI,
                   scoreType = 'std',
                   minSize = 0,
                   maxSize = Inf)
p = plotEnrichment(pathway2[["tids"]], rank_QKI) +
  labs(title = paste0("TIDs enclosed genes ranking enrichment QKI KD vs SCR\nNES = ", signif(fgsea_QKI$NES[1], digits = 3), "\np-val = ", signif(fgsea_QKI$pval[1], digits = 3)))
p[["layers"]][[1]][["aes_params"]]$colour = "#009E73"
as.ggplot(print(p))
ggsave("/home/shared_folder/output_plots/QKI_GSEA_tids_all.pdf", width = 8, height = 5)

rank_SAFB = SAFB_dge$logFC
names(rank_SAFB) = SAFB_dge$gene_ID
plotGseaTable(pathway, stats = rank_SAFB, fgseaRes = fgsea_SAFB, gseaParam = 0.5)
dev.off()
fgsea_SAFB <- fgsea(pathways = pathway,
                    stats = rank_SAFB,
                    scoreType = 'std',
                    minSize = 0,
                    maxSize = Inf)
p = plotEnrichment(pathway2[["tids"]], rank_SAFB) +
  labs(title = paste0("TIDs enclosed genes ranking enrichment SAFB KD vs SCR\nNES = ", signif(fgsea_SAFB$NES[1], digits = 3), "\np-val = ", signif(fgsea_SAFB$pval[1], digits = 3)))
p[["layers"]][[1]][["aes_params"]]$colour = "#DEA0FD"
as.ggplot(print(p))
ggsave("/home/shared_folder/output_plots/SAFB_GSEA_tids_all.pdf", width = 8, height = 5)

