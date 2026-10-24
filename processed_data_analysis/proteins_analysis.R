library(readxl)
library(reshape2)
library(ggplotify)
library(ggpubr)
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

##proteomics data GSEA
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

#core genes in pulldown
res_i = res_gsea_BP %>% dplyr::filter(Description == i)
core_genes = res_i$core_enrichment
core_genes = strsplit(core_genes, "/")[[1]]
ggplot(proteomics_pulldown %>% dplyr::filter(Gene %in% core_genes) %>% dplyr::mutate(SAFB = ifelse(Gene == "SAFB", "SAFB", "no"), sign = ifelse(qvalue < 0.05, Gene, "")), aes(x = -Log2FC_KO_WT, y = -log10(qvalue))) +
  geom_point(aes(fill = SAFB, color = SAFB), size = 5, shape = 21) +  scale_fill_manual(values = c("grey", "magenta")) + scale_color_manual(values = c("grey", "black")) +
  ggtitle(i) + theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0))) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = sign), max.overlaps = 5,	size = 3, hjust = 2) + xlab("log2 fold change WT/RBM20KO") 
ggsave("/home/shared_folder/output_plots/volcano_plot_SAFB.pdf")


#SAFB KD
SAFB_dge <- read.csv("/home/shared_folder/gene_expression_data/DGE_KD_data/SAFB.csv", row.names=1) %>%
  dplyr::select(c(gene_ID, gene_name, logFC))
colnames(SAFB_dge) = c("gene_ID", "gene_name", "log2FC_SAFB_KD")

down_SAFB = (SAFB_dge %>% dplyr::filter(logFC < -2 & adj.P.Val < 0.05))$gene_ID
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

factory_genes_RUES2 = factory_genes_RUES2 %>%
dplyr::mutate(SAFB_KD = ifelse(gene_ID %in% (SAFB_dge %>% dplyr::filter(col_to_plot == "down_SCR"))$gene_ID, "down_SAFB_KD", ifelse(gene_ID %in% (SAFB_dge %>% dplyr::filter(col_to_plot == "up_SAFB_KD"))$gene_ID, "up_SAFB_KD", "")))

ggplot() +
  geom_point(aes(x = factory_genes_RUES2$logFC_RBM20KO, y = factory_genes_RUES2$logFC_TTNpromKO), size = 2, color = "grey", alpha = 0.5) +
  stat_density_2d(aes(x = (factory_genes_RUES2 %>% dplyr::filter(SAFB_KD == "down_SAFB_KD"))$logFC_RBM20KO, y = (factory_genes_RUES2 %>% dplyr::filter(SAFB_KD == "down_SAFB_KD"))$logFC_TTNpromKO), color = "#6db6ff") +
  stat_density_2d(aes(x = (factory_genes_RUES2 %>% dplyr::filter(SAFB_KD == "up_SAFB_KD"))$logFC_RBM20KO, y = (factory_genes_RUES2 %>% dplyr::filter(SAFB_KD == "up_SAFB_KD"))$logFC_TTNpromKO), color = "#ff6db6") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_x_continuous(limits = c(-5, 5)) + scale_y_continuous(limits = c(-5, 5)) + xlab("RBM20KO vs WT logFC") + ylab("TTNpromKO vs WT logFC")
ggsave("/home/shared_folder/output_plots/SAFB_dge_factory_density.pdf")

#TIDs GSEA
act_genes <- read.csv("/home/shared_folder/chromatin_data/ChromHMM_data/act_genes.csv", row.names=1)
rep_genes <- read.csv("/home/shared_folder/chromatin_data/ChromHMM_data/DGE_RUES2_data/rep_genes.csv", row.names=1)
biv_genes <- read.csv("/home/shared_folder/chromatin_data/ChromHMM_data/DGE_RUES2_data/biv_genes.csv", row.names=1)

pathway = list(tids = c(act_genes$geneID, rep_genes$geneID, biv_genes$geneID))

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

