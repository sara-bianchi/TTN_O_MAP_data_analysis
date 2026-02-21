library(dplyr)
library(ggplot2)
library(reshape2)
library(readxl)
library(VennDiagram)
library(gridExtra)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size=15, face="bold", margin = margin(10,0,10,0))))

##O-MAP pulldown transcripts
#volcano plots - supplemntary figure 6 panel A
o_map_WT_KO <- read_excel("/home/shared_folder/pulldown_data/transcripts_data/o_map_WT_KO.xlsx")
o_map_WT_KO = o_map_WT_KO %>% dplyr::mutate(eCLIP = ifelse(Symbol %in% transC_lists$Gene, "RBM20_target", "no_target")) %>% 
  dplyr::mutate(eCLIP = ifelse(FDR < 0.05 & abs(log2FC) > log2(1.5), eCLIP, "insignificant"))
ggplot(o_map_WT_KO %>% dplyr::arrange(desc(eCLIP))) + geom_point(aes(x = log2FC, y = -log10(FDR), color = eCLIP, fill = eCLIP), size = 3, shape = 21) +
  scale_fill_manual(values = c("grey", "black", "magenta")) +
  scale_color_manual(values = c("grey", "black", "purple")) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") + geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
  ylab("-Log10 (p-value)") + xlab("Log2 (fold change O-MAP/Input)")
ggsave("/home/shared_folder/output_plots/RBM20_targets_HC_volcano_SE.pdf")

WT_For_Volcano <- read.csv("/home/shared_folder/pulldown_data/transcripts_data/WT_For_Volcano.csv") %>%
  dplyr::mutate(eCLIP = ifelse(GeneName %in% transC_lists$Gene, "RBM20_target", "no_target")) %>% 
  dplyr::mutate(eCLIP = ifelse(padj < 0.05 & abs(log2FoldChange) > log2(1.5), eCLIP, "insignificant"))
ggplot(WT_For_Volcano %>% dplyr::arrange(desc(eCLIP))) + geom_point(aes(x = log2FoldChange, y = -log10(padj), color = eCLIP, fill = eCLIP), size = 3, shape = 21) +
  scale_fill_manual(values = c("grey", "black", "magenta")) +
  scale_color_manual(values = c("grey", "black", "purple")) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") + geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
  ylab("-Log10 (p-value)") + xlab("Log2 (fold change O-MAP/Input)")
ggsave("/home/shared_folder/output_plots/RBM20_targets_HC_volcano_WT.pdf")

WT_For_Volcano <- WT_For_Volcano %>%
  dplyr::mutate(chr2 = ifelse(Chr == "2", "chr2", "other")) %>% 
  dplyr::mutate(chr = ifelse(padj < 0.05 & abs(log2FoldChange) > log2(1.5), chr2, "insignificant"))
ggplot(WT_For_Volcano %>% dplyr::arrange(desc(chr))) + geom_point(aes(x = log2FoldChange, y = -log10(padj), color = chr, fill = chr), size = 3, shape = 21) +
  scale_fill_manual(values = c("magenta", "grey", "black")) +
  scale_color_manual(values = c("purple", "grey", "black")) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") + geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
  ylab("-Log10 (p-value)") + xlab("Log2 (fold change O-MAP/Input)")
ggsave("/home/shared_folder/output_plots/chr2_volcano_WT.pdf")

WT_For_Volcano <- WT_For_Volcano %>%
  dplyr::mutate(TIDs = ifelse(EnsemblName %in% TIDs_genes, "TID", "no_target")) %>% 
  dplyr::mutate(TIDs = ifelse(padj < 0.05 & abs(log2FoldChange) > log2(1.5), TIDs, "insignificant"))
ggplot(WT_For_Volcano %>% dplyr::arrange(desc(TIDs))) + geom_point(aes(x = log2FoldChange, y = -log10(padj), color = TIDs, fill = TIDs), size = 3, shape = 21) +
  scale_fill_manual(values = c("grey", "black", "magenta")) +
  scale_color_manual(values = c("grey", "black", "purple")) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") + geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
  ylab("-Log10 (p-value)") + xlab("Log2 (fold change O-MAP/Input)")
ggsave("/home/shared_folder/output_plots/TIDs_volcano_WT.pdf")


KO_For_Volcano <- read.csv("/home/shared_folder/pulldown_data/transcripts_data/KO_For_Volcano.csv") %>%
  dplyr::mutate(eCLIP = ifelse(GeneName %in% transC_lists$Gene, "RBM20_target", "no_target")) %>% 
  dplyr::mutate(eCLIP = ifelse(padj < 0.05 & abs(log2FoldChange) > log2(1.5), eCLIP, "insignificant"))
ggplot(KO_For_Volcano %>% dplyr::arrange(desc(eCLIP))) + geom_point(aes(x = log2FoldChange, y = -log10(padj), color = eCLIP, fill = eCLIP), size = 3, shape = 21) +
  scale_fill_manual(values = c("grey", "black", "magenta")) +
  scale_color_manual(values = c("grey", "black", "purple")) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") + geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
  ylab("-Log10 (p-value)") + xlab("Log2 (fold change O-MAP/Input)")
ggsave("/home/shared_folder/output_plots/RBM20_targets_HC_volcano_KO.pdf")

KO_For_Volcano <- KO_For_Volcano %>%
  dplyr::mutate(chr2 = ifelse(Chr == "2", "chr2", "other")) %>% 
  dplyr::mutate(chr = ifelse(padj < 0.05 & abs(log2FoldChange) > log2(1.5), chr2, "insignificant"))
ggplot(KO_For_Volcano %>% dplyr::arrange(desc(chr))) + geom_point(aes(x = log2FoldChange, y = -log10(padj), color = chr, fill = chr), size = 3, shape = 21) +
  scale_fill_manual(values = c("magenta", "grey", "black")) +
  scale_color_manual(values = c("purple", "grey", "black")) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") + geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
  ylab("-Log10 (p-value)") + xlab("Log2 (fold change O-MAP/Input)")
ggsave("/home/shared_folder/output_plots/chr2_volcano_KO.pdf")

KO_For_Volcano <- KO_For_Volcano %>%
  dplyr::mutate(TIDs = ifelse(EnsemblName %in% TIDs_genes, "TID", "no_target")) %>% 
  dplyr::mutate(TIDs = ifelse(padj < 0.05 & abs(log2FoldChange) > log2(1.5), TIDs, "insignificant"))
ggplot(KO_For_Volcano %>% dplyr::arrange(desc(TIDs))) + geom_point(aes(x = log2FoldChange, y = -log10(padj), color = TIDs, fill = TIDs), size = 3, shape = 21) +
  scale_fill_manual(values = c("grey", "black", "magenta")) +
  scale_color_manual(values = c("grey", "black", "purple")) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") + geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
  ylab("-Log10 (p-value)") + xlab("Log2 (fold change O-MAP/Input)")
ggsave("/home/shared_folder/output_plots/TIDs_volcano_KO.pdf")

#venn diagrams - supplementary figure 6 panel B
RBM20_eCLIP_peaks <- read_excel("/home/shared_folder/chromatin_data/eCLIP_data/RBM20_eCLIP_peaks.xlsx")
genes = RBM20_eCLIP_peaks$`Gene Symbol(s)`
genes_f = c()
for(i in genes){
  if(length(strsplit(i, split = "\\|")[[1]]) == 1){
    genes_f = append(genes_f, i)
  }
}
genes_peaks = data.frame(geneSymbol = genes_f)
genes_peaks$n = 1
genes_peaks = genes_peaks %>% dplyr::group_by(geneSymbol) %>% dplyr::summarise(peaks = sum(n))

geneInfo <- read.table("/home/shared_folder/pulldown_data/transcripts_data/geneInfo.tab", quote="\"", comment.char="", skip = 1)
colnames(geneInfo) = c("Gene", "Symbol", "gene_type")

WT_OMAP <- read_excel("/home/shared_folder/pulldown_data/transcripts_data/o_map_WT.xlsx") %>% dplyr::left_join(geneInfo, by = "Gene")
KO_OMAP <- read_excel("/home/shared_folder/pulldown_data/transcripts_data/o_map_RBM20KO.xlsx") %>% dplyr::left_join(geneInfo, by = "Gene")

venn = venn.diagram(x = list(na.omit(WT_OMAP$Symbol), na.omit(KO_OMAP$Symbol), na.omit(genes_peaks$geneSymbol)),
                    category.names = c("WT\nO-MAP", "RBM20KO\nO-MAP", "RBM20 eCLIP targets"),
                    fill = c("#3283FE", "#F6222E", "grey"),
                    filename = NULL)
tab = geneInfo %>% dplyr::mutate(eCLIP = ifelse(Symbol %in% na.omit(genes_peaks$geneSymbol), "eCLIP", "z"),
                                 OMAP_WT = ifelse(Symbol %in% na.omit(WT_OMAP$Symbol), "WT_pd", "z"),
                                 OMAP_KO = ifelse(Symbol %in% na.omit(KO_OMAP$Symbol), "KO_pd", "z"))
test_WT = fisher.test(table(tab %>% dplyr::select(eCLIP, OMAP_WT)))
test_KO = fisher.test(table(tab %>% dplyr::select(eCLIP, OMAP_KO)))

pdf("/home/shared_folder/output_plots/VennDIagram_transcripts.pdf")
grid.draw(venn)
grid.text(paste0("p-val = ", signif(test_WT[["p.value"]], digits = 3)), x = 0.3, y = 0.4)
grid.text(paste0("p-val = ", signif(test_KO[["p.value"]], digits = 3)), x = 0.7, y = 0.4)
dev.off()

transC_lists <- read_excel("/home/shared_folder/pulldown_data/transcripts_data/transC_lists.xlsx") %>%
  dplyr::filter(HC_eCLIP == "Yes")

venn = venn.diagram(x = list(na.omit(WT_OMAP$Symbol), na.omit(KO_OMAP$Symbol), na.omit(transC_lists$Gene)),
                    category.names = c("WT\nO-MAP", "RBM20KO\nO-MAP", "RBM20 eCLIP HC targets"),
                    fill = c("#3283FE", "#F6222E", "grey"),
                    filename = NULL)
tab = geneInfo %>% dplyr::mutate(eCLIP = ifelse(Symbol %in% na.omit(transC_lists$Gene), "eCLIP", "z"),
                                 OMAP_WT = ifelse(Symbol %in% na.omit(WT_OMAP$Symbol), "WT_pd", "z"),
                                 OMAP_KO = ifelse(Symbol %in% na.omit(KO_OMAP$Symbol), "KO_pd", "z"))
test_WT = fisher.test(table(tab %>% dplyr::select(eCLIP, OMAP_WT)))
test_KO = fisher.test(table(tab %>% dplyr::select(eCLIP, OMAP_KO)))

pdf("/home/shared_folder/output_plots/VennDIagram_transcripts_HC_eCLIP.pdf")
grid.draw(venn)
grid.text(paste0("p-val = ", signif(test_WT[["p.value"]], digits = 3)), x = 0.15, y = 0.43)
grid.text(paste0("p-val = ", signif(test_KO[["p.value"]], digits = 3)), x = 0.85, y = 0.43)
dev.off()

genes_chr2 <- read.csv("/home/shared_folder/chromatin_data/Chrom_HMM_data/tss.csv") %>% dplyr::filter(seqnames == "2")
venn = venn.diagram(x = list(na.omit(WT_OMAP$Gene), na.omit(KO_OMAP$Gene), na.omit(genes_chr2$X)),
                    category.names = c("WT\nO-MAP", "RBM20KO\nO-MAP", "genes on chromosome 2"),
                    fill = c("#3283FE", "#F6222E", "grey"),
                    filename = NULL)
tab = geneInfo %>% dplyr::mutate(chr2 = ifelse(Gene %in% na.omit(genes_chr2$X), "chr2", "z"),
                                 OMAP_WT = ifelse(Symbol %in% na.omit(WT_OMAP$Symbol), "WT_pd", "z"),
                                 OMAP_KO = ifelse(Symbol %in% na.omit(KO_OMAP$Symbol), "KO_pd", "z"))
test_WT = fisher.test(table(tab %>% dplyr::select(chr2, OMAP_WT)))
test_KO = fisher.test(table(tab %>% dplyr::select(chr2, OMAP_KO)))

pdf("/home/shared_folder/output_plots/VennDIagram_transcripts_chr2.pdf")
grid.draw(venn)
grid.text(paste0("p-val = ", signif(test_WT[["p.value"]], digits = 3)), x = 0.12, y = 0.37)
grid.text(paste0("p-val = ", signif(test_KO[["p.value"]], digits = 3)), x = 0.88, y = 0.37)
dev.off()

act_genes <- read.csv("/home/shared_folder/chromatin_data/Chrom_HMM_data/act_genes.csv")
rep_genes <- read.csv("/home/shared_folder/chromatin_data/Chrom_HMM_data/rep_genes.csv")
biv_genes <- read.csv("/home/shared_folder/chromatin_data/Chrom_HMM_data/biv_genes.csv")
TIDs_genes = unique(c(act_genes$geneID, rep_genes$geneID, biv_genes$geneID))

venn = venn.diagram(x = list(na.omit(WT_OMAP$Gene), na.omit(KO_OMAP$Gene), na.omit(TIDs_genes)),
                    category.names = c("WT\nO-MAP", "RBM20KO\nO-MAP", "TIDs enclosed genes"),
                    fill = c("#3283FE", "#F6222E", "grey"),
                    filename = NULL)
tab = geneInfo %>% dplyr::mutate(tids = ifelse(Gene %in% na.omit(TIDs_genes), "tids", "z"),
                                 OMAP_WT = ifelse(Symbol %in% na.omit(WT_OMAP$Symbol), "WT_pd", "z"),
                                 OMAP_KO = ifelse(Symbol %in% na.omit(KO_OMAP$Symbol), "KO_pd", "z"))
test_WT = fisher.test(table(tab %>% dplyr::select(tids, OMAP_WT)))
test_KO = fisher.test(table(tab %>% dplyr::select(tids, OMAP_KO)))

pdf("/home/shared_folder/output_plots/VennDIagram_transcripts_TIDs.pdf")
grid.draw(venn)
grid.text(paste0("p-val = ", signif(test_WT[["p.value"]], digits = 3)), x = 0.12, y = 0.37)
grid.text(paste0("p-val = ", signif(test_KO[["p.value"]], digits = 3)), x = 0.88, y = 0.37)
dev.off()