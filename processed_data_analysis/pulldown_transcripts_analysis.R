library(dplyr)
library(ggplot2)
library(reshape2)
library(readxl)
library(VennDiagram)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg38)
library(transite)
library(biomaRt)
library(plyranges)
library(ggrepel)
library(ADImpute)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size=15, face="bold", margin = margin(10,0,10,0))))

##O-MAP pulldown transcripts
#volcano plots
transC_lists <- read_excel("/home/shared_folder/processed_data_analysis/chromatin_data/eCLIP_data/transC_lists.xlsx") %>%
  dplyr::filter(HC_eCLIP == "Yes")

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

#venn diagrams
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

#GC content analysis
Homo_sapiens.GRCh38.113 = read.delim("/home/shared_folder/genome_track_files/Homo_sapiens.GRCh38.113.gtf", header=FALSE, comment.char="#")
Homo_sapiens.GRCh38.113$geneID = sapply(Homo_sapiens.GRCh38.113$V9, function(x){strsplit(strsplit(x, "gene_id ")[[1]][2], ";")[[1]][1]})

#WT
WT_coords = Homo_sapiens.GRCh38.113 %>% dplyr::filter(geneID %in% WT_OMAP$Gene)
WT_coords = WT_coords %>% dplyr::filter(!is.na(as.numeric(V1)))
WT_coords$V1 = paste0("chr", WT_coords$V1)
WT_exons = WT_coords %>% dplyr::filter(V3 == "exon") %>% dplyr::select(c(V1, V4, V5, geneID)) %>% dplyr::distinct() %>% dplyr::arrange(V1) %>% dplyr::arrange(V4)
WT_pos = data.frame(V1 = c(), V4 = c(), V5 = c(), geneID = c())

for(i in unique(WT_exons$geneID)){
  WT_i = WT_exons %>% dplyr::filter(geneID == i)
  WT_pos = rbind(WT_pos, data.frame(V1 = c(WT_i$V1[1]), V4 = c(min(WT_i$V4)), V5 = c(max(WT_i$V5)), geneID = c(i)))
}

WT_5UTR = WT_coords %>% dplyr::filter(V3 == "five_prime_utr") %>% dplyr::select(c(V1, V4, V5, geneID)) %>% dplyr::distinct()
WT_3UTR = WT_coords %>% dplyr::filter(V3 == "three_prime_utr") %>% dplyr::select(c(V1, V4, V5, geneID)) %>% dplyr::distinct()

genome <- BSgenome.Hsapiens.UCSC.hg38
colnames(WT_exons) = c("chr", "start", "end", "geneID")
colnames(WT_pos) = c("chr", "start", "end", "geneID")
exons = makeGRangesFromDataFrame(WT_exons, keep.extra.columns = TRUE)
genes = makeGRangesFromDataFrame(WT_pos, keep.extra.columns = TRUE)
genes_list_WT <- split(genes, genes@elementMetadata@listData[["geneID"]])
exons_list_WT <- split(exons, exons@elementMetadata@listData[["geneID"]])

colnames(WT_5UTR) = c("chr", "start", "end", "geneID")
colnames(WT_3UTR) = c("chr", "start", "end", "geneID")
UTR5 = makeGRangesFromDataFrame(WT_5UTR, keep.extra.columns = TRUE)
UTR3 = makeGRangesFromDataFrame(WT_3UTR, keep.extra.columns = TRUE)
UTR5_list_WT <- split(UTR5, UTR5@elementMetadata@listData[["geneID"]])
UTR3_list_WT <- split(UTR3, UTR3@elementMetadata@listData[["geneID"]])

introns_list_WT = mapply(IRanges::setdiff, genes_list_WT, exons_list_WT)
exons_list_WT = mapply(IRanges::setdiff, exons_list_WT, UTR5_list_WT)
exons_list_WT = mapply(IRanges::setdiff, exons_list_WT, UTR3_list_WT)

GC_exons_WT = c()
GC_introns_WT = c()
exons_df_WT = data.frame(name = c(), seq = c())
introns_df_WT = data.frame(name = c(), seq = c())

for(i in 1:length(genes_list_WT)){
  seq <- getSeq(genome, exons_list_WT[[i]])
  if(length(seq) != 0){gc <- letterFrequency(seq, letters="GC", as.prob=TRUE)
  GC_exons_WT = append(GC_exons_WT, mean(gc))
  exons_df_WT = rbind(exons_df_WT, data.frame(name = genes_list_WT[[i]]$geneID, seq = seq))}
  
  seq <- getSeq(genome, introns_list_WT[[i]])
  if(length(seq) != 0){gc <- letterFrequency(seq, letters="GC", as.prob=TRUE)
  GC_introns_WT = append(GC_introns_WT, mean(gc))
  introns_df_WT = rbind(introns_df_WT, data.frame(name = genes_list_WT[[i]]$geneID, seq = seq))}
  
}

GC_5UTR_WT = c()
UTR5_df_WT = data.frame(name = c(), seq = c())
for(i in 1:length(UTR5_list_WT)){
  seq <- getSeq(genome, UTR5_list_WT[[i]])
  if(length(seq) != 0){gc <- letterFrequency(seq, letters="GC", as.prob=TRUE)
  GC_5UTR_WT = append(GC_5UTR_WT, mean(gc))
  UTR5_df_WT = rbind(UTR5_df_WT, data.frame(name = UTR5_list_WT[[i]]$geneID, seq = seq))}
}

GC_3UTR_WT = c()
UTR3_df_WT = data.frame(name = c(), seq = c())
for(i in 1:length(UTR3_list_WT)){
  seq <- getSeq(genome, UTR3_list_WT[[i]])
  if(length(seq) != 0){gc <- letterFrequency(seq, letters="GC", as.prob=TRUE)
  GC_3UTR_WT = append(GC_3UTR_WT, mean(gc))
  UTR3_df_WT = rbind(UTR3_df_WT, data.frame(name = UTR3_list_WT[[i]]$geneID, seq = seq))}
}

summary_WT = rbind(data.frame(GC = GC_3UTR_WT, type = rep("3_UTR", n = length(GC_3UTR_WT))),
                   data.frame(GC = GC_5UTR_WT, type = rep("5_UTR", n = length(GC_5UTR_WT))),
                   data.frame(GC = GC_introns_WT, type = rep("introns", n = length(GC_introns_WT))),
                   data.frame(GC = GC_exons_WT, type = rep("exons", n = length(GC_exons_WT))))

#KO
KO_coords = Homo_sapiens.GRCh38.113 %>% dplyr::filter(geneID %in% KO_OMAP$Gene)
KO_coords = KO_coords %>% dplyr::filter(!is.na(as.numeric(V1)))
KO_coords$V1 = paste0("chr", KO_coords$V1)
KO_exons = KO_coords %>% dplyr::filter(V3 == "exon") %>% dplyr::select(c(V1, V4, V5, geneID)) %>% dplyr::distinct() %>% dplyr::arrange(V1) %>% dplyr::arrange(V4)
KO_pos = data.frame(V1 = c(), V4 = c(), V5 = c(), geneID = c())

for(i in unique(KO_exons$geneID)){
  KO_i = KO_exons %>% dplyr::filter(geneID == i)
  KO_pos = rbind(KO_pos, data.frame(V1 = c(KO_i$V1[1]), V4 = c(min(KO_i$V4)), V5 = c(max(KO_i$V5)), geneID = c(i)))
}

KO_5UTR = KO_coords %>% dplyr::filter(V3 == "five_prime_utr") %>% dplyr::select(c(V1, V4, V5, geneID)) %>% dplyr::distinct()
KO_3UTR = KO_coords %>% dplyr::filter(V3 == "three_prime_utr") %>% dplyr::select(c(V1, V4, V5, geneID)) %>% dplyr::distinct()

colnames(KO_exons) = c("chr", "start", "end", "geneID")
colnames(KO_pos) = c("chr", "start", "end", "geneID")
exons = makeGRangesFromDataFrame(KO_exons, keep.extra.columns = TRUE)
genes = makeGRangesFromDataFrame(KO_pos, keep.extra.columns = TRUE)
genes_list_KO <- split(genes, genes@elementMetadata@listData[["geneID"]])
exons_list_KO <- split(exons, exons@elementMetadata@listData[["geneID"]])

colnames(KO_5UTR) = c("chr", "start", "end", "geneID")
colnames(KO_3UTR) = c("chr", "start", "end", "geneID")
UTR5 = makeGRangesFromDataFrame(KO_5UTR, keep.extra.columns = TRUE)
UTR3 = makeGRangesFromDataFrame(KO_3UTR, keep.extra.columns = TRUE)
UTR5_list_KO <- split(UTR5, UTR5@elementMetadata@listData[["geneID"]])
UTR3_list_KO <- split(UTR3, UTR3@elementMetadata@listData[["geneID"]])

introns_list_KO = mapply(IRanges::setdiff, genes_list_KO, exons_list_KO)
exons_list_KO = mapply(IRanges::setdiff, exons_list_KO, UTR5_list_KO)
exons_list_KO = mapply(IRanges::setdiff, exons_list_KO, UTR3_list_KO)

GC_exons_KO = c()
GC_introns_KO = c()
exons_df_KO = data.frame(name = c(), seq = c())
introns_df_KO = data.frame(name = c(), seq = c())
for(i in 1:length(genes_list_KO)){
  seq <- getSeq(genome, exons_list_KO[[i]])
  if(length(seq) != 0){gc <- letterFrequency(seq, letters="GC", as.prob=TRUE)
  GC_exons_KO = append(GC_exons_KO, mean(gc))
  exons_df_KO = rbind(exons_df_KO, data.frame(name = genes_list_KO[[i]]$geneID, seq = seq))}
  
  seq <- getSeq(genome, introns_list_KO[[i]])
  if(length(seq) != 0){gc <- letterFrequency(seq, letters="GC", as.prob=TRUE)
  GC_introns_KO = append(GC_introns_KO, mean(gc))
  introns_df_KO = rbind(introns_df_KO, data.frame(name = genes_list_KO[[i]]$geneID, seq = seq))}
}

GC_5UTR_KO = c()
UTR5_df_KO = data.frame(name = c(), seq = c())
for(i in 1:length(UTR5_list_KO)){
  seq <- getSeq(genome, UTR5_list_KO[[i]])
  if(length(seq) != 0){gc <- letterFrequency(seq, letters="GC", as.prob=TRUE)
  GC_5UTR_KO = append(GC_5UTR_KO, mean(gc))
  UTR5_df_KO = rbind(UTR5_df_KO, data.frame(name = UTR5_list_KO[[i]]$geneID, seq = seq))}
}

GC_3UTR_KO = c()
UTR3_df_KO = data.frame(name = c(), seq = c())
for(i in 1:length(UTR3_list_KO)){
  seq <- getSeq(genome, UTR3_list_KO[[i]])
  if(length(seq) != 0){gc <- letterFrequency(seq, letters="GC", as.prob=TRUE)
  GC_3UTR_KO = append(GC_3UTR_KO, mean(gc))
  UTR3_df_KO = rbind(UTR3_df_KO, data.frame(name = UTR3_list_KO[[i]]$geneID, seq = seq))}
}

summary_KO = rbind(data.frame(GC = GC_3UTR_KO, type = rep("3_UTR", n = length(GC_3UTR_KO))),
                   data.frame(GC = GC_5UTR_KO, type = rep("5_UTR", n = length(GC_5UTR_KO))),
                   data.frame(GC = GC_introns_KO, type = rep("introns", n = length(GC_introns_KO))),
                   data.frame(GC = GC_exons_KO, type = rep("exons", n = length(GC_exons_KO))))

summary = rbind(summary_WT %>% dplyr::mutate(cell = "WT"), summary_KO %>% dplyr::mutate(cell = "KO"))

ggplot(summary, aes(x = cell, y = GC)) + geom_boxplot(aes(fill = cell)) +
  facet_grid(~type) + theme(strip.background = element_blank()) +
  scale_fill_manual(values = c("red", "#3283FE"))
ggsave("/home/shared_folder/output_plots/boxplot_WT_KO.pdf")

int_p = ks.test(GC_introns_KO, GC_introns_WT)["p.value"][[1]]
ggplot(summary %>% dplyr::filter(type == "introns")) + stat_ecdf(aes(x = GC, color = cell)) + scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("red", "#3283FE")) + ggtitle(paste0("p-val introns = ", signif(int_p, 4))) + geom_vline(xintercept = 0.5, linetype = "dashed")
ggsave("/home/shared_folder/output_plots/GC_introns.pdf", width = 5, height = 5)

ex_p = ks.test(GC_exons_KO, GC_exons_WT)["p.value"][[1]]
ggplot(summary %>% dplyr::filter(type == "exons")) + stat_ecdf(aes(x = GC, color = cell)) + scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("red", "#3283FE")) + ggtitle(paste0("p-val exons = ", signif(ex_p, 4))) + geom_vline(xintercept = 0.5, linetype = "dashed")
ggsave("/home/shared_folder/output_plots/GC_exons.pdf", width = 5, height = 5)

u3_p = ks.test(GC_3UTR_KO, GC_3UTR_WT)["p.value"][[1]]
ggplot(summary %>% dplyr::filter(type == "3_UTR")) + stat_ecdf(aes(x = GC, color = cell)) + scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("red", "#3283FE")) + ggtitle(paste0("p-val 3'UTR = ", signif(u3_p, 4))) + geom_vline(xintercept = 0.5, linetype = "dashed")
ggsave("/home/shared_folder/output_plots/GC_3UTR.pdf", width = 5, height = 5)

u5_p = ks.test(GC_5UTR_KO, GC_5UTR_WT)["p.value"][[1]]
ggplot(summary %>% dplyr::filter(type == "5_UTR")) + stat_ecdf(aes(x = GC, color = cell)) + scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("red", "#3283FE")) + ggtitle(paste0("p-val 5'UTR = ", signif(u5_p, 4))) + geom_vline(xintercept = 0.5, linetype = "dashed")
ggsave("/home/shared_folder/output_plots/GC_5UTR.pdf", width = 5, height = 5)

#motifs
motifs = get_motifs_meta_info()
m_QKI_ELAVL1 = get_motif_by_id(c("M046_0.6", "149_16041388", "M031_0.6", "M108_0.6", "M112_0.6", "M127_0.6", "M232_0.6"))

bases = c("A", "G", "C", "T")
check_bases <- function(str){
  val = strsplit(str, "")[[1]]
  sum = 0
  for(i in val){
    if(!i %in% bases){
      sum = sum + 1
    }
  }
  return(sum)
}

background_WT = data.frame(name = c(), seq = c())
for(i in 1:length(genes_list_WT)){
  seq <- getSeq(genome, genes_list_WT[[i]])
  background_WT = rbind(background_WT, data.frame(name = genes_list_WT[[i]]$geneID, seq = seq))
}
background_WT$keep = as.character(lapply(background_WT$seq, check_bases))
background_WT = background_WT %>% dplyr::filter(keep == 0)
back_WT = background_WT$seq
back_WT = chartr("T", "U", back_WT)
names(back_WT) = background_WT$name
back_WT = list(back_WT)

exons_df_WT$keep = as.character(lapply(exons_df_WT$seq, check_bases))
exons_df_WT = exons_df_WT %>% dplyr::filter(keep == 0)
exons_df_WT$nchar = sapply(exons_df_WT$seq, nchar)
exons_df_WT = exons_df_WT[exons_df_WT$nchar > 7,]
exons_WT = exons_df_WT$seq
exons_WT = chartr("T", "U", exons_WT)
exons_df_WT$name = paste0(exons_df_WT$name, "|exon", seq(1, nrow(exons_df_WT)))
names(exons_WT) = exons_df_WT$name
exons_WT = list(exons_WT)

introns_df_WT$keep = as.character(lapply(introns_df_WT$seq, check_bases))
introns_df_WT = introns_df_WT %>% dplyr::filter(keep == 0)
introns_df_WT$nchar = sapply(introns_df_WT$seq, nchar)
introns_df_WT = introns_df_WT[introns_df_WT$nchar > 7,]
introns_WT = introns_df_WT$seq
introns_WT = chartr("T", "U", introns_WT)
introns_df_WT$name = paste0(introns_df_WT$name, "|intron", seq(1, nrow(introns_df_WT)))
names(introns_WT) = introns_df_WT$name
introns_WT = list(introns_WT)

UTR3_df_WT$keep = as.character(lapply(UTR3_df_WT$seq, check_bases))
UTR3_df_WT = UTR3_df_WT %>% dplyr::filter(keep == 0)
UTR3_df_WT$nchar = sapply(UTR3_df_WT$seq, nchar)
UTR3_df_WT = UTR3_df_WT[UTR3_df_WT$nchar > 7,]
UTR3_WT = UTR3_df_WT$seq
UTR3_WT = chartr("T", "U", UTR3_WT)
UTR3_df_WT$name = paste0(UTR3_df_WT$name, "|3UTR", seq(1, nrow(UTR3_df_WT)))
names(UTR3_WT) = UTR3_df_WT$name
UTR3_WT = list(UTR3_WT)

UTR5_df_WT$keep = as.character(lapply(UTR5_df_WT$seq, check_bases))
UTR5_df_WT = UTR5_df_WT %>% dplyr::filter(keep == 0)
UTR5_df_WT$nchar = sapply(UTR5_df_WT$seq, nchar)
UTR5_df_WT = UTR5_df_WT[UTR5_df_WT$nchar > 7,]
UTR5_WT = UTR5_df_WT$seq
UTR5_WT = chartr("T", "U", UTR5_WT)
UTR5_df_WT$name = paste0(UTR5_df_WT$name, "|5UTR", seq(1, nrow(UTR5_df_WT)))
names(UTR5_WT) = UTR5_df_WT$name
UTR5_WT = list(UTR5_WT)

back_WT = c(exons_df_WT$seq, introns_df_WT$seq, UTR3_df_WT$seq, UTR5_df_WT$seq)
names(back_WT) = c(exons_df_WT$name, introns_df_WT$name, UTR3_df_WT$name, UTR5_df_WT$name)
back_WT = as.list(back_WT)

UTR5_df_KO$nchar = sapply(UTR5_df_KO$seq, nchar)
UTR5_df_KO = UTR5_df_KO[UTR5_df_KO$nchar > 7,]
UTR5_KO = UTR5_df_KO$seq
UTR5_KO = chartr("T", "U", UTR5_KO)
names(UTR5_KO) = UTR5_df_KO$name
UTR5_KO = list(UTR5_KO)

enrichment_exons_WT = run_matrix_tsma(
  exons_WT,
  back_WT,
  motifs = m_QKI_ELAVL1,
  max_hits = Inf,
  threshold_method = "p_value",
  threshold_value = 0.05,
  cache = FALSE)

enrichment_introns_WT = run_matrix_tsma(
  introns_WT,
  back_WT,
  motifs = m_QKI_ELAVL1,
  max_hits = Inf,
  threshold_method = "p_value",
  threshold_value = 0.05,
  cache = FALSE)

enrichment_3UTR_WT = run_matrix_tsma(
  UTR3_WT,
  back_WT,
  motifs = m_QKI_ELAVL1,
  max_hits = Inf,
  threshold_method = "p_value",
  threshold_value = 0.05,
  cache = FALSE)

enrichment_5UTR_WT = run_matrix_tsma(
  UTR5_WT,
  back_WT,
  motifs = m_QKI_ELAVL1,
  max_hits = Inf,
  threshold_method = "p_value",
  threshold_value = 0.05,
  cache = FALSE)

#plots
tab = enrichment_introns_WT[["enrichment_dfs"]][[1]]
tab$region = "intron"
tab$mean_enrichment = mean(tab$enrichment)
tab = rbind(tab, enrichment_exons_WT[["enrichment_dfs"]][[1]] %>% dplyr::mutate(region = "exon", mean_enrichment = mean(enrichment)),
            enrichment_5UTR_WT[["enrichment_dfs"]][[1]] %>% dplyr::mutate(region = "5UTR", mean_enrichment = mean(enrichment)),
            enrichment_3UTR_WT[["enrichment_dfs"]][[1]] %>% dplyr::mutate(region = "3UTR", mean_enrichment = mean(enrichment)))
tab = tab %>% dplyr::mutate(sign = ifelse(adj_p_value < 0.05, "*", ""))
ggplot(tab) + geom_violin(aes(x = region, y = enrichment), fill = "lightgrey") +
  geom_jitter(aes(x = region, y = enrichment, color = motif_rbps, fill = motif_rbps), size = 2, shape = 21) +
  geom_hline(yintercept = 1, linetype = "dashed") + geom_text(aes(x = region, y = mean_enrichment, label = sign)) +
  scale_color_manual(values = c("#009E73", "#AA0DFE")) + scale_fill_manual(values = c("lightgreen", "#DEA0FD"))
ggsave("/home/shared_folder/output_plots/motifs_enrichment_WT_violin.pdf")


