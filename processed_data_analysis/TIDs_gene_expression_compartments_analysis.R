library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(ensembldb)
library(topGO)
library(limma)
library(org.Hs.eg.db)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(gtools)
library(VennDiagram)
library(gridExtra)
library(readxl)
library(GO.db)
library(AnnotationDbi)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0))))

##TIDs/CIDs compartments analysis
#TIDs compartment changes
TIDs = read.delim("/home/shared_folder/pulldown_data/chromatin_data/GSE234808_WT_Interacting_Domains.bed", header=FALSE)
colnames(TIDs) = c("chr", "start", "end", "id")
trans_TIDs = TIDs %>% dplyr::filter(chr != "chr2")

PC_WT = read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_WT_100000.bedGraph", header=FALSE, skip = 4)
colnames(PC_WT) = c("chr", "start", "end", "WT")
PC_RBM20KO = read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_KO_100000.bedGraph", header=FALSE, skip = 4)
colnames(PC_RBM20KO) = c("chr", "start", "end", "RBM20KO")

PC_combined = PC_WT %>% dplyr::left_join(PC_RBM20KO, by = c("chr", "start", "end")) %>%
  dplyr::filter(chr != "chr2")

tids_PC_comb =  data.frame(chr = c(), start = c(), end = c(), WT = c(), RBM20KO = c(), id = c())
for(i in 1:length(rownames(trans_TIDs))){
  pc_i = PC_combined %>% dplyr::filter(chr == trans_TIDs$chr[i] & start > trans_TIDs$start[i] & end < trans_TIDs$end[i])
  pc_i$id = as.character(rep(c(trans_TIDs$id[i]), n = length(rownames(pc_i))))
  tids_PC_comb = rbind(tids_PC_comb, pc_i)
}
PC_combined = PC_combined %>% dplyr::filter(! chr %in% c("chrX", "chr21", "chr22"))
tids_PC_comb_f = tids_PC_comb %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>%
  dplyr::filter(! chr %in% c("chrX", "chr21", "chr22"))
wt_all = wilcox.test(PC_combined$WT, PC_combined$RBM20KO, paired = TRUE, alternative = "two.sided") #0.147
wt = wilcox.test(tids_PC_comb_f$WT, tids_PC_comb_f$RBM20KO, paired = TRUE, alternative = "two.sided") #2.402e-05

PC_combined = PC_combined %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>%
  dplyr::mutate(TIDs = ifelse(bin %in% tids_PC_comb_f$bin, "TIDs", ""))
ggplot(PC_combined %>% dplyr::arrange(TIDs)) + geom_point(aes(x = WT, y = RBM20KO, fill = TIDs, color = TIDs), shape = 21, size = 2) +
  scale_fill_manual(values = c("white", "grey")) + scale_color_manual(values = c("grey", "black")) +
  scale_y_continuous(limits = c(-3, 3)) + scale_x_continuous(limits = c(-3, 3)) +
  ggtitle(paste0("WT vs RBM20KO PC values \np-value TIDs = ", signif(wt[["p.value"]], digits = 3))) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") + geom_hline(aes(yintercept = 0), linetype = "dotted") +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/compartment_changes_all_TIDs.pdf")

#CIDs compartment changes
CIDs = TIDs %>% dplyr::filter(chr == "chr2")
PC_combined_cis = PC_WT %>% dplyr::left_join(PC_RBM20KO, by = c("chr", "start", "end")) %>%
  dplyr::filter(chr == "chr2")

cids_PC_comb =  data.frame(chr = c(), start = c(), end = c(), WT = c(), RBM20KO = c(), id = c())
for(i in 1:length(rownames(CIDs))){
  pc_i = PC_combined_cis  %>% dplyr::filter(chr == CIDs$chr[i] & start > CIDs$start[i] & end < CIDs$end[i])
  pc_i$id = as.character(rep(c(CIDs$id[i]), n = length(rownames(pc_i))))
  cids_PC_comb = rbind(cids_PC_comb, pc_i)
}
wt_all = wilcox.test(PC_combined_cis$WT, PC_combined_cis$RBM20KO, paired = TRUE, alternative = "two.sided") #0.4908
wt = wilcox.test(cids_PC_comb$WT, cids_PC_comb$RBM20KO, paired = TRUE, alternative = "two.sided") #0.1433

cids_PC_comb = cids_PC_comb %>% dplyr::mutate(bin = paste0(chr, "_", start))
PC_combined_cis = PC_combined_cis %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>%
  dplyr::mutate(CIDs = ifelse(bin %in% cids_PC_comb$bin, "CIDs", ""))
ggplot(PC_combined_cis %>% dplyr::arrange(CIDs)) + geom_point(aes(x = WT, y = RBM20KO, fill = CIDs, color = CIDs), shape = 21, size = 2) +
  scale_fill_manual(values = c("white", "grey")) + scale_color_manual(values = c("grey", "black")) +
  scale_y_continuous(limits = c(-3, 3)) + scale_x_continuous(limits = c(-3, 3)) +
  ggtitle(paste0("WT vs RBM20KO PC values \np-value CIDs = ", signif(wt[["p.value"]], digits = 3))) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") + geom_hline(aes(yintercept = 0), linetype = "dotted") +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/compartment_changes_all_CIDs.pdf")

##ChromHMM data integration
tss <- read.csv("/home/shared_folder/chromatin_data/ChromHMM_data/tss.csv")
colnames(tss) = c("geneID", "chr", "start", "end", "width", "strand", "description")
tss = tss %>% dplyr::select(c("geneID", "chr", "start", "end", "width", "strand"))
tss = tss %>% dplyr::mutate(chr = paste0("chr", chr))
geneInfo <- read.table("/home/shared_folder/gene_expression_data/DGE_RUES2_data/geneInfo.tab", quote="\"", comment.char="", skip = 1)
colnames(geneInfo) = c("geneID", "gene_name", "gene_type")

PC_combined = PC_WT %>% dplyr::left_join(PC_RBM20KO, by = c("chr", "start", "end"))
ChromHMM_data <- read.delim("/home/shared_folder/chromatin_data/ChromHMM_data/DFF334WAB_RUES2ChromHMM1.bed", header=FALSE) %>% dplyr::select(c("V1", "V2", "V3", "V4"))
colnames(ChromHMM_data) = c("chr", "start", "end", "signature")
ChromHMM_data = ChromHMM_data %>% dplyr::mutate(length = end - start)
PC_combined = PC_combined %>% dplyr::mutate(bin = paste0(chr, "_", start))
chrom_sign = data.frame(EnhA1 = c(), EnhA2 = c(), EnhBiv = c(), EnhG1 = c(), EnhG2 = c(), EnhWk = c(), Het = c(), Quies = c(), ReprPC = c(), ReprPCWk = c(), TssA = c(), TssBiv = c(), TssFlnk = c(), TssFlnkD = c(), TssFlnkU = c(), Tx = c(), TxWk = c(), ZNFRpts = c(), unannotated = c(), bin = c()) 
ChromHMM_data = ChromHMM_data %>% dplyr::mutate(signature = ifelse(signature == "ZNF/Rpts", "ZNFRpts", signature))
for(i in 1:length(rownames(PC_combined))){
  Chrom_i = ChromHMM_data %>% dplyr::filter(chr == PC_combined$chr[i])
  Chrom_i = Chrom_i %>% dplyr::filter((start <= PC_combined$start[i] & end >= PC_combined$end[i]) |
                                        (start >= PC_combined$start[i] & end <= PC_combined$end[i]) |
                                        (start >= PC_combined$start[i] & start <= PC_combined$end[i]) |
                                        (end >= PC_combined$start[i] & end <= PC_combined$end[i]))
  if(length(rownames(Chrom_i)) != 0){
    Chrom_i = Chrom_i %>% dplyr::arrange(start)
    Chrom_i = Chrom_i %>% dplyr::mutate(start = ifelse(start < PC_combined$start[i], PC_combined$start[i], start))
    Chrom_i = Chrom_i %>% dplyr::mutate(end = ifelse(end > PC_combined$end[i], PC_combined$end[i], end))
    Chrom_i = Chrom_i %>% dplyr::mutate(length = end - start)
    Chrom_i = Chrom_i %>% dplyr::mutate(perc = length/100000 * 100)
    summ_i = Chrom_i %>% dplyr::group_by(signature) %>% dplyr::summarise(perc = sum(perc))
    summ_i = as.data.frame(summ_i)
  } else{summ_i = data.frame(signature = c(), perc = c())}
  
  for(j in unique(ChromHMM_data$signature)){
    if(! j %in% summ_i$signature){
      j_perc = data.frame(signature = j, perc = 0)
      summ_i = rbind(summ_i, j_perc)
    }
  }
  summ_i = summ_i %>% dplyr::arrange(signature)
  rownames(summ_i) = summ_i$signature
  summ_i = summ_i %>% dplyr::select(perc)
  tot = sum(summ_i$perc)
  summ_t = as.data.frame(t(summ_i))
  summ_t$unannotated = 100-tot
  summ_t$bin = PC_combined$bin[i]
  chrom_sign = rbind(chrom_sign, summ_t)
}
tids_PC_comb = tids_PC_comb %>% dplyr::mutate(bin = paste0(chr, "_", start))
tids_info = tids_PC_comb %>% dplyr::select(bin, id, WT, RBM20KO)
tids_sign = tids_PC_comb %>% dplyr::select(bin, id) %>% dplyr::left_join(chrom_sign, by = "bin")
melt_tids = melt(tids_sign) %>% dplyr::left_join(tids_info)
melt_tids = melt_tids %>% dplyr::mutate(delta = WT - RBM20KO)

chrom_tot = chrom_sign %>% dplyr::select(-bin)
chrom_tot = as.data.frame(colSums(chrom_tot)/length(rownames(chrom_tot)))
chrom_tot$variable = rownames(chrom_tot)
colnames(chrom_tot) = c("average", "variable")

melt_tids = melt_tids %>% dplyr::left_join(chrom_tot, by = "variable")
melt_tids = melt_tids %>% dplyr::mutate(value_norm = value/average)

melt_rep = melt_tids %>% dplyr::filter(variable %in% c("Het", "Quies", "ReprPCWk", "ReprPC")) %>% dplyr::group_by(bin) %>% dplyr::summarise(rep = mean(value_norm))
melt_act = melt_tids %>% dplyr::filter(variable %in% c("EnhA1", "EnhA2", "Tx", "TssA")) %>% dplyr::group_by(bin) %>% dplyr::summarise(act = mean(value_norm))
melt_biv = melt_tids %>% dplyr::filter(variable %in% c("EnhBiv", "TssBiv")) %>% dplyr::group_by(bin) %>% dplyr::summarise(biv = mean(value_norm))

tids_annot = tids_PC_comb %>% dplyr::left_join(melt_rep, by = "bin") %>% dplyr::left_join(melt_act, by = "bin") %>% dplyr::left_join(melt_biv, by = "bin")
tids_annot = tids_annot %>% dplyr::mutate(sign = ifelse(rep > act & biv < 1, "rep",
                                                        ifelse(act > rep & biv < 1, "act",
                                                               ifelse(biv >= 1, "biv", ""))))

TIDs_genes = data.frame(state = c(), dens = c(), n_genes = c(), n_bins = c())
act_tids = tids_annot %>% dplyr::filter(sign == "act")
act_tids_genes = data.frame(geneID = c())
for(i in 1:length(rownames(act_tids))){
  tss_i = tss %>% dplyr::filter(chr == act_tids$chr[i])
  tss_i = tss_i %>% dplyr::mutate(tss = ifelse(strand == "+", start, end)) %>% dplyr::filter(tss > act_tids$start[i] & tss < act_tids$end[i]) %>% dplyr::select(c(geneID, width))
  act_tids_genes = rbind(act_tids_genes, tss_i)
}
TIDs_genes = rbind(TIDs_genes, data.frame(state = "act", dens = length(rownames(act_tids_genes))/length(rownames(act_tids)), n_genes = length(rownames(act_tids_genes)), n_bins = length(rownames(act_tids))))
write.csv(act_tids_genes %>% dplyr::select(geneID), "/home/shared_folder/chromatin_data/ChromHMM_data/act_genes.csv")

rep_tids = tids_annot %>% dplyr::filter(sign == "rep")
rep_tids_genes = data.frame(geneID = c())
for(i in 1:length(rownames(rep_tids))){
  tss_i = tss %>% dplyr::filter(chr == rep_tids$chr[i])
  tss_i = tss_i %>% dplyr::mutate(tss = ifelse(strand == "+", start, end)) %>% dplyr::filter(tss > rep_tids$start[i] & tss < rep_tids$end[i]) %>% dplyr::select(c(geneID, width))
  rep_tids_genes = rbind(rep_tids_genes, tss_i)
}
TIDs_genes = rbind(TIDs_genes, data.frame(state = "rep", dens = length(rownames(rep_tids_genes))/length(rownames(rep_tids)), n_genes = length(rownames(rep_tids_genes)), n_bins = length(rownames(rep_tids))))
write.csv(rep_tids_genes %>% dplyr::select(geneID), "/home/shared_folder/chromatin_data/ChromHMM_data/rep_genes.csv")

biv_tids = tids_annot %>% dplyr::filter(sign == "biv")
biv_tids_genes = data.frame(geneID = c())
for(i in 1:length(rownames(biv_tids))){
  tss_i = tss %>% dplyr::filter(chr == biv_tids$chr[i])
  tss_i = tss_i %>% dplyr::mutate(tss = ifelse(strand == "+", start, end)) %>% dplyr::filter(tss > biv_tids$start[i] & tss < biv_tids$end[i]) %>% dplyr::select(c(geneID, width))
  biv_tids_genes = rbind(biv_tids_genes, tss_i)
}
dens = length(rownames(biv_tids_genes))/length(rownames(biv_tids))
TIDs_genes = rbind(TIDs_genes, data.frame(state = "biv", dens = length(rownames(biv_tids_genes))/length(rownames(biv_tids)), n_genes = length(rownames(biv_tids_genes)), n_bins = length(rownames(biv_tids))))
write.csv(biv_tids_genes %>% dplyr::select(geneID), "/home/shared_folder/chromatin_data/ChromHMM_data/biv_genes.csv")
write.csv(TIDs_genes, "/home/shared_folder/chromatin_data/ChromHMM_data/gene_denstiy.csv")

#TIDs scatterplots
PC_combined = PC_combined %>% dplyr::filter(! chr %in% c("chrX", "chr21", "chr22")) 

tids_PC_comb_a = tids_PC_comb %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>%
  dplyr::filter(! chr %in% c("chrX", "chr21", "chr22")) %>% dplyr::filter(bin %in% act_tids$bin)
wt = wilcox.test(tids_PC_comb_a$WT, tids_PC_comb_a$RBM20KO, paired = TRUE, alternative = "two.sided") #0.9278
PC_combined = PC_combined %>%
  dplyr::mutate(active = ifelse(bin %in% act_tids$bin, "active_TIDs", ""))
ggplot(PC_combined %>% dplyr::arrange(active)) + geom_point(aes(x = WT, y = RBM20KO, fill = active, color = active), shape = 21, size = 2) +
  scale_fill_manual(values = c("white", "#FEAF16")) + scale_color_manual(values = c("grey", "black")) +
  scale_y_continuous(limits = c(-3, 3)) + scale_x_continuous(limits = c(-3, 3)) +
  ggtitle(paste0("WT vs RBM20KO PC values \np-value active TIDs = ", signif(wt[["p.value"]], digits = 3))) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") + geom_hline(aes(yintercept = 0), linetype = "dotted") +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/compartment_changes_active_TIDs.pdf")

tids_PC_comb_r = tids_PC_comb %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>%
  dplyr::filter(! chr %in% c("chrX", "chr21", "chr22")) %>% dplyr::filter(bin %in% rep_tids$bin)
wt = wilcox.test(tids_PC_comb_r$WT, tids_PC_comb_r$RBM20KO, paired = TRUE, alternative = "two.sided") #1.553e-05
PC_combined = PC_combined %>%
  dplyr::mutate(repressive = ifelse(bin %in% rep_tids$bin, "repressive_TIDs", ""))
ggplot(PC_combined %>% dplyr::arrange(repressive)) + geom_point(aes(x = WT, y = RBM20KO, fill = repressive, color = repressive), shape = 21, size = 2) +
  scale_fill_manual(values = c("white", "#56B4E9")) + scale_color_manual(values = c("grey", "black")) +
  scale_y_continuous(limits = c(-3, 3)) + scale_x_continuous(limits = c(-3, 3)) +
  ggtitle(paste0("WT vs RBM20KO PC values \np-value repressive TIDs = ", signif(wt[["p.value"]], digits = 3))) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") + geom_hline(aes(yintercept = 0), linetype = "dotted") +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/compartment_changes_repressive_TIDs.pdf")

tids_PC_comb_b = tids_PC_comb %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>%
  dplyr::filter(! chr %in% c("chrX", "chr21", "chr22")) %>% dplyr::filter(bin %in% biv_tids$bin)
wt = wilcox.test(tids_PC_comb_b$WT, tids_PC_comb_b$RBM20KO, paired = TRUE, alternative = "two.sided") #0.00135
PC_combined = PC_combined %>%
  dplyr::mutate(bivalent = ifelse(bin %in% biv_tids$bin, "bivalent_TIDs", ""))
ggplot(PC_combined %>% dplyr::arrange(bivalent)) + geom_point(aes(x = WT, y = RBM20KO, fill = bivalent, color = bivalent), shape = 21, size = 2) +
  scale_fill_manual(values = c("white", "#1CBE4F")) + scale_color_manual(values = c("grey", "black")) +
  scale_y_continuous(limits = c(-3, 3)) + scale_x_continuous(limits = c(-3, 3)) +
  ggtitle(paste0("WT vs RBM20KO PC values \np-value bivalent TIDs = ", signif(wt[["p.value"]], digits = 3))) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") + geom_hline(aes(yintercept = 0), linetype = "dotted") +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/compartment_changes_bivalent_TIDs.pdf")

#TIDs statistics
TIDs_genes = TIDs_genes %>% dplyr::select(-dens)
TIDs_genes = melt(TIDs_genes)
ggplot(TIDs_genes) + geom_bar(aes(x = variable, y = value, fill = state), stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#FEAF16", "#1CBE4F", "#56B4E9")) + xlab("") + ylab("relative proportion")
ggsave("/home/shared_folder/output_plots/stats_TIDs_ChromHMM.pdf")

#gene length distribution
tids_genes = rbind(act_tids_genes %>% dplyr::mutate(state = "active"), rep_tids_genes %>%
                     dplyr::mutate(state = "repressive"), biv_tids_genes %>%
                     dplyr::mutate(state = "bivalent"))
ggplot(tids_genes) + stat_ecdf(aes(x = log2(width), color = state), size = 1) +
  scale_color_manual(values = c("#FEAF16", "#1CBE4F", "#56B4E9")) +
  ylab("cumulative frequency") + xlab("log2(gene length)")
ggsave("/home/shared_folder/output_plots/gene_length_distribution.pdf")

#gene expression of chromatin states
TPMnormalized_counts_modules <- read.csv("/home/shared_folder/gene_expression_data/DGE_RUES2_data/TPMnormalized_counts_modules.csv") %>% dplyr::select(-X)
TPMcounts = melt(TPMnormalized_counts_modules)
table <- read_excel("/home/shared_folder/gene_expression_data/DGE_RUES2_data/table.xlsx")
colnames(TPMcounts) = c("gene_ID", "color", "gene_name", "gene_type", "sample_ID", "TPM")
TPMcounts = TPMcounts %>% dplyr::left_join(table, by = "sample_ID")
TPMtids = TPMcounts %>% dplyr::mutate(rep = ifelse(gene_ID %in% rep_tids_genes$geneID, "rep", ""),
                                      act = ifelse(gene_ID %in% act_tids_genes$geneID, "act", ""),
                                      biv = ifelse(gene_ID %in% biv_tids_genes$geneID, "biv", "")) %>%
  dplyr::filter(rep != "" | act != "" | biv != "") %>% dplyr::mutate(signature = paste0(rep, act, biv)) %>% dplyr::filter(condition == "WT")

ggplot(TPMtids) + geom_violin(aes(x = signature, y = log2(TPM))) + geom_boxplot(aes(x = signature, y = log2(TPM),  fill = signature)) +
  scale_fill_manual(values = c("#FEAF16", "#1CBE4F", "#56B4E9"))
ggsave("/home/shared_folder/output_plots/gene_expression_distribution_RUES2_WT.pdf")

TPMnormalized_counts_modules <- read.csv("/home/shared_folder/gene_expression_data/DGE_WTC11_data/TPMnormalized_counts_modules.csv") %>% dplyr::select(-X)
TPMcounts = melt(TPMnormalized_counts_modules)
table <- read_excel("/home/shared_folder/gene_expression_data/DGE_WTC11_data/table.xlsx")
colnames(TPMcounts) = c("gene_ID", "color", "gene_name", "gene_type", "sample_ID", "TPM")
TPMcounts = TPMcounts %>% dplyr::left_join(table, by = "sample_ID")
TPMtids = TPMcounts %>% dplyr::mutate(rep = ifelse(gene_ID %in% rep_tids_genes$geneID, "rep", ""),
                                      act = ifelse(gene_ID %in% act_tids_genes$geneID, "act", ""),
                                      biv = ifelse(gene_ID %in% biv_tids_genes$geneID, "biv", "")) %>%
  dplyr::filter(rep != "" | act != "" | biv != "") %>% dplyr::mutate(signature = paste0(rep, act, biv)) %>% dplyr::filter(condition == "WT")

ggplot(TPMtids) + geom_violin(aes(x = signature, y = log2(TPM))) + geom_boxplot(aes(x = signature, y = log2(TPM),  fill = signature)) +
  scale_fill_manual(values = c("#FEAF16", "#1CBE4F", "#56B4E9"))
ggsave("/home/shared_folder/output_plots/gene_expression_distribution_WTC11_WT.pdf")

##gene expression in TIDs
WTvsRBM20.KO <- read.csv("/home/shared_folder/chromatin_data/DGE_RUES2_data/WTvsRBM20.KO.csv")
WTvsTTN.PROM.KO <- read.csv("/home/shared_folder/chromatin_data/DGE_RUES2_data/WTvsTTN.PROM.KO.csv")
WTvsTTN.KO <- read.csv("/home/shared_folder/chromatin_data/DGE_RUES2_data/WTvsTTN.KO.csv")

factory_genes_RUES2 = WTvsRBM20.KO %>% dplyr::select(c("gene_name", "gene_ID", "logFC", "adj.P.Val")) %>% dplyr::mutate(logFC = -logFC)
colnames(factory_genes_RUES2) = c("gene_name", "gene_ID", "logFC_RBM20KO", "adj.P.Val_RBM20KO")
ttnpromKO_genes = WTvsTTN.PROM.KO %>% dplyr::select(c("gene_name", "gene_ID", "logFC", "adj.P.Val")) %>% dplyr::mutate(logFC = -logFC)
colnames(ttnpromKO_genes) = c("gene_name", "gene_ID", "logFC_TTNpromKO", "adj.P.Val_TTNpromKO")
factory_genes_RUES2 = factory_genes_RUES2 %>% dplyr::left_join(ttnpromKO_genes, by = c("gene_name", "gene_ID"))
ttnKO_genes = WTvsTTN.KO %>% dplyr::select(c("gene_name", "gene_ID", "logFC", "adj.P.Val")) %>% dplyr::mutate(logFC = -logFC)
colnames(ttnKO_genes) = c("gene_name", "gene_ID", "logFC_TTNKO", "adj.P.Val_TTNKO")
factory_genes_RUES2 = factory_genes_RUES2 %>% dplyr::left_join(ttnKO_genes, by = c("gene_name", "gene_ID"))

p.val_table1 = data.frame(state = c(), RBM20KO = c(), pTTNKO = c(), perc = c(), ref = c(), p.val = c())
for(i in c("active", "repressive", "bivalent")){
  if(i == "active"){genes = act_tids_genes$geneID}
  else if(i == "repressive"){genes = rep_tids_genes$geneID}
  else if(i == "bivalent"){genes = biv_tids_genes$geneID}
  
  intersection = factory_genes_RUES2 %>% dplyr::mutate(sign = ifelse(logFC_RBM20KO > 0 & logFC_TTNpromKO > 0, "up", "z")) %>%
    dplyr::mutate(tids = ifelse(factory_genes_RUES2$gene_ID %in% genes, i, "z")) %>% dplyr::select(c("sign", "tids"))
  table = table(intersection)
  perc = table[1,1] / (table[1,1]+table[2,1]) * 100
  test <- fisher.test(table, alternative = "two.sided")
  p.val_table1 = rbind(p.val_table1, data.frame(state = c(i), RBM20KO = c("RBM20KO_upregulated"), pTTNKO = ("pTTNKO_upregulated"),
                                                perc = c(perc), ref = c(table[1,2] / (table[1,2]+table[2,2]) * 100), p.val = c(test[["p.value"]])))
  
  intersection = factory_genes_RUES2 %>% dplyr::mutate(sign = ifelse(logFC_RBM20KO > 0 & logFC_TTNpromKO < 0, "up", "z")) %>%
    dplyr::mutate(tids = ifelse(factory_genes_RUES2$gene_ID %in% genes, i, "z")) %>% dplyr::select(c("sign", "tids"))
  table = table(intersection)
  perc = table[1,1] / (table[1,1]+table[2,1]) * 100
  test <- fisher.test(table, alternative = "two.sided")
  p.val_table1 = rbind(p.val_table1, data.frame(state = c(i), RBM20KO = c("RBM20KO_upregulated"), pTTNKO = ("pTTNKO_downregulated"),
                                                perc = c(perc), ref = c(table[1,2] / (table[1,2]+table[2,2]) * 100), p.val = c(test[["p.value"]])))
  
  intersection = factory_genes_RUES2 %>% dplyr::mutate(sign = ifelse(logFC_RBM20KO < 0 & logFC_TTNpromKO > 0, "up", "z")) %>%
    dplyr::mutate(tids = ifelse(factory_genes_RUES2$gene_ID %in% genes, i, "z")) %>% dplyr::select(c("sign", "tids"))
  table = table(intersection)
  perc = table[1,1] / (table[1,1]+table[2,1]) * 100
  test <- fisher.test(table, alternative = "two.sided")
  p.val_table1 = rbind(p.val_table1, data.frame(state = c(i), RBM20KO = c("RBM20KO_downregulated"), pTTNKO = ("pTTNKO_upregulated"),
                                                perc = c(perc), ref = c(table[1,2] / (table[1,2]+table[2,2]) * 100), p.val = c(test[["p.value"]])))
  
  intersection = factory_genes_RUES2 %>% dplyr::mutate(sign = ifelse(logFC_RBM20KO < 0 & logFC_TTNpromKO < 0, "up", "z")) %>%
    dplyr::mutate(tids = ifelse(factory_genes_RUES2$gene_ID %in% genes, i, "z")) %>% dplyr::select(c("sign", "tids"))
  table = table(intersection)
  perc = table[1,1] / (table[1,1]+table[2,1]) * 100
  test <- fisher.test(table, alternative = "two.sided")
  p.val_table1 = rbind(p.val_table1, data.frame(state = c(i), RBM20KO = c("RBM20KO_downregulated"), pTTNKO = ("pTTNKO_downregulated"),
                                                perc = c(perc), ref = c(table[1,2] / (table[1,2]+table[2,2]) * 100), p.val = c(test[["p.value"]])))
  
}

p.val_table1 = p.val_table1 %>% dplyr::mutate(enrichment = perc/ref)
p.val_table1 = p.val_table1 %>% dplyr::mutate(q.val = p.adjust(p.val, method = "BH"))

p.val_table2 = data.frame(state = c(), RBM20KO = c(), TTNKO = c(), perc = c(), ref = c(), p.val = c())
for(i in c("active", "repressive", "bivalent")){
  if(i == "active"){genes = act_tids_genes$geneID}
  else if(i == "repressive"){genes = rep_tids_genes$geneID}
  else if(i == "bivalent"){genes = biv_tids_genes$geneID}
  
  intersection = factory_genes_RUES2 %>% dplyr::mutate(sign = ifelse(logFC_RBM20KO > 0 & logFC_TTNKO > 0, "up", "z")) %>%
    dplyr::mutate(tids = ifelse(factory_genes_RUES2$gene_ID %in% genes, i, "z")) %>% dplyr::select(c("sign", "tids"))
  table = table(intersection)
  perc = table[1,1] / (table[1,1]+table[2,1]) * 100
  test <- fisher.test(table, alternative = "two.sided")
  p.val_table2 = rbind(p.val_table2, data.frame(state = c(i), RBM20KO = c("RBM20KO_upregulated"), TTNKO = ("TTNKO_upregulated"),
                                                perc = c(perc), ref = c(table[1,2] / (table[1,2]+table[2,2]) * 100), p.val = c(test[["p.value"]])))
  
  intersection = factory_genes_RUES2 %>% dplyr::mutate(sign = ifelse(logFC_RBM20KO > 0 & logFC_TTNKO < 0, "up", "z")) %>%
    dplyr::mutate(tids = ifelse(factory_genes_RUES2$gene_ID %in% genes, i, "z")) %>% dplyr::select(c("sign", "tids"))
  table = table(intersection)
  perc = table[1,1] / (table[1,1]+table[2,1]) * 100
  test <- fisher.test(table, alternative = "two.sided")
  p.val_table2 = rbind(p.val_table2, data.frame(state = c(i), RBM20KO = c("RBM20KO_upregulated"), TTNKO = ("TTNKO_downregulated"),
                                                perc = c(perc), ref = c(table[1,2] / (table[1,2]+table[2,2]) * 100), p.val = c(test[["p.value"]])))
  
  intersection = factory_genes_RUES2 %>% dplyr::mutate(sign = ifelse(logFC_RBM20KO < 0 & logFC_TTNKO > 0, "up", "z")) %>%
    dplyr::mutate(tids = ifelse(factory_genes_RUES2$gene_ID %in% genes, i, "z")) %>% dplyr::select(c("sign", "tids"))
  table = table(intersection)
  perc = table[1,1] / (table[1,1]+table[2,1]) * 100
  test <- fisher.test(table, alternative = "two.sided")
  p.val_table2 = rbind(p.val_table2, data.frame(state = c(i), RBM20KO = c("RBM20KO_downregulated"), TTNKO = ("TTNKO_upregulated"),
                                                perc = c(perc), ref = c(table[1,2] / (table[1,2]+table[2,2]) * 100), p.val = c(test[["p.value"]])))
  
  intersection = factory_genes_RUES2 %>% dplyr::mutate(sign = ifelse(logFC_RBM20KO < 0 & logFC_TTNKO < 0, "up", "z")) %>%
    dplyr::mutate(tids = ifelse(factory_genes_RUES2$gene_ID %in% genes, i, "z")) %>% dplyr::select(c("sign", "tids"))
  table = table(intersection)
  perc = table[1,1] / (table[1,1]+table[2,1]) * 100
  test <- fisher.test(table, alternative = "two.sided")
  p.val_table2 = rbind(p.val_table2, data.frame(state = c(i), RBM20KO = c("RBM20KO_downregulated"), TTNKO = ("TTNKO_downregulated"),
                                                perc = c(perc), ref = c(table[1,2] / (table[1,2]+table[2,2]) * 100), p.val = c(test[["p.value"]])))
  
}

p.val_table2 = p.val_table2 %>% dplyr::mutate(enrichment = perc/ref)
p.val_table2 = p.val_table2 %>% dplyr::mutate(q.val = p.adjust(p.val, method = "BH"))

enrichment1_act = dcast(p.val_table1 %>% dplyr::filter(state == "active"), pTTNKO ~ RBM20KO+state, value.var = "enrichment") 
rownames(enrichment1_act) = enrichment1_act$pTTNKO
enrichment1_act = enrichment1_act %>% dplyr::select(-pTTNKO)
enrichment1_rep = dcast(p.val_table1 %>% dplyr::filter(state == "repressive"), pTTNKO ~ RBM20KO+state, value.var = "enrichment")
rownames(enrichment1_rep) = enrichment1_rep$pTTNKO
enrichment1_rep = enrichment1_rep %>% dplyr::select(-pTTNKO)
enrichment1_biv = dcast(p.val_table1 %>% dplyr::filter(state == "bivalent"), pTTNKO ~ RBM20KO+state, value.var = "enrichment")
rownames(enrichment1_biv) = enrichment1_biv$pTTNKO
enrichment1_biv = enrichment1_biv %>% dplyr::select(-pTTNKO)

enrichment2_act = dcast(p.val_table2 %>% dplyr::filter(state == "active"), TTNKO ~ RBM20KO+state, value.var = "enrichment")
rownames(enrichment2_act) = enrichment2_act$TTNKO
enrichment2_act = enrichment2_act %>% dplyr::select(-TTNKO)
enrichment2_rep = dcast(p.val_table2 %>% dplyr::filter(state == "repressive"), TTNKO ~ RBM20KO+state, value.var = "enrichment")
rownames(enrichment2_rep) = enrichment2_rep$TTNKO
enrichment2_rep = enrichment2_rep %>% dplyr::select(-TTNKO)
enrichment2_biv = dcast(p.val_table2 %>% dplyr::filter(state == "bivalent"), TTNKO ~ RBM20KO+state, value.var = "enrichment")
rownames(enrichment2_biv) = enrichment2_biv$TTNKO
enrichment2_biv = enrichment2_biv %>% dplyr::select(-TTNKO)

enrichment_act = rbind(enrichment1_act, enrichment2_act)
enrichment_rep = rbind(enrichment1_rep, enrichment2_rep)
enrichment_biv = rbind(enrichment1_biv, enrichment2_biv)
enrichment = cbind(enrichment_act, enrichment_rep, enrichment_biv)
enrichment = enrichment[c(2, 1, 4, 3),]

annotation_col = data.frame(name = colnames(enrichment), chrom_state = c("active", "active", "repressive", "repressive", "bivalent", "bivalent"))
rownames(annotation_col) = annotation_col$name
annotation_col = annotation_col %>% dplyr::select(-name)

annotation_colors = list(chrom_state = c(active = "#FEAF16", repressive = "#56B4E9", bivalent = "#1CBE4F"))

p.val_table1 = p.val_table1 %>% dplyr::mutate(sign = ifelse(q.val < 0.05, " *", "")) %>% dplyr::mutate(n = paste0(signif(enrichment, digits = 3), sign))
n1_act = dcast(p.val_table1 %>% dplyr::filter(state == "active"), pTTNKO ~ RBM20KO+state, value.var = "n") 
rownames(n1_act) = n1_act$pTTNKO
n1_act = n1_act %>% dplyr::select(-pTTNKO)
n1_rep = dcast(p.val_table1 %>% dplyr::filter(state == "repressive"), pTTNKO ~ RBM20KO+state, value.var = "n")
rownames(n1_rep) = n1_rep$pTTNKO
n1_rep = n1_rep %>% dplyr::select(-pTTNKO)
n1_biv = dcast(p.val_table1 %>% dplyr::filter(state == "bivalent"), pTTNKO ~ RBM20KO+state, value.var = "n")
rownames(n1_biv) = n1_biv$pTTNKO
n1_biv = n1_biv %>% dplyr::select(-pTTNKO)

p.val_table2 = p.val_table2 %>% dplyr::mutate(sign = ifelse(q.val < 0.05, " *", "")) %>% dplyr::mutate(n = paste0(signif(enrichment, digits = 3), sign))
n2_act = dcast(p.val_table2 %>% dplyr::filter(state == "active"), TTNKO ~ RBM20KO+state, value.var = "n") 
rownames(n2_act) = n2_act$TTNKO
n2_act = n2_act %>% dplyr::select(-TTNKO)
n2_rep = dcast(p.val_table2 %>% dplyr::filter(state == "repressive"), TTNKO ~ RBM20KO+state, value.var = "n")
rownames(n2_rep) = n2_rep$TTNKO
n2_rep = n2_rep %>% dplyr::select(-TTNKO)
n2_biv = dcast(p.val_table2 %>% dplyr::filter(state == "bivalent"), TTNKO ~ RBM20KO+state, value.var = "n")
rownames(n2_biv) = n2_biv$TTNKO
n2_biv = n2_biv %>% dplyr::select(-TTNKO)

n_act = rbind(n1_act, n2_act)
n_rep = rbind(n1_rep, n2_rep)
n_biv = rbind(n1_biv, n2_biv)
numbers = cbind(n_act, n_rep, n_biv)
numbers = numbers[c(2, 1, 4, 3),]

#heatmap enrichment
as.ggplot(pheatmap(enrichment, annotation_col = annotation_col,
                   cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = numbers,
                   labels_col = rep(c("RBM20KO_downregulated", "RBM20KO_upregulated"), times = 3),
                   lebels_row = c("pTTNKO_downregulated", "pTTNKO_upregulated", "TTNKO_upregulated", "TTNKO_downregulated"),
                   fontsize_row = 7, fontsize_col = 7, gaps_row = 2, gaps_col = c(2, 4), annotation_colors = annotation_colors,
                   color = viridis::viridis(n = 25, option = "H"), number_color = "white"))
ggsave("/home/shared_folder/output_plots/enrichment_table_RUES2.pdf")

#density of double upregulated genes pTTNKO
factory_genes_RUES2 = factory_genes_RUES2 %>% dplyr::mutate(name = ifelse(gene_ID %in% biv_tids_genes$geneID, gene_name, "")) %>%
  dplyr::mutate(biv = ifelse(name != "", "biv", "")) %>% dplyr::arrange(biv)
ggplot(factory_genes_RUES2 %>% dplyr::filter(biv == "biv"), aes(x = logFC_RBM20KO, y = logFC_TTNpromKO)) + stat_density_2d(aes(color = biv), size = 1, color = "#1CBE4F") + geom_point(size = 1) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + geom_hline(aes(yintercept = 0), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = name), max.overlaps = 10,	size = 2, hjust = 1.2) +
  scale_x_continuous(limits = c(-2.5, 4)) + scale_y_continuous(limits = c(-2.5, 4)) +
  ggtitle(paste0("q-val RBM20KO up TTN prom KO up = ", signif(p.val_table1$q.val[9], digits = 3), "\nenrichment score = ", signif(p.val_table1$enrichment[9], digits = 3))) +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/biv_density_genes_pTTNKO_RUES2.pdf")

#density of double upregulated genes pTTNKO
ggplot(factory_genes_RUES2 %>% dplyr::filter(biv == "biv"), aes(x = logFC_RBM20KO, y = logFC_TTNKO)) + stat_density_2d(aes(color = biv), size = 1, color = "#1CBE4F") + geom_point(size = 1) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") + geom_hline(aes(yintercept = 0), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = name), max.overlaps = 10,	size = 2, hjust = 1.2) +
  scale_x_continuous(limits = c(-2.5, 4)) + scale_y_continuous(limits = c(-2.5, 4)) +
  ggtitle(paste0("q-val RBM20KO up TTN KO up = ", signif(p.val_table2$q.val[9], digits = 3), "\nenrichment score = ", signif(p.val_table2$enrichment[9], digits = 3))) +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/biv_density_genes_TTNKO_RUES2.pdf")

