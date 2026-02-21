library(dplyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(readxl)
library(rtracklayer)
library(ggplotify)
library(pheatmap)
library(fgsea)
library(transite)
library(BSgenome)
library(biomaRt)


theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size=15, face="bold", margin = margin(10,0,10,0))))

#functions for splicing variants importation
get_JCcounts <- function(str){
  val = strsplit(str, ",")[[1]]
  sum = 0
  for(i in 1:length(val)){
    sum = sum + as.numeric(val[i])
  }
  return(sum)
  
}

get_PSI_n <- function(str, n){
  val = strsplit(str, ",")[[1]]
  return(val[n])
}

get_PSI <- function(str){
  val = strsplit(str, ",")[[1]]
  sum = 0
  n = 0
  for(i in 1:length(val)){
    if(val[i] != "NA"){
      n = n + 1
      sum = sum + as.numeric(val[i])
    }
  }
  return(sum/n)
}

bases = c("A", "T", "C", "G")
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

##splicing alterations in pulldown targets according to RUES2 CMs data
geneInfo <- read.table("/home/shared_folder/pulldown_data/transcripts_data/geneInfo.tab", quote="\"", comment.char="", skip = 1)
colnames(geneInfo) = c("Gene", "Symbol", "gene_type")

transC_lists <- read_excel("/home/shared_folder/chromatin_data/eCLIP_data/transC_lists.xlsx") %>% dplyr::filter(HC_eCLIP == "Yes")

WT_OMAP <- read_excel("/home/shared_folder/pulldown_data/transcripts_data/o_map_WT.xlsx") %>%
  dplyr::left_join(geneInfo, by = "Gene") %>% dplyr::mutate(eCLIP = ifelse(Symbol %in% transC_lists$Gene, "eCLIP", "z"))
KO_OMAP <- read_excel("/home/shared_folder/pulldown_data/transcripts_data/o_map_RBM20KO.xlsx") %>%
  dplyr::left_join(geneInfo, by = "Gene") %>% dplyr::mutate(eCLIP = ifelse(Symbol %in% transC_lists$Gene, "eCLIP", "z"))

A3SS.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/A3SS.MATS.JC.txt"))
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A3SS.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(A3SS.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
A3SS.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(A3SS.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(A3SS.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")
A5SS.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/A5SS.MATS.JC.txt"))
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A5SS.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(A5SS.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
A5SS.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(A5SS.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(A5SS.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")
RI.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/RI.MATS.JC.txt"))
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_riE:", riExonStart_0base, "-", riExonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
RI.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(RI.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
RI.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(RI.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(RI.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")
MXE.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/MXE.MATS.JC.txt"))
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E1:", X1stExonStart_0base, "-", X1stExonEnd, "_E2:", X2ndExonStart_0base, "-", X2ndExonEnd, "_uE:", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
MXE.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(MXE.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
MXE.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(MXE.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(MXE.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")
SE.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/SE.MATS.JC.txt"))
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E:", exonStart_0base, "-", exonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
SE.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(SE.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
SE.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(SE.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(SE.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

merge = rbind(A3SS.WT_vs_RBM20KO %>% dplyr::mutate(event = "A3SS"),
              A5SS.WT_vs_RBM20KO %>% dplyr::mutate(event = "A5SS"),
              RI.WT_vs_RBM20KO %>% dplyr::mutate(event = "RI"),
              MXE.WT_vs_RBM20KO %>% dplyr::mutate(event = "MXE"),
              SE.WT_vs_RBM20KO %>% dplyr::mutate(event = "SE"))
merge[is.na(merge)] = ""
merge = merge %>% dplyr::mutate(eCLIP = ifelse(geneSymbol %in% transC_lists$Gene & geneSymbol != "TTN", "RBM20_target", "non_target")) %>%
  dplyr::mutate(eCLIP = ifelse(abs(delta_PSI_RBM20KOvsWT) > 0.10 & FDR_RBM20KOvsWT < 0.05, eCLIP, "insignificant")) %>%
  dplyr::mutate(WT = ifelse(geneID %in% WT_OMAP$Gene, "pulldown", "non_target")) %>%
  dplyr::mutate(WT = ifelse(abs(delta_PSI_RBM20KOvsWT) > 0.10 & FDR_RBM20KOvsWT < 0.05, WT, "insignificant")) %>% 
  dplyr::mutate(KO = ifelse(geneID %in% KO_OMAP$Gene, "pulldown", "non_target")) %>%
  dplyr::mutate(KO = ifelse(abs(delta_PSI_RBM20KOvsWT) > 0.10 & FDR_RBM20KOvsWT < 0.05, KO, "insignificant"))

#pulldown targets in differentially splicied genes - figure 5 panel D
ggplot(merge %>% dplyr::arrange(WT), aes(x = delta_PSI_RBM20KOvsWT, y = -log10(FDR_RBM20KOvsWT))) + geom_point(aes(fill = WT), shape = 21, size = 3, alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("grey", "white", "#3283FE")) + scale_x_continuous(limits = c(-1.5,1.5)) + ggtitle("WT O-MAP pulldown transcripts") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0))) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = 0.10), linetype = "dashed") + geom_vline(aes(xintercept = -0.10), linetype = "dashed")
ggsave("/home/shared_folder/output_plots/volcano_RUES2_WT_OMAP_targets.pdf")

ggplot(merge %>% dplyr::arrange(KO), aes(x = delta_PSI_RBM20KOvsWT, y = -log10(FDR_RBM20KOvsWT))) + geom_point(aes(fill = KO), shape = 21, size = 3, alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("grey", "white", "#F6222E")) + scale_x_continuous(limits = c(-1.5, 1.5)) + ggtitle("RBM20KO O-MAP pulldown transcripts") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0))) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = 0.10), linetype = "dashed") + geom_vline(aes(xintercept = -0.10), linetype = "dashed")
ggsave("/home/shared_folder/output_plots/volcano_RUES2_RBM20KO_OMAP_targets.pdf")

#RBM20 targtes in pulldown transcripts - figure 5 panel E
merge = merge %>% dplyr::filter(geneSymbol != "TTN") %>% dplyr::mutate(name = ifelse(eCLIP == "RBM20_target", geneSymbol, ""))
ggplot(merge %>% dplyr::filter(WT == "pulldown") %>% dplyr::arrange(desc(eCLIP)), aes(x = delta_PSI_RBM20KOvsWT, y = -log10(FDR_RBM20KOvsWT))) + geom_point(aes(fill = eCLIP, color = eCLIP), shape = 21, size = 3, alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("white", "#009E73")) + scale_x_continuous(limits = c(-1.5,1.5)) + ggtitle("RBM20 targets in WT O-MAP pulldown transcripts") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0))) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = 0.10), linetype = "dashed") + geom_vline(aes(xintercept = -0.10), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = name), max.overlaps = 30,	size = 2, hjust = 1.2)
ggsave("/home/shared_folder/output_plots/volcano_RUES2_WT_eCLIP.pdf")

ggplot(merge %>% dplyr::filter(KO == "pulldown") %>% dplyr::arrange(desc(eCLIP)), aes(x = delta_PSI_RBM20KOvsWT, y = -log10(FDR_RBM20KOvsWT))) + geom_point(aes(fill = eCLIP, color = eCLIP), shape = 21, size = 3, alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("white", "#009E73")) + scale_x_continuous(limits = c(-1.5,1.5)) + ggtitle("RBM20 targets in RBM20KO O-MAP pulldown transcripts") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0))) + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_vline(aes(xintercept = 0.10), linetype = "dashed") + geom_vline(aes(xintercept = -0.10), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = name), max.overlaps = 30,	size = 2, hjust = 1.2)
ggsave("/home/shared_folder/output_plots/volcano_RUES2_KO_eCLIP.pdf")

##pulldown splicing
A3SS.MATS.JC <- read.delim("/home/shared_folder/splicing_data/pulldown_rMATS_data/A3SS.MATS.JC.txt") %>% dplyr::select(c(geneSymbol, IncLevelDifference, FDR)) %>% dplyr::mutate(event = "A3SS")
A5SS.MATS.JC <- read.delim("/home/shared_folder/splicing_data/pulldown_rMATS_data/A5SS.MATS.JC.txt") %>% dplyr::select(c(geneSymbol, IncLevelDifference, FDR)) %>% dplyr::mutate(event = "A5SS")
RI.MATS.JC <- read.delim("/home/shared_folder/splicing_data/pulldown_rMATS_data/RI.MATS.JC.txt") %>% dplyr::select(c(geneSymbol, IncLevelDifference, FDR)) %>% dplyr::mutate(event = "RI")
MXE.MATS.JC <- read.delim("/home/shared_folder/splicing_data/pulldown_rMATS_data/MXE.MATS.JC.txt") %>% dplyr::select(c(geneSymbol, IncLevelDifference, FDR)) %>% dplyr::mutate(event = "MXE")
SE.MATS.JC <- read.delim("/home/shared_folder/splicing_data/pulldown_rMATS_data/SE.MATS.JC.txt") %>% dplyr::select(c(geneSymbol, IncLevelDifference, FDR)) %>% dplyr::mutate(event = "SE")
merge_pulldown = rbind(A3SS.MATS.JC, A5SS.MATS.JC, RI.MATS.JC, MXE.MATS.JC, SE.MATS.JC)
merge_pulldown[is.na(merge_pulldown)] = ""
merge_pulldown = merge_pulldown %>%
  dplyr::mutate(pulldown = ifelse(geneSymbol %in% WT_OMAP$Symbol & !geneSymbol %in% KO_OMAP$Symbol, "WT_specific",
                                  ifelse(geneSymbol %in% KO_OMAP$Symbol & !geneSymbol %in% WT_OMAP$Symbol, "RBM20KO_specific", ""))) %>%
  dplyr::mutate(pulldown = ifelse(FDR < 0.05 & abs(IncLevelDifference) > 0.10, pulldown, "insignificant")) %>%
  dplyr::mutate(eCLIP = ifelse(geneSymbol == "TTN", "TTN", ifelse(geneSymbol %in% transC_lists$Gene, "RBM20_target", "no_target"))) %>% 
  dplyr::mutate(eCLIP = ifelse(FDR < 0.05 & abs(IncLevelDifference) > 0.10, eCLIP, "insignificant")) %>%
  dplyr::mutate(event_type = ifelse(FDR < 0.05 & abs(IncLevelDifference) > 0.10, event, "insignificant"))

ggplot(merge_pulldown %>% dplyr::arrange(pulldown)) + geom_point(aes(x = IncLevelDifference, y = -log10(FDR), fill = event_type), shape = 21, size = 3, alpha = 0.7) +
  scale_fill_manual(values = c("#E69F00", "#DEA0FD", "grey", "#ff6db6", "#009E73", "#56B4E9")) +
  geom_vline(aes(xintercept = -0.1), linetype = "dashed") + geom_vline(aes(xintercept = 0.1), linetype = "dashed") + 
  geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed")
ggsave("/home/shared_folder/output_plots/all_events_volcano.pdf")


##WTC11 splicing analysis
A3SS.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KO/A3SS.MATS.JC.txt"))
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A3SS.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(A3SS.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
A3SS.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(A3SS.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
A3SS.WT_vs_RBM20KO$WT_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel1, get_PSI))
A3SS.WT_vs_RBM20KO$RBM20KO_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel2, get_PSI))
A3SS.WT_vs_RBM20KO$WT_n1_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 1))
A3SS.WT_vs_RBM20KO$WT_n3_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 2))
A3SS.WT_vs_RBM20KO$WT_n4_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 3))
A3SS.WT_vs_RBM20KO$WT_n5_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 4))
A3SS.WT_vs_RBM20KO$WT_n6_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 5))
A3SS.WT_vs_RBM20KO$RBM20KO_n3_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 1))
A3SS.WT_vs_RBM20KO$RBM20KO_n4_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 2))
A3SS.WT_vs_RBM20KO$RBM20KO_n5_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 3))
A3SS.WT_vs_RBM20KO$RBM20KO_n6_PSI = as.character(lapply(A3SS.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 4))
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                            RBM20KO_PSI, RBM20KO_n3_PSI, RBM20KO_n4_PSI, RBM20KO_n5_PSI, RBM20KO_n6_PSI, IncLevelDifference, FDR))
colnames(A3SS.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                 "RBM20KO_PSI", "RBM20KO_n3_PSI", "RBM20KO_n4_PSI", "RBM20KO_n5_PSI", "RBM20KO_n6_PSI", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

A3SS.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_TTNKO/A3SS.MATS.JC.txt")
A3SS.WT_vs_TTNpromKO = A3SS.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A3SS.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(A3SS.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
A3SS.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(A3SS.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
A3SS.WT_vs_TTNpromKO$WT_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel1, get_PSI))
A3SS.WT_vs_TTNpromKO$TTNpromKO_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel2, get_PSI))
A3SS.WT_vs_TTNpromKO$WT_n1_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 1))
A3SS.WT_vs_TTNpromKO$WT_n3_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 2))
A3SS.WT_vs_TTNpromKO$WT_n4_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 3))
A3SS.WT_vs_TTNpromKO$WT_n5_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 4))
A3SS.WT_vs_TTNpromKO$WT_n6_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 5))
A3SS.WT_vs_TTNpromKO$TTNpromKO_n1_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 1))
A3SS.WT_vs_TTNpromKO$TTNpromKO_n3_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 2))
A3SS.WT_vs_TTNpromKO$TTNpromKO_n4_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 3))
A3SS.WT_vs_TTNpromKO$TTNpromKO_n5_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 4))
A3SS.WT_vs_TTNpromKO$TTNpromKO_n6_PSI = as.character(lapply(A3SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 5))
A3SS.WT_vs_TTNpromKO = A3SS.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A3SS.WT_vs_TTNpromKO = A3SS.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                                TTNpromKO_PSI, TTNpromKO_n1_PSI, TTNpromKO_n3_PSI, TTNpromKO_n4_PSI, TTNpromKO_n5_PSI, TTNpromKO_n6_PSI, IncLevelDifference, FDR))
colnames(A3SS.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                   "TTNpromKO_PSI", "TTNpromKO_n1_PSI", "TTNpromKO_n3_PSI", "TTNpromKO_n4_PSI", "TTNpromKO_n5_PSI", "TTNpromKO_n6_PSI", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

A3SS.WT_vs_doubleKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KOTTNKO/A3SS.MATS.JC.txt") 
A3SS.WT_vs_doubleKO = A3SS.WT_vs_doubleKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A3SS.WT_vs_doubleKO$JC_count_WT = as.character(lapply(A3SS.WT_vs_doubleKO$IJC_SAMPLE_1, get_JCcounts))
A3SS.WT_vs_doubleKO$JC_count_doubleKO = as.character(lapply(A3SS.WT_vs_doubleKO$IJC_SAMPLE_2, get_JCcounts))
A3SS.WT_vs_doubleKO$WT_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel1, get_PSI))
A3SS.WT_vs_doubleKO$doubleKO_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel2, get_PSI))
A3SS.WT_vs_doubleKO$WT_n1_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 1))
A3SS.WT_vs_doubleKO$WT_n3_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 2))
A3SS.WT_vs_doubleKO$WT_n4_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 3))
A3SS.WT_vs_doubleKO$WT_n5_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 4))
A3SS.WT_vs_doubleKO$WT_n6_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 5))
A3SS.WT_vs_doubleKO$doubleKO_n1_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 1))
A3SS.WT_vs_doubleKO$doubleKO_n3_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 2))
A3SS.WT_vs_doubleKO$doubleKO_n4_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 3))
A3SS.WT_vs_doubleKO$doubleKO_n5_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 4))
A3SS.WT_vs_doubleKO$doubleKO_n6_PSI = as.character(lapply(A3SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 5))
A3SS.WT_vs_doubleKO = A3SS.WT_vs_doubleKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A3SS.WT_vs_doubleKO = A3SS.WT_vs_doubleKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_doubleKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                              doubleKO_PSI, doubleKO_n1_PSI, doubleKO_n3_PSI, doubleKO_n4_PSI, doubleKO_n5_PSI, doubleKO_n6_PSI, IncLevelDifference, FDR))
colnames(A3SS.WT_vs_doubleKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_doubleKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                  "doubleKO_PSI", "doubleKO_n1_PSI", "doubleKO_n3_PSI", "doubleKO_n4_PSI", "doubleKO_n5_PSI", "doubleKO_n6_PSI", "delta_PSI_doubleKOvsWT", "FDR_doubleKOvsWT")

A3SS_intersect = A3SS.WT_vs_RBM20KO %>% dplyr::inner_join(A3SS.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT")) %>%
  dplyr::inner_join(A3SS.WT_vs_doubleKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI","JC_WT"))
A3SS_intersect = A3SS_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_doubleKO)) %>% dplyr::filter(JC_sum > 19)

A5SS.WT_vs_RBM20KO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KO/A5SS.MATS.JC.txt") 
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A5SS.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(A5SS.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
A5SS.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(A5SS.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
A5SS.WT_vs_RBM20KO$WT_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel1, get_PSI))
A5SS.WT_vs_RBM20KO$RBM20KO_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel2, get_PSI))
A5SS.WT_vs_RBM20KO$WT_n1_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 1))
A5SS.WT_vs_RBM20KO$WT_n3_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 2))
A5SS.WT_vs_RBM20KO$WT_n4_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 3))
A5SS.WT_vs_RBM20KO$WT_n5_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 4))
A5SS.WT_vs_RBM20KO$WT_n6_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 5))
A5SS.WT_vs_RBM20KO$RBM20KO_n3_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 1))
A5SS.WT_vs_RBM20KO$RBM20KO_n4_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 2))
A5SS.WT_vs_RBM20KO$RBM20KO_n5_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 3))
A5SS.WT_vs_RBM20KO$RBM20KO_n6_PSI = as.character(lapply(A5SS.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 4))
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                            RBM20KO_PSI, RBM20KO_n3_PSI, RBM20KO_n4_PSI, RBM20KO_n5_PSI, RBM20KO_n6_PSI, IncLevelDifference, FDR))
colnames(A5SS.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                 "RBM20KO_PSI", "RBM20KO_n3_PSI", "RBM20KO_n4_PSI", "RBM20KO_n5_PSI", "RBM20KO_n6_PSI", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

A5SS.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_TTNKO/A5SS.MATS.JC.txt")
A5SS.WT_vs_TTNpromKO = A5SS.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A5SS.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(A5SS.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
A5SS.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(A5SS.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
A5SS.WT_vs_TTNpromKO$WT_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel1, get_PSI))
A5SS.WT_vs_TTNpromKO$TTNpromKO_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel2, get_PSI))
A5SS.WT_vs_TTNpromKO$WT_n1_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 1))
A5SS.WT_vs_TTNpromKO$WT_n3_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 2))
A5SS.WT_vs_TTNpromKO$WT_n4_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 3))
A5SS.WT_vs_TTNpromKO$WT_n5_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 4))
A5SS.WT_vs_TTNpromKO$WT_n6_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 5))
A5SS.WT_vs_TTNpromKO$TTNpromKO_n1_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 1))
A5SS.WT_vs_TTNpromKO$TTNpromKO_n3_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 2))
A5SS.WT_vs_TTNpromKO$TTNpromKO_n4_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 3))
A5SS.WT_vs_TTNpromKO$TTNpromKO_n5_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 4))
A5SS.WT_vs_TTNpromKO$TTNpromKO_n6_PSI = as.character(lapply(A5SS.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 5))
A5SS.WT_vs_TTNpromKO = A5SS.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A5SS.WT_vs_TTNpromKO = A5SS.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                                TTNpromKO_PSI, TTNpromKO_n1_PSI, TTNpromKO_n3_PSI, TTNpromKO_n4_PSI, TTNpromKO_n5_PSI, TTNpromKO_n6_PSI, IncLevelDifference, FDR))
colnames(A5SS.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                   "TTNpromKO_PSI", "TTNpromKO_n1_PSI", "TTNpromKO_n3_PSI", "TTNpromKO_n4_PSI", "TTNpromKO_n5_PSI", "TTNpromKO_n6_PSI", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

A5SS.WT_vs_doubleKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KOTTNKO/A5SS.MATS.JC.txt") 
A5SS.WT_vs_doubleKO = A5SS.WT_vs_doubleKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A5SS.WT_vs_doubleKO$JC_count_WT = as.character(lapply(A5SS.WT_vs_doubleKO$IJC_SAMPLE_1, get_JCcounts))
A5SS.WT_vs_doubleKO$JC_count_doubleKO = as.character(lapply(A5SS.WT_vs_doubleKO$IJC_SAMPLE_2, get_JCcounts))
A5SS.WT_vs_doubleKO$WT_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel1, get_PSI))
A5SS.WT_vs_doubleKO$doubleKO_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel2, get_PSI))
A5SS.WT_vs_doubleKO$WT_n1_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 1))
A5SS.WT_vs_doubleKO$WT_n3_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 2))
A5SS.WT_vs_doubleKO$WT_n4_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 3))
A5SS.WT_vs_doubleKO$WT_n5_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 4))
A5SS.WT_vs_doubleKO$WT_n6_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 5))
A5SS.WT_vs_doubleKO$doubleKO_n1_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 1))
A5SS.WT_vs_doubleKO$doubleKO_n3_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 2))
A5SS.WT_vs_doubleKO$doubleKO_n4_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 3))
A5SS.WT_vs_doubleKO$doubleKO_n5_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 4))
A5SS.WT_vs_doubleKO$doubleKO_n6_PSI = as.character(lapply(A5SS.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 5))
A5SS.WT_vs_doubleKO = A5SS.WT_vs_doubleKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A5SS.WT_vs_doubleKO = A5SS.WT_vs_doubleKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_doubleKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                              doubleKO_PSI, doubleKO_n1_PSI, doubleKO_n3_PSI, doubleKO_n4_PSI, doubleKO_n5_PSI, doubleKO_n6_PSI, IncLevelDifference, FDR))
colnames(A5SS.WT_vs_doubleKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_doubleKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                  "doubleKO_PSI", "doubleKO_n1_PSI", "doubleKO_n3_PSI", "doubleKO_n4_PSI", "doubleKO_n5_PSI", "doubleKO_n6_PSI", "delta_PSI_doubleKOvsWT", "FDR_doubleKOvsWT")

A5SS_intersect = A5SS.WT_vs_RBM20KO %>% dplyr::inner_join(A5SS.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT")) %>%
  dplyr::inner_join(A5SS.WT_vs_doubleKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT"))
A5SS_intersect = A5SS_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_doubleKO)) %>% dplyr::filter(JC_sum > 19)

RI.WT_vs_RBM20KO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KO/RI.MATS.JC.txt") 
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_riE:", riExonStart_0base, "-", riExonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
RI.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(RI.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
RI.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(RI.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
RI.WT_vs_RBM20KO$WT_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel1, get_PSI))
RI.WT_vs_RBM20KO$RBM20KO_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel2, get_PSI))
RI.WT_vs_RBM20KO$WT_n1_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 1))
RI.WT_vs_RBM20KO$WT_n3_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 2))
RI.WT_vs_RBM20KO$WT_n4_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 3))
RI.WT_vs_RBM20KO$WT_n5_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 4))
RI.WT_vs_RBM20KO$WT_n6_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 5))
RI.WT_vs_RBM20KO$RBM20KO_n3_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 1))
RI.WT_vs_RBM20KO$RBM20KO_n4_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 2))
RI.WT_vs_RBM20KO$RBM20KO_n5_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 3))
RI.WT_vs_RBM20KO$RBM20KO_n6_PSI = as.character(lapply(RI.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 4))
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                        RBM20KO_PSI, RBM20KO_n3_PSI, RBM20KO_n4_PSI, RBM20KO_n5_PSI, RBM20KO_n6_PSI, IncLevelDifference, FDR))
colnames(RI.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                               "RBM20KO_PSI", "RBM20KO_n3_PSI", "RBM20KO_n4_PSI", "RBM20KO_n5_PSI", "RBM20KO_n6_PSI", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

RI.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_TTNKO/RI.MATS.JC.txt")
RI.WT_vs_TTNpromKO = RI.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_riE:", riExonStart_0base, "-", riExonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
RI.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(RI.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
RI.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(RI.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
RI.WT_vs_TTNpromKO$WT_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel1, get_PSI))
RI.WT_vs_TTNpromKO$TTNpromKO_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel2, get_PSI))
RI.WT_vs_TTNpromKO$WT_n1_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 1))
RI.WT_vs_TTNpromKO$WT_n3_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 2))
RI.WT_vs_TTNpromKO$WT_n4_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 3))
RI.WT_vs_TTNpromKO$WT_n5_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 4))
RI.WT_vs_TTNpromKO$WT_n6_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 5))
RI.WT_vs_TTNpromKO$TTNpromKO_n1_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 1))
RI.WT_vs_TTNpromKO$TTNpromKO_n3_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 2))
RI.WT_vs_TTNpromKO$TTNpromKO_n4_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 3))
RI.WT_vs_TTNpromKO$TTNpromKO_n5_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 4))
RI.WT_vs_TTNpromKO$TTNpromKO_n6_PSI = as.character(lapply(RI.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 5))
RI.WT_vs_TTNpromKO = RI.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
RI.WT_vs_TTNpromKO = RI.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                            TTNpromKO_PSI, TTNpromKO_n1_PSI, TTNpromKO_n3_PSI, TTNpromKO_n4_PSI, TTNpromKO_n5_PSI, TTNpromKO_n6_PSI, IncLevelDifference, FDR))
colnames(RI.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                 "TTNpromKO_PSI", "TTNpromKO_n1_PSI", "TTNpromKO_n3_PSI", "TTNpromKO_n4_PSI", "TTNpromKO_n5_PSI", "TTNpromKO_n6_PSI", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

RI.WT_vs_doubleKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KOTTNKO/RI.MATS.JC.txt") 
RI.WT_vs_doubleKO = RI.WT_vs_doubleKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_riE:", riExonStart_0base, "-", riExonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
RI.WT_vs_doubleKO$JC_count_WT = as.character(lapply(RI.WT_vs_doubleKO$IJC_SAMPLE_1, get_JCcounts))
RI.WT_vs_doubleKO$JC_count_doubleKO = as.character(lapply(RI.WT_vs_doubleKO$IJC_SAMPLE_2, get_JCcounts))
RI.WT_vs_doubleKO$WT_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel1, get_PSI))
RI.WT_vs_doubleKO$doubleKO_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel2, get_PSI))
RI.WT_vs_doubleKO$WT_n1_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 1))
RI.WT_vs_doubleKO$WT_n3_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 2))
RI.WT_vs_doubleKO$WT_n4_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 3))
RI.WT_vs_doubleKO$WT_n5_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 4))
RI.WT_vs_doubleKO$WT_n6_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 5))
RI.WT_vs_doubleKO$doubleKO_n1_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 1))
RI.WT_vs_doubleKO$doubleKO_n3_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 2))
RI.WT_vs_doubleKO$doubleKO_n4_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 3))
RI.WT_vs_doubleKO$doubleKO_n5_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 4))
RI.WT_vs_doubleKO$doubleKO_n6_PSI = as.character(lapply(RI.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 5))
RI.WT_vs_doubleKO = RI.WT_vs_doubleKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
RI.WT_vs_doubleKO = RI.WT_vs_doubleKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_doubleKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                          doubleKO_PSI, doubleKO_n1_PSI, doubleKO_n3_PSI, doubleKO_n4_PSI, doubleKO_n5_PSI, doubleKO_n6_PSI, IncLevelDifference, FDR))
colnames(RI.WT_vs_doubleKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_doubleKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                "doubleKO_PSI", "doubleKO_n1_PSI", "doubleKO_n3_PSI", "doubleKO_n4_PSI", "doubleKO_n5_PSI", "doubleKO_n6_PSI", "delta_PSI_doubleKOvsWT", "FDR_doubleKOvsWT")

RI_intersect = RI.WT_vs_RBM20KO %>% dplyr::inner_join(RI.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT")) %>%
  dplyr::inner_join(RI.WT_vs_doubleKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT"))
RI_intersect = RI_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_doubleKO)) %>% dplyr::filter(JC_sum > 19)

MXE.WT_vs_RBM20KO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KO/MXE.MATS.JC.txt") 
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E1:", X1stExonStart_0base, "-", X1stExonEnd, "_E2:", X2ndExonStart_0base, "-", X2ndExonEnd, "_uE:", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
MXE.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(MXE.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
MXE.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(MXE.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
MXE.WT_vs_RBM20KO$WT_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel1, get_PSI))
MXE.WT_vs_RBM20KO$RBM20KO_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel2, get_PSI))
MXE.WT_vs_RBM20KO$WT_n1_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 1))
MXE.WT_vs_RBM20KO$WT_n3_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 2))
MXE.WT_vs_RBM20KO$WT_n4_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 3))
MXE.WT_vs_RBM20KO$WT_n5_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 4))
MXE.WT_vs_RBM20KO$WT_n6_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 5))
MXE.WT_vs_RBM20KO$RBM20KO_n3_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 1))
MXE.WT_vs_RBM20KO$RBM20KO_n4_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 2))
MXE.WT_vs_RBM20KO$RBM20KO_n5_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 3))
MXE.WT_vs_RBM20KO$RBM20KO_n6_PSI = as.character(lapply(MXE.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 4))
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                          RBM20KO_PSI, RBM20KO_n3_PSI, RBM20KO_n4_PSI, RBM20KO_n5_PSI, RBM20KO_n6_PSI, IncLevelDifference, FDR))
colnames(MXE.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                "RBM20KO_PSI", "RBM20KO_n3_PSI", "RBM20KO_n4_PSI", "RBM20KO_n5_PSI", "RBM20KO_n6_PSI", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

MXE.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_TTNKO/MXE.MATS.JC.txt")
MXE.WT_vs_TTNpromKO = MXE.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E1:", X1stExonStart_0base, "-", X1stExonEnd, "_E2:", X2ndExonStart_0base, "-", X2ndExonEnd, "_uE:", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
MXE.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(MXE.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
MXE.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(MXE.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
MXE.WT_vs_TTNpromKO$WT_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel1, get_PSI))
MXE.WT_vs_TTNpromKO$TTNpromKO_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel2, get_PSI))
MXE.WT_vs_TTNpromKO$WT_n1_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 1))
MXE.WT_vs_TTNpromKO$WT_n3_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 2))
MXE.WT_vs_TTNpromKO$WT_n4_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 3))
MXE.WT_vs_TTNpromKO$WT_n5_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 4))
MXE.WT_vs_TTNpromKO$WT_n6_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 5))
MXE.WT_vs_TTNpromKO$TTNpromKO_n1_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 1))
MXE.WT_vs_TTNpromKO$TTNpromKO_n3_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 2))
MXE.WT_vs_TTNpromKO$TTNpromKO_n4_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 3))
MXE.WT_vs_TTNpromKO$TTNpromKO_n5_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 4))
MXE.WT_vs_TTNpromKO$TTNpromKO_n6_PSI = as.character(lapply(MXE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 5))
MXE.WT_vs_TTNpromKO = MXE.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
MXE.WT_vs_TTNpromKO = MXE.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                              TTNpromKO_PSI, TTNpromKO_n1_PSI, TTNpromKO_n3_PSI, TTNpromKO_n4_PSI, TTNpromKO_n5_PSI, TTNpromKO_n6_PSI, IncLevelDifference, FDR))
colnames(MXE.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                  "TTNpromKO_PSI", "TTNpromKO_n1_PSI", "TTNpromKO_n3_PSI", "TTNpromKO_n4_PSI", "TTNpromKO_n5_PSI", "TTNpromKO_n6_PSI", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

MXE.WT_vs_doubleKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KOTTNKO/MXE.MATS.JC.txt") 
MXE.WT_vs_doubleKO = MXE.WT_vs_doubleKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E1:", X1stExonStart_0base, "-", X1stExonEnd, "_E2:", X2ndExonStart_0base, "-", X2ndExonEnd, "_uE:", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
MXE.WT_vs_doubleKO$JC_count_WT = as.character(lapply(MXE.WT_vs_doubleKO$IJC_SAMPLE_1, get_JCcounts))
MXE.WT_vs_doubleKO$JC_count_doubleKO = as.character(lapply(MXE.WT_vs_doubleKO$IJC_SAMPLE_2, get_JCcounts))
MXE.WT_vs_doubleKO$WT_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel1, get_PSI))
MXE.WT_vs_doubleKO$doubleKO_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel2, get_PSI))
MXE.WT_vs_doubleKO$WT_n1_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 1))
MXE.WT_vs_doubleKO$WT_n3_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 2))
MXE.WT_vs_doubleKO$WT_n4_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 3))
MXE.WT_vs_doubleKO$WT_n5_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 4))
MXE.WT_vs_doubleKO$WT_n6_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 5))
MXE.WT_vs_doubleKO$doubleKO_n1_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 1))
MXE.WT_vs_doubleKO$doubleKO_n3_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 2))
MXE.WT_vs_doubleKO$doubleKO_n4_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 3))
MXE.WT_vs_doubleKO$doubleKO_n5_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 4))
MXE.WT_vs_doubleKO$doubleKO_n6_PSI = as.character(lapply(MXE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 5))
MXE.WT_vs_doubleKO = MXE.WT_vs_doubleKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
MXE.WT_vs_doubleKO = MXE.WT_vs_doubleKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_doubleKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                            doubleKO_PSI, doubleKO_n1_PSI, doubleKO_n3_PSI, doubleKO_n4_PSI, doubleKO_n5_PSI, doubleKO_n6_PSI, IncLevelDifference, FDR))
colnames(MXE.WT_vs_doubleKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_doubleKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                 "doubleKO_PSI", "doubleKO_n1_PSI", "doubleKO_n3_PSI", "doubleKO_n4_PSI", "doubleKO_n5_PSI", "doubleKO_n6_PSI", "delta_PSI_doubleKOvsWT", "FDR_doubleKOvsWT")

MXE_intersect = MXE.WT_vs_RBM20KO %>% dplyr::inner_join(MXE.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT")) %>%
  dplyr::inner_join(MXE.WT_vs_doubleKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT"))
MXE_intersect = MXE_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_doubleKO)) %>% dplyr::filter(JC_sum > 19)

SE.WT_vs_RBM20KO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KO/SE.MATS.JC.txt") 
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E:", exonStart_0base, "-", exonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
SE.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(SE.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
SE.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(SE.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
SE.WT_vs_RBM20KO$WT_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel1, get_PSI))
SE.WT_vs_RBM20KO$RBM20KO_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel2, get_PSI))
SE.WT_vs_RBM20KO$WT_n1_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 1))
SE.WT_vs_RBM20KO$WT_n3_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 2))
SE.WT_vs_RBM20KO$WT_n4_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 3))
SE.WT_vs_RBM20KO$WT_n5_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 4))
SE.WT_vs_RBM20KO$WT_n6_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel1, get_PSI_n, n = 5))
SE.WT_vs_RBM20KO$RBM20KO_n3_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 1))
SE.WT_vs_RBM20KO$RBM20KO_n4_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 2))
SE.WT_vs_RBM20KO$RBM20KO_n5_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 3))
SE.WT_vs_RBM20KO$RBM20KO_n6_PSI = as.character(lapply(SE.WT_vs_RBM20KO$IncLevel2, get_PSI_n, n = 4))
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                        RBM20KO_PSI, RBM20KO_n3_PSI, RBM20KO_n4_PSI, RBM20KO_n5_PSI, RBM20KO_n6_PSI, IncLevelDifference, FDR))
colnames(SE.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                               "RBM20KO_PSI", "RBM20KO_n3_PSI", "RBM20KO_n4_PSI", "RBM20KO_n5_PSI", "RBM20KO_n6_PSI", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

SE.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_TTNKO/SE.MATS.JC.txt")
SE.WT_vs_TTNpromKO = SE.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E:", exonStart_0base, "-", exonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
SE.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(SE.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
SE.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(SE.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
SE.WT_vs_TTNpromKO$WT_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel1, get_PSI))
SE.WT_vs_TTNpromKO$TTNpromKO_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel2, get_PSI))
SE.WT_vs_TTNpromKO$WT_n1_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 1))
SE.WT_vs_TTNpromKO$WT_n3_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 2))
SE.WT_vs_TTNpromKO$WT_n4_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 3))
SE.WT_vs_TTNpromKO$WT_n5_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 4))
SE.WT_vs_TTNpromKO$WT_n6_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel1, get_PSI_n, n = 5))
SE.WT_vs_TTNpromKO$TTNpromKO_n1_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 1))
SE.WT_vs_TTNpromKO$TTNpromKO_n3_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 2))
SE.WT_vs_TTNpromKO$TTNpromKO_n4_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 3))
SE.WT_vs_TTNpromKO$TTNpromKO_n5_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 4))
SE.WT_vs_TTNpromKO$TTNpromKO_n6_PSI = as.character(lapply(SE.WT_vs_TTNpromKO$IncLevel2, get_PSI_n, n = 5))
SE.WT_vs_TTNpromKO = SE.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
SE.WT_vs_TTNpromKO = SE.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                            TTNpromKO_PSI, TTNpromKO_n1_PSI, TTNpromKO_n3_PSI, TTNpromKO_n4_PSI, TTNpromKO_n5_PSI, TTNpromKO_n6_PSI, IncLevelDifference, FDR))
colnames(SE.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                 "TTNpromKO_PSI", "TTNpromKO_n1_PSI", "TTNpromKO_n3_PSI", "TTNpromKO_n4_PSI", "TTNpromKO_n5_PSI", "TTNpromKO_n6_PSI", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

SE.WT_vs_doubleKO <- read.delim("/home/shared_folder/splicing_data/WTC11_rMATS_data/WT_vs_RBM20KOTTNKO/SE.MATS.JC.txt") 
SE.WT_vs_doubleKO = SE.WT_vs_doubleKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E:", exonStart_0base, "-", exonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
SE.WT_vs_doubleKO$JC_count_WT = as.character(lapply(SE.WT_vs_doubleKO$IJC_SAMPLE_1, get_JCcounts))
SE.WT_vs_doubleKO$JC_count_doubleKO = as.character(lapply(SE.WT_vs_doubleKO$IJC_SAMPLE_2, get_JCcounts))
SE.WT_vs_doubleKO$WT_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel1, get_PSI))
SE.WT_vs_doubleKO$doubleKO_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel2, get_PSI))
SE.WT_vs_doubleKO$WT_n1_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 1))
SE.WT_vs_doubleKO$WT_n3_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 2))
SE.WT_vs_doubleKO$WT_n4_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 3))
SE.WT_vs_doubleKO$WT_n5_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 4))
SE.WT_vs_doubleKO$WT_n6_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel1, get_PSI_n, n = 5))
SE.WT_vs_doubleKO$doubleKO_n1_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 1))
SE.WT_vs_doubleKO$doubleKO_n3_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 2))
SE.WT_vs_doubleKO$doubleKO_n4_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 3))
SE.WT_vs_doubleKO$doubleKO_n5_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 4))
SE.WT_vs_doubleKO$doubleKO_n6_PSI = as.character(lapply(SE.WT_vs_doubleKO$IncLevel2, get_PSI_n, n = 5))
SE.WT_vs_doubleKO = SE.WT_vs_doubleKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
SE.WT_vs_doubleKO = SE.WT_vs_doubleKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_doubleKO, WT_PSI, WT_n1_PSI, WT_n3_PSI, WT_n4_PSI, WT_n5_PSI, WT_n6_PSI,
                                                          doubleKO_PSI, doubleKO_n1_PSI, doubleKO_n3_PSI, doubleKO_n4_PSI, doubleKO_n5_PSI, doubleKO_n6_PSI, IncLevelDifference, FDR))
colnames(SE.WT_vs_doubleKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_doubleKO", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI",
                                "doubleKO_PSI", "doubleKO_n1_PSI", "doubleKO_n3_PSI", "doubleKO_n4_PSI", "doubleKO_n5_PSI", "doubleKO_n6_PSI", "delta_PSI_doubleKOvsWT", "FDR_doubleKOvsWT")

SE_intersect = SE.WT_vs_RBM20KO %>% dplyr::inner_join(SE.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT")) %>%
  dplyr::inner_join(SE.WT_vs_doubleKO, by = c("variant_ID", "geneSymbol", "geneID", "WT_PSI", "WT_n1_PSI", "WT_n3_PSI", "WT_n4_PSI", "WT_n5_PSI", "WT_n6_PSI", "JC_WT"))
SE_intersect = SE_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_doubleKO)) %>% dplyr::filter(JC_sum > 19)

merge = rbind(A3SS_intersect %>% dplyr::mutate(event = "A3SS"),
              A5SS_intersect %>% dplyr::mutate(event = "A5SS"),
              RI_intersect %>% dplyr::mutate(event = "RI"),
              MXE_intersect %>% dplyr::mutate(event = "MXE"),
              SE_intersect %>% dplyr::mutate(event = "SE"))

merge[is.na(merge)] = ""
merge = merge %>% dplyr::mutate(eCLIP = ifelse(geneSymbol %in% transC_lists$Gene, "RBM20_target", "non_target")) %>%
  dplyr::mutate(eCLIP = ifelse(geneSymbol == "TTN", "TTN", eCLIP)) %>%
  dplyr::mutate(WT_OMAP = ifelse(geneID %in% WT_OMAP$Gene, "pulldown", "non_target")) %>% 
  dplyr::mutate(RBM20KO_OMAP = ifelse(geneID %in% KO_OMAP$Gene, "pulldown", "non_target"))

merge = merge %>%
  dplyr::mutate(RBM20KO = ifelse(FDR_RBM20KOvsWT < 0.05 & delta_PSI_RBM20KOvsWT > 0.10, "upregulated", ifelse(FDR_RBM20KOvsWT < 0.05 & delta_PSI_RBM20KOvsWT < -0.10, "downregulated", "NS"))) %>%
  dplyr::mutate(TTNpromKO = ifelse(FDR_TTNpromKOvsWT < 0.05 & delta_PSI_TTNpromKOvsWT > 0.10, "upregulated", ifelse(FDR_TTNpromKOvsWT < 0.05 & delta_PSI_TTNpromKOvsWT < -0.10, "downregulated", "NS"))) %>%
  dplyr::mutate(doubleKO = ifelse(FDR_doubleKOvsWT < 0.05 & delta_PSI_doubleKOvsWT > 0.10, "upregulated", ifelse(FDR_doubleKOvsWT < 0.05 & delta_PSI_doubleKOvsWT < -0.10, "downregulated", "NS")))

#scatterplots splicing events - figure 5 panel H
merge_sign = merge %>% dplyr::filter(RBM20KO != "NS" | TTNpromKO != "NS" | doubleKO != "NS") %>% dplyr::arrange(doubleKO)
lin_all = lm(as.numeric((merge_sign)$delta_PSI_RBM20KOvsWT) ~ as.numeric((merge_sign)$delta_PSI_TTNpromKOvsWT))
ggplot(merge_sign, aes(x = delta_PSI_RBM20KOvsWT, y = delta_PSI_TTNpromKOvsWT)) + geom_point(aes(fill = event), shape = 21, size = 3, alpha = 0.5, color = "black") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_fill_manual(values = c("#E69F00", "#DEA0FD", "#ff6db6", "#009E73", "#56B4E9")) +
  scale_y_continuous(limits = c(-1, 1)) + scale_x_continuous(limits = c(-1, 1)) +
  ggtitle(paste0("R ^2 linear correlation = ", signif(summary(lin_all)$r.squared, digits = 3))) +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/scatterplot_spicing_events.pdf")

#scatterplots splicing events - supplementary figure 7 panel C
lin_eCLIP = lm(as.numeric((merge_sign %>% dplyr::filter(eCLIP == "RBM20_target" & geneSymbol != "TTN"))$delta_PSI_RBM20KOvsWT) ~ as.numeric((merge_sign %>% dplyr::filter(eCLIP == "RBM20_target" & geneSymbol != "TTN"))$delta_PSI_TTNpromKOvsWT))
merge_sign = merge_sign %>% dplyr::mutate(gene_name = ifelse(eCLIP == "RBM20_target", geneSymbol, ""))
ggplot(merge_sign %>% dplyr::arrange(desc(eCLIP)), aes(x = delta_PSI_RBM20KOvsWT, y = delta_PSI_TTNpromKOvsWT)) + geom_point(aes(fill = eCLIP), shape = 21, size = 3, alpha = 0.5, color = "black") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_fill_manual(values = c("white", "#1C8356", "#ff6db6")) + ggrepel::geom_text_repel(aes(label = gene_name), max.overlaps = 50,	size = 2, hjust = 1.2) +
  scale_y_continuous(limits = c(-1, 1)) + scale_x_continuous(limits = c(-1, 1)) +
  ggtitle(paste0("R ^2 linear correlation = ", signif(summary(lin_eCLIP)$r.squared, digits = 3))) +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/scatterplot_eCLIP.pdf")

lin_WT = lm(as.numeric((merge_sign %>% dplyr::filter(WT_OMAP == "pulldown"))$delta_PSI_RBM20KOvsWT) ~ as.numeric((merge_sign %>% dplyr::filter(WT_OMAP == "pulldown"))$delta_PSI_TTNpromKOvsWT))
merge_sign = merge_sign %>% dplyr::mutate(gene_name = ifelse(WT_OMAP == "pulldown", geneSymbol, ""))
ggplot(merge_sign %>% dplyr::arrange(WT_OMAP), aes(x = delta_PSI_RBM20KOvsWT, y = delta_PSI_TTNpromKOvsWT)) + geom_point(aes(fill = WT_OMAP), shape = 21, size = 3, alpha = 0.5, color = "black") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_fill_manual(values = c("white", "#3283FE")) + ggrepel::geom_text_repel(aes(label = gene_name), max.overlaps = 20,	size = 2, hjust = 1.2) +
  scale_y_continuous(limits = c(-1, 1)) + scale_x_continuous(limits = c(-1, 1)) +
  ggtitle(paste0("R ^2 linear correlation = ", signif(summary(lin_WT)$r.squared, digits = 3))) +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/scatterplot_OMAP_WT.pdf")

lin_KO = lm(as.numeric((merge_sign %>% dplyr::filter(RBM20KO_OMAP == "pulldown"))$delta_PSI_RBM20KOvsWT) ~ as.numeric((merge_sign %>% dplyr::filter(RBM20KO_OMAP == "pulldown"))$delta_PSI_TTNpromKOvsWT))
merge_sign = merge_sign %>% dplyr::mutate(gene_name = ifelse(RBM20KO_OMAP == "pulldown", geneSymbol, ""))
ggplot(merge_sign %>% dplyr::arrange(RBM20KO_OMAP), aes(x = delta_PSI_RBM20KOvsWT, y = delta_PSI_TTNpromKOvsWT)) + geom_point(aes(fill = RBM20KO_OMAP), shape = 21, size = 3, alpha = 0.5, color = "black") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_fill_manual(values = c("white", "#F6222E")) + ggrepel::geom_text_repel(aes(label = gene_name), max.overlaps = 20,	size = 2, hjust = 1.2) +
  scale_y_continuous(limits = c(-1, 1)) + scale_x_continuous(limits = c(-1, 1)) +
  ggtitle(paste0("R ^2 linear correlation = ", signif(summary(lin_KO)$r.squared, digits = 3))) +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/scatterplot_OMAP_RBM20KO.pdf")

#heatmap of splicing events with double KO effect - figure 5 panel I
RBM20KO_diff = merge %>% dplyr::filter(FDR_RBM20KOvsWT < 0.05 & abs(delta_PSI_RBM20KOvsWT) > 0.10) %>%
  dplyr::arrange(eCLIP) %>% dplyr::arrange(WT_OMAP) %>% dplyr::arrange(RBM20KO_OMAP) %>%
  dplyr::arrange(TTNpromKO) %>% dplyr::arrange(doubleKO) %>% dplyr::arrange(RBM20KO) %>% dplyr::arrange(event) 
rownames(RBM20KO_diff) = RBM20KO_diff$variant_ID
matrix = RBM20KO_diff %>% dplyr::select(c(WT_PSI, RBM20KO_PSI, TTNpromKO_PSI, doubleKO_PSI))
matrix = as.data.frame(lapply(matrix, as.numeric))
annotation_row = RBM20KO_diff %>% dplyr::select(event, doubleKO, TTNpromKO, RBM20KO, eCLIP, WT_OMAP, RBM20KO_OMAP)
rownames(matrix) = rownames(annotation_row)
annotation_colors = list(eCLIP = c(RBM20_target = "black", TTN = "black", non_target = "white"),
                         WT_OMAP = c(pulldown = "black", non_target = "white"),
                         RBM20KO_OMAP = c(pulldown = "black", non_target = "white"),
                         RBM20KO = c(upregulated = "#F6222E", downregulated = "#3283FE", NS = "grey"), 
                         TTNpromKO = c(upregulated = "#F6222E", downregulated = "#3283FE", NS = "grey"),
                         doubleKO = c(upregulated = "#F6222E", downregulated = "#3283FE", NS = "grey"),
                         event = c(A3SS = "#E69F00", A5SS = "#DEA0FD", RI = "#009E73", MXE = "#ff6db6", SE = "#56B4E9"))
gaps = c()
for(i in 2:nrow(RBM20KO_diff)){
  if(RBM20KO_diff$event[i] != RBM20KO_diff$event[i - 1]){gaps = append(gaps, i - 1)}
}
as.ggplot(pheatmap(matrix, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_row = annotation_row,
                   scale = "row", annotation_colors = annotation_colors, gaps_row = gaps, color = viridis::viridis(n = 25, option = "H")))
ggsave("/home/shared_folder/output_plots/heatmap_RBM20KO_sign_PSI.pdf")

#GSEA of O-MAP enrichment of recovered splicing variants - supplementary figure 7 panel B
delta_OMAP <- read_excel("/home/shared_folder/pulldown_data/transcripts_data/o_map_WT_KO.xlsx")
colnames(delta_OMAP) = c("Gene", "Symbol", "logFC_WTvsKO", "padj")
rank = delta_OMAP$logFC_WTvsKO
names(rank) = delta_OMAP$Gene
genes_list = RBM20KO_diff %>% dplyr::filter(doubleKO == "NS")

genes_A3SS = (genes_list %>% dplyr::filter(event == "A3SS"))$geneID
genes_A5SS = (genes_list %>% dplyr::filter(event == "A5SS"))$geneID
genes_RI = (genes_list %>% dplyr::filter(event == "RI"))$geneID
genes_MXE = (genes_list %>% dplyr::filter(event == "MXE"))$geneID
genes_SE = (genes_list %>% dplyr::filter(event == "SE"))$geneID

pathway = list(A3SS = genes_A3SS, A5SS = genes_A5SS, RI = genes_RI, MXE = genes_MXE, SE = genes_SE)

fgsea <- fgsea(pathways = pathway, 
               stats    = rank,
               scoreType = 'std',
               minSize  = 0,
               maxSize  = Inf)
p = plotEnrichment(pathway[["A3SS"]], rank) + 
  labs(title = paste0("A3SS ranking enrichment\nNES = ", fgsea$NES[1], "\np-val = ", fgsea$pval[1]))
p[["layers"]][[1]][["aes_params"]]$colour = "#E69F00"
as.ggplot(print(p))
ggsave("/home/shared_folder/output_plots/A3SS_gsea.pdf")

p = plotEnrichment(pathway[["A5SS"]], rank) + 
  labs(title = paste0("A5SS ranking enrichment\nNES = ", fgsea$NES[2], "\np-val = ", fgsea$pval[2]))
p[["layers"]][[1]][["aes_params"]]$colour = "#DEA0FD"
as.ggplot(print(p))
ggsave("/home/shared_folder/panelE/A5SS_gsea.pdf")

p = plotEnrichment(pathway[["RI"]], rank) + 
  labs(title = paste0("RI ranking enrichment\nNES = ", fgsea$NES[4], "\np-val = ", fgsea$pval[4]))
p[["layers"]][[1]][["aes_params"]]$colour = "#009E73"
as.ggplot(print(p))
ggsave("/home/shared_folder/output_plots/RI_gsea.pdf")

p = plotEnrichment(pathway[["MXE"]], rank) + 
  labs(title = paste0("MXE ranking enrichment\nNES = ", fgsea$NES[3], "\np-val = ", fgsea$pval[3]))
p[["layers"]][[1]][["aes_params"]]$colour = "#ff6db6"
as.ggplot(print(p))
ggsave("/home/shared_folder/output_plots/MXE_gsea.pdf")

p = plotEnrichment(pathway[["SE"]], rank) + 
  labs(title = paste0("SE ranking enrichment\nNES = ", fgsea$NES[5], "\np-val = ", fgsea$pval[5]))
p[["layers"]][[1]][["aes_params"]]$colour = "#56B4E9"
as.ggplot(print(p))
ggsave("/home/shared_folder/output_plots/SE_gsea.pdf")

#linear correlation of double KOs with signle KOs - supplementary figure 10 panel A
R_filtered = data.frame(event = c(), TTNpromKO_effect = c(), RBM20KO_effect = c())

lin1 = lm(as.numeric((A3SS_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((A3SS_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$RBM20KO_PSI))
lin2 = lm(as.numeric((A3SS_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((A3SS_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$TTNpromKO_PSI))
R_filtered = rbind(R_filtered, data.frame(event = c("A3SS"), TTNpromKO_effect = c(summary(lin1)$r.squared), RBM20KO_effect = c(summary(lin2)$r.squared)))

lin1 = lm(as.numeric((A5SS_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((A5SS_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$RBM20KO_PSI))
lin2 = lm(as.numeric((A5SS_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((A5SS_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$TTNpromKO_PSI))
R_filtered = rbind(R_filtered, data.frame(event = c("A5SS"), TTNpromKO_effect = c(summary(lin1)$r.squared), RBM20KO_effect = c(summary(lin2)$r.squared)))

lin1 = lm(as.numeric((RI_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((RI_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$RBM20KO_PSI))
lin2 = lm(as.numeric((RI_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((RI_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$TTNpromKO_PSI))
R_filtered = rbind(R_filtered, data.frame(event = c("RI"), TTNpromKO_effect = c(summary(lin1)$r.squared), RBM20KO_effect = c(summary(lin2)$r.squared)))

lin1 = lm(as.numeric((MXE_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((MXE_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$RBM20KO_PSI))
lin2 = lm(as.numeric((MXE_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((MXE_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$TTNpromKO_PSI))
R_filtered = rbind(R_filtered, data.frame(event = c("MXE"), TTNpromKO_effect = c(summary(lin1)$r.squared), RBM20KO_effect = c(summary(lin2)$r.squared)))

lin1 = lm(as.numeric((SE_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((SE_intersect %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05))$RBM20KO_PSI))
lin2 = lm(as.numeric((SE_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$doubleKO_PSI) ~ as.numeric((SE_intersect %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05))$TTNpromKO_PSI))
R_filtered = rbind(R_filtered, data.frame(event = c("SE"), TTNpromKO_effect = c(summary(lin1)$r.squared), RBM20KO_effect = c(summary(lin2)$r.squared)))

melt_R_filtered = melt(R_filtered)
colnames(melt_R_filtered) = c("event", "effect", "R_squared")
ggplot(melt_R_filtered) + geom_col(aes(x = event, y = R_squared, fill = effect), position = "dodge") + scale_fill_manual(values = c("#56B4E9", "#D55E00"))
ggsave("/home/shared_folder/output_plots/summary_effect.pdf")

breaks_log2 = 2^seq(-10, 0, by = 1)
ggplot(merge %>% dplyr::filter(FDR_RBM20KOvsWT > 0.05)) +
  geom_density_2d_filled(aes(x = as.numeric(TTNpromKO_PSI), y = as.numeric(doubleKO_PSI)), breaks = breaks_log2, contour_var = "ndensity") +
  theme(legend.position = "none") + geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed")
ggsave("/home/shared_folder/output_plots/RBM20KO_undirect_effect_density.pdf")

ggplot(merge %>% dplyr::filter(FDR_TTNpromKOvsWT > 0.05)) +
  geom_density_2d_filled(aes(x = as.numeric(RBM20KO_PSI), y = as.numeric(doubleKO_PSI)), breaks = breaks_log2, contour_var = "ndensity") +
  theme(legend.position = "none") + geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed")
ggsave("/home/shared_folder/output_plots/TTNpromKO_undirect_effect_density.pdf")

#motifs enrichment in rescued transcripts
mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
transcripts_all = geneInfo$Symbol
background = getSequence(id = c(transcripts_all), 
                         type = "hgnc_symbol", 
                         seqType = "transcript_exon_intron", 
                         mart = mart)
background$keep = as.character(lapply(background$transcript_exon_intron, check_bases))
background = background %>% dplyr::filter(keep == 0)
background = background %>% dplyr::mutate(name = paste0(hgnc_symbol, "|mRNA"))
len = c()
for(i in 1:nrow(background)){
  len = append(len, nchar(background$transcript_exon_intron[i]))
}
background$len = len
background$n = 1:nrow(background)
stat = background %>% dplyr::group_by(hgnc_symbol) %>% dplyr::summarise(len = max(len)) %>% dplyr::mutate(k = "yes")
background_filtered = background %>% dplyr::left_join(stat, by = c("hgnc_symbol", "len")) %>% dplyr::filter(k == "yes") %>% dplyr::distinct()
transcripts_background = background_filtered$transcript_exon_intron
names(transcripts_background) = background_filtered$name

transcripts_list = unique(c((genes_list %>% dplyr::filter(event == "A3SS"))$geneSymbol,
                            (genes_list %>% dplyr::filter(event == "A5SS"))$geneSymbol,
                            (genes_list %>% dplyr::filter(event == "RI"))$geneSymbol,
                            (genes_list %>% dplyr::filter(event == "MXE"))$geneSymbol,
                            (genes_list %>% dplyr::filter(event == "SE"))$geneSymbol))
foreground = getSequence(id = c(transcripts_list), 
                         type = "hgnc_symbol", 
                         seqType = "transcript_exon_intron", 
                         mart = mart)
foreground$keep = as.character(lapply(foreground$transcript_exon_intron, check_bases))
foreground = foreground %>% dplyr::filter(keep == 0)
foreground = foreground %>% dplyr::mutate(name = paste0(hgnc_symbol, "|mRNA"))
len = c()
for(i in 1:nrow(foreground)){
  len = append(len, nchar(foreground$transcript_exon_intron[i]))
}
foreground$len = len
foreground$n = 1:nrow(foreground)
stat = foreground %>% dplyr::group_by(hgnc_symbol) %>% dplyr::summarise(len = max(len)) %>% dplyr::mutate(k = "yes")
foreground_filtered = foreground %>% dplyr::left_join(stat, by = c("hgnc_symbol", "len")) %>% dplyr::filter(k == "yes") %>% dplyr::distinct()
transcripts_foreground = foreground_filtered$transcript_exon_intron
names(transcripts_foreground) = foreground_filtered$name
foreground_set = list(transcripts_foreground)

RBP_Information_all_motifs <- read.delim("/home/shared_folder/RBP_data/RBP_Information_all_motifs.txt")
motifs = list()
for(i in 1:nrow(RBP_Information_all_motifs)){
  if(RBP_Information_all_motifs$Motif_ID[i] != "" & file.exists(paste0("/home/shared_folder/RBP_data/pwms_all_motifs/pwms_all_motifs/", RBP_Information_all_motifs$Motif_ID[i], ".txt"))){
    motif_matrix <- read.delim(paste0("/home/shared_folder/input_data/pwms_all_motifs/pwms_all_motifs/", RBP_Information_all_motifs$Motif_ID[i], ".txt")) %>% dplyr::select(-Pos)
    motif = create_matrix_motif(paste0("id", i), RBP_Information_all_motifs$RBP_Name[i], motif_matrix, "RBPDB", "hs", "RBDP_database")
    motifs = append(motifs, motif)
  }
}

df_enrichment = data.frame(foreground_abs = c(), background_abs = c(),
                           foreground_rel = c(), background_rel = c(),
                           ratio = c(), pval = c(), rbp = c())

motifs_short = list()
for(i in 1:length(motifs)){
  if(motifs[[i]]@length <= 20){motifs_short = append(motifs_short, motifs[i])}
}

len = length(motifs)
motifs_long = list()
for(i in 1:length(motifs_short)){
  enrichment_i = run_matrix_tsma(
    foreground_set,
    transcripts_background,
    motifs = motifs[i],
    max_hits = Inf,
    threshold_method = "p_value",
    threshold_value = 0.05,
    n_cores = 20,
    cache = FALSE)
  df_enrichment = rbind(df_enrichment,
                        data.frame(foreground_abs = c(enrichment_i[["foreground_scores"]][[1]][["df"]][["absolute_hits"]]),
                                   background_abs = c(enrichment_i[["background_scores"]][["df"]][["absolute_hits"]]),
                                   foreground_rel = c(enrichment_i[["foreground_scores"]][[1]][["df"]][["relative_hits"]]),
                                   background_rel = c(enrichment_i[["background_scores"]][["df"]][["relative_hits"]]),
                                   ratio = c(enrichment_i[["enrichment_dfs"]][[1]][["enrichment"]]),
                                   pval = c(enrichment_i[["enrichment_dfs"]][[1]][["p_value"]]),
                                   rbp = c(enrichment_i[["enrichment_dfs"]][[1]][["motif_rbps"]])))
}

proteomics_pulldown <- read_excel("/home/shared_folder/pulldown_data/proteomics_data/proteomics_pulldown.xlsx", skip = 2)[, 1:19]
proteomics_pulldown[is.na(proteomics_pulldown)] = ""
pulldown_WT = (proteomics_pulldown %>% dplyr::filter(T_significance == "+" & Log2FC_KO_WT < 0))$Gene
pulldown_KO = (proteomics_pulldown %>% dplyr::filter(T_significance == "+" & Log2FC_KO_WT > 0))$Gene
enrichment$qval = p.adjust(enrichment$pval, method = "BH")
enrichment = enrichment %>% dplyr::mutate(relative_enrichment = ratio - 1) %>%
  dplyr::mutate(rbp_name = ifelse(qval < 0.05 & abs(relative_enrichment) > 0.02, rbp, ""))
enrichment = enrichment %>% dplyr::mutate(pulldown = ifelse(qval > 0.05 | abs(relative_enrichment) < 0.02, "insignificant", ifelse(rbp %in% pulldown_WT, "WT", ifelse(rbp %in% pulldown_KO, "RBM20KO", ""))))


#tomtom tpo evaluate motifs similarity
dir.create("/home/shared_folder/memes")
system2("chmod", "777 -R /home/shared_folder/memes")
for(i in 1:nrow(RBP_Information_all_motifs)){
  if(RBP_Information_all_motifs$Motif_ID[i] != "" & file.exists(paste0("/home/shared_folder/RBP_data/pwms_all_motifs/pwms_all_motifs/", RBP_Information_all_motifs$Motif_ID[i], ".txt"))){
    motif_matrix <- read.delim(paste0("/home/shared_folder/RBP_data/pwms_all_motifs/pwms_all_motifs/", RBP_Information_all_motifs$Motif_ID[i], ".txt")) %>% dplyr::select(-Pos)
    system2("cp", paste0("/home/shared_folder/RBP_data/template.meme /home/shared_folder/memes/", RBP_Information_all_motifs$Motif_ID[i], ".meme"))
    system2("sed", paste0("-i 's|w= 19|w= ", nrow(motif_matrix), "|g' ", "/home/shared_folder/memes/", RBP_Information_all_motifs$Motif_ID[i], ".meme"))
    
    for(j in 1:nrow(motif_matrix)){
      string = paste0(motif_matrix$A[j], " ", motif_matrix$C[j], " " , motif_matrix$G[j], " ", motif_matrix$U[j])
      write(string, file = paste0("/home/shared_folder/memes/", RBP_Information_all_motifs$Motif_ID[i], ".meme"), append = TRUE)
    }
  }
}

ids_valid = c()
length = c()
shannon = c()
for(i in 1:nrow(RBP_Information_all_motifs)){
  if(RBP_Information_all_motifs$Motif_ID[i] != "" & file.exists(paste0("/home/shared_folder/RBP_data/pwms_all_motifs/pwms_all_motifs/", RBP_Information_all_motifs$Motif_ID[i], ".txt"))){
    ids_valid = append(ids_valid, RBP_Information_all_motifs$Motif_ID[i])
    motif_matrix <- read.delim(paste0("/home/shared_folder/RBP_data/pwms_all_motifs/pwms_all_motifs/", RBP_Information_all_motifs$Motif_ID[i], ".txt")) %>% dplyr::select(-Pos)
    entropy = 0
    for(j in 1:nrow(motif_matrix)){
      entropy = entropy + 2
      for(k in 1:ncol(motif_matrix)){
        if(motif_matrix[j, k] != 0){entropy = entropy + motif_matrix[j, k]*log2(motif_matrix[j, k])}
      }
    }
    length = append(length, nrow(motif_matrix))
    shannon = append(shannon, entropy/nrow(motif_matrix))
  }
}

RBPs_filtered = RBP_Information_all_motifs %>% dplyr::filter(Motif_ID %in% ids_valid) %>% dplyr::select(c(Motif_ID, RBP_Name))
RBPs_filtered$motif_length = length
RBPs_filtered$shannon_entropy = shannon
RBPs_filtered = RBPs_filtered %>% dplyr::filter(length <= 20)

RBPs_filtered = RBPs_filtered %>% dplyr::distinct()
tomtom_all = data.frame(Query_ID = c(), Target_ID = c(), Optimal_offset = c(), p.value = c(), E.value = c(), q.value = c(), Overlap = c(), Query_consensus = c(), Target_consensus = c(), Orientation = c(), RBP = c())
for(i in unique(RBPs_filtered$RBP_Name)){
  print(i)
  RBP_i = RBPs_filtered %>% dplyr::filter(RBP_Name == i)
  for(j in 1:nrow(RBP_i)){
    for(k in j:nrow(RBP_i)){
      if(j != k){
        system2("/usr/local/anaconda/bin/tomtom", paste0("/home/shared_folder/memes/", RBP_i$Motif_ID[j], ".meme /home/shared_folder/memes/",RBP_i$Motif_ID[k], ".meme -oc /home/shared_folder/tomtom"))
        tomtom <- read.delim("/home/shared_folder/tomtom/tomtom.tsv", comment.char="#")
        if(nrow(tomtom)!= 0){
          tomtom$Query_ID = RBP_i$Motif_ID[j]
          tomtom$Target_ID = RBP_i$Motif_ID[k]
          tomtom$RBP = i
          tomtom_all = rbind(tomtom_all, tomtom)
        }
      }
    }
  }
}


motifs = data.frame(Motif_ID = unique(RBPs_filtered$Motif_ID))
overlaps = data.frame(Motif_ID = unique(RBPs_filtered$Motif_ID))
meta = data.frame(Motif_ID = c(), proteins = c())
for(i in 1:nrow(motifs)){
  proteins_i = unique((RBPs_filtered %>% dplyr::filter(Motif_ID == motifs$Motif_ID[i]))$RBP_Name)
  proteins_list = ""
  for(j in 1:length(proteins_i)){
    if(j != 1){proteins_list = paste0(proteins_list, "_", proteins_i[j])}
    else{proteins_list = paste0(proteins_list, proteins_i[j])}
  }
  meta = rbind(meta, data.frame(Motif_ID = c(motifs$Motif_ID[i]), proteins = c(proteins_list)))
  motifs_j = data.frame(Motif_ID = c(), pval = c())
  overlaps_j = data.frame(Motif_ID = c(), overlaps = c())
  for(j in (i+1):nrow(motifs)){
    system2("/usr/local/anaconda/bin/tomtom", paste0("/home/shared_folder/memes/", motifs$Motif_ID[i], ".meme /home/shared_folder/memes/", motifs$Motif_ID[j], ".meme -oc /home/shared_folder/panelG/tomtom"))
    tomtom <- read.delim("/home/shared_folder/tomtom/tomtom.tsv", comment.char="#")
    if(nrow(tomtom) != 0){
        motifs_j = rbind(motifs_j, data.frame(Motif_ID = c(motifs$Motif_ID[j]), pval = tomtom$p.value[1]))
        overlaps_j = rbind(overlaps_j, data.frame(Motif_ID = c(motifs$Motif_ID[j]), pval = tomtom$Overlap[1]))
    }
    else{
      motifs_j = rbind(motifs_j, data.frame(Motif_ID = c(motifs$Motif_ID[j]), pval = NA))
      overlaps_j = rbind(overlaps_j, data.frame(Motif_ID = c(motifs$Motif_ID[j]), pval = NA))
    }
  }
  motifs = motifs %>% dplyr::left_join(motifs_j, by = "Motif_ID")
  overlaps = overlaps %>% dplyr::left_join(overlaps_j, by = "Motif_ID")
}
rownames(motifs) = motifs$Motif_ID
motifs = motifs %>% dplyr::select(-Motif_ID)
colnames(motifs) = c(unique(RBPs_filtered$Motif_ID))
rownames(overlaps) = overlaps$Motif_ID
overlaps = overlaps %>% dplyr::select(-Motif_ID)
colnames(overlaps) = c(unique(RBPs_filtered$Motif_ID))

motifs[is.na(motifs)] = 1
total = unique(RBPs_filtered$RBP_Name)
results = data.frame(pval_threshold = c(), dup_fraction = c(), n_clusters = c())
for(pval in c(0.05, 0.01, 0.001, 0.0001, 0.00001)){
  motifs_filtering = motifs
  meta_filtered = meta
  motifs_keep = c()
  motifs_discard = c()
  while(any(unlist(motifs_filtering) < pval)){
    for(i in 1:nrow(motifs_filtering)){
      for(j in 1:ncol(motifs_filtering)){
        restart = FALSE
        pval_i = motifs_filtering[i, j]
        if(pval_i < pval){
          motif_i = rownames(motifs_filtering)[i]
          motif_j = colnames(motifs_filtering)[j]
          shannon_i = as.numeric((RBPs_filtered %>% dplyr::filter(Motif_ID == motif_i))$shannon_entropy[1])
          shannon_j = as.numeric((RBPs_filtered %>% dplyr::filter(Motif_ID == motif_j))$shannon_entropy[1])
          if(shannon_i > shannon_j){
            motif_discard = motif_j
            motif_keep = motif_i
          } else{
            motif_discard = motif_i
            motif_keep = motif_j
          }
          motifs_filtering = motifs_filtering[, (colnames(motifs_filtering) != motif_discard)]
          motifs_filtering = motifs_filtering[(rownames(motifs_filtering) != motif_discard), ]
          print(dim(motifs_filtering))
          meta_keep = meta_filtered %>% dplyr::filter(Motif_ID == motif_keep)
          print(meta_keep$proteins[1])
          print(strsplit(meta_keep$proteins[1], "_")[[1]])
          meta_discard = meta_filtered %>% dplyr::filter(Motif_ID == motif_discard)
          print(meta_discard$proteins[1])
          print(strsplit(meta_discard$proteins[1], "_")[[1]])
          proteins = sort(unique(c(strsplit(meta_keep$proteins[1], "_")[[1]], strsplit(meta_discard$proteins[1], "_")[[1]])))
          proteins_list = ""
          for(k in 1:length(proteins)){
            if(k != 1){proteins_list = paste0(proteins_list, "_", proteins[k])}
            else{proteins_list = paste0(proteins_list, proteins[k])}
          }
          print(proteins_list)
          meta_filtered = meta_filtered %>% dplyr::filter(!Motif_ID %in% c(motif_discard, motif_keep))
          meta_filtered = rbind(meta_filtered, data.frame(Motif_ID = c(motif_keep), proteins = c(proteins_list)))
          print(nrow(meta_filtered))
          print(length(unique(c(rownames(motifs_filtering), colnames(motifs_filtering)))))
          restart = TRUE
        }
        if(restart){break}
      }
      if(restart){break}
    }
  }
  
  enrichment_perMotif = meta_filtered %>% dplyr::left_join(enrichment, by = "Motif_ID") %>% dplyr::group_by(Motif_ID) %>% dplyr::summarise(relative_enrichment = mean(relative_enrichment), pval = mean(pval)) %>%
    dplyr::distinct() %>% dplyr::left_join(meta_filtered, by = "Motif_ID")
  enrichment_perMotif$qval = p.adjust(enrichment_perMotif$pval, method = "BH")
  enrichment_perMotif = enrichment_perMotif %>% dplyr::mutate(color = ifelse(qval < 0.05 & relative_enrichment > 0, "enriched", ifelse(qval < 0.05 & relative_enrichment < 0, "depleted", "insignificant"))) %>%
    dplyr::mutate(rbp_name = ifelse(color != "insignificant", proteins, ""))
  enrichment_perMotif = enrichment_perMotif %>% dplyr::left_join(RBPs_filtered %>% dplyr::select(-RBP_Name) %>% dplyr::distinct(), by = "Motif_ID")
  ggplot(enrichment_perMotif %>% dplyr::arrange(pval), aes(x = relative_enrichment, y = -log10(qval))) + geom_point(aes(color = color, fill = color), size = 4, shape = 21) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") + ggrepel::geom_text_repel(aes(label = rbp_name), max.overlaps = 40,	size = 3, hjust = 1.2) +
    scale_fill_manual(values = c("#56B4E9", "#FEAF16", "white")) + scale_color_manual(values = c("#0072B2", "#D55E00", "grey"))
  ggsave(paste0("/home/shared_folder/panelG/motifs_aggregate_volcano_threhsold_pval", pval, ".pdf"), width = 12, height = 8)
  ggplot(enrichment_perMotif) + geom_density(aes(x = shannon_entropy), fill = "#B6DBFF") + xlab("information_content") + scale_x_continuous(limits = c(0, 2))
  ggsave(paste0("/home/shared_folder/panelG/information_density_threhsold_pval", pval, ".pdf"))
  write.csv(enrichment_perMotif, paste0("/home/shared_folder/output_plots/motifs_aggregates_threshold_pval", pval, ".csv"))
  
  n_motifs = nrow(enrichment_perMotif)
  proteins = c()
  for(i in enrichment_perMotif$proteins){proteins = append(proteins, unique(strsplit(i, "_")[[1]]))}
  proteins = na.omit(proteins)
  results = rbind(results, data.frame(pval_threshold = c(pval), dup_fraction = c(length(proteins)/length(total)), n_clusters = c(n_motifs)))

}

ggplot(results, aes(x = dup_fraction, y = n_clusters)) + geom_line() + geom_point(size = 3, fill = "#3283FE", shape = 21) + geom_label(aes(label = pval_threshold))
ggsave("/home/shared_folder/output_plots/duplication_fraction.pdf")
