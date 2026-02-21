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
library(gridExtra)
library(readxl)
library(GO.db)
library(AnnotationDbi)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0))))

##tracks for KO validation - supplementary figure 8 panel B
rawcounts_filtered <- read.csv("/home/shared_folder/gene_expression_data/DGE_RUES2_data/rawcounts_filtered.csv") %>% dplyr::select(-c(X, gene_ID, gene_name, gene_type, gene_biotype, chromosome, chr_start, gc_content))
sums = as.data.frame(t(colSums(rawcounts_filtered)))
bed <- read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_E8_TTN.bed", header=FALSE) %>% dplyr::select(-V4) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_mTeSR1.1_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_mTeSR1.3_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_mTeSR1.4_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_mTeSR1.5_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3"))
bed[is.na(bed)] = 0
bed = bed %>% dplyr::mutate(chr = paste0("chr", V1)) %>% dplyr::mutate(value = (V5.x/sums$WT_E8[1] + V5.y/sums$WT_mTeSR1.1[1] + V5.x.x/sums$WT_mTeSR1.3[1] + V5.y.y/sums$WT_mTeSR1.4[1] + V5/sums$WT_mTeSR1.5[1])/5*100000000) %>% dplyr::select(c(chr, V2, V3, value))
write.table(bed, "/home/shared_folder/gene_expression_data/bedgraph_RUES2/WT_TTN.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

bed <- read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/RBM20.KO_E8_TTN.bed", header=FALSE) %>% dplyr::select(-V4) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/RBM20.KO_mTeSR1.1_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/RBM20.KO_mTeSR1.3_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3"))
bed[is.na(bed)] = 0
bed = bed %>% dplyr::mutate(chr = paste0("chr", V1)) %>% dplyr::mutate(value = (V5.x/sums$RBM20.KO_E8[1] + V5.y/sums$RBM20.KO_mTeSR1.1[1] + V5/sums$RBM20.KO_mTeSR1.3[1])/3*100000000) %>% dplyr::select(c(chr, V2, V3, value))
write.table(bed, "/home/shared_folder/gene_expression_data/bedgraph_RUES2/RBM20.KO_TTN.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

bed <- read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_E8_TTN.bed", header=FALSE) %>% dplyr::select(-V4) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_mTeSR1.3_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_mTeSR1.4_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_mTeSR1.5_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_mTeSR1.6_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3"))
bed[is.na(bed)] = 0
bed = bed %>% dplyr::mutate(chr = paste0("chr", V1)) %>% dplyr::mutate(value = (V5.x/sums$TTN.PROM.KO_E8[1] + V5.y/sums$TTN.PROM.KO_mTeSR1.3[1] + V5.x.x/sums$TTN.PROM.KO_mTeSR1.4[1] + V5.y.y/sums$TTN.PROM.KO_mTeSR1.5[1] + V5/sums$TTN.PROM.KO_mTeSR1.6[1])/5*100000000) %>% dplyr::select(c(chr, V2, V3, value))
write.table(bed, "/home/shared_folder/gene_expression_data/bedgraph_RUES2/TTN.PROM.KO_TTN.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

bed <- read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.KO100_E8_TTN.bed", header=FALSE) %>% dplyr::select(-V4) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.KO100_mTeSR1.3_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.KO115_mTeSR1.3_TTN.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3"))
bed[is.na(bed)] = 0
bed = bed %>% dplyr::mutate(chr = paste0("chr", V1)) %>% dplyr::mutate(value = (V5.x/sums$TTN.KO100_E8[1] + V5.y/sums$TTN.KO100_mTeSR1.3[1] + V5/sums$TTN.KO115_mTeSR1.3[1])/3*100000000) %>% dplyr::select(c(chr, V2, V3, value))
write.table(bed, "/home/shared_folder/gene_expression_data/bedgraph_RUES2/TTN.KO_TTN.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

bed <- read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_E8_RBM20.bed", header=FALSE) %>% dplyr::select(-V4) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_mTeSR1.1_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_mTeSR1.3_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_mTeSR1.4_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/WT_mTeSR1.5_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3"))
bed[is.na(bed)] = 0
bed = bed %>% dplyr::mutate(chr = paste0("chr", V1)) %>% dplyr::mutate(value = (V5.x/sums$WT_E8[1] + V5.y/sums$WT_mTeSR1.1[1] + V5.x.x/sums$WT_mTeSR1.3[1] + V5.y.y/sums$WT_mTeSR1.4[1] + V5/sums$WT_mTeSR1.5[1])/5*100000000) %>% dplyr::select(c(chr, V2, V3, value))
write.table(bed, "/home/shared_folder/gene_expression_data/bedgraph_RUES2/WT_RBM20.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

bed <- read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/RBM20.KO_E8_RBM20.bed", header=FALSE) %>% dplyr::select(-V4) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/RBM20.KO_mTeSR1.1_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/RBM20.KO_mTeSR1.3_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3"))
bed[is.na(bed)] = 0
bed = bed %>% dplyr::mutate(chr = paste0("chr", V1)) %>% dplyr::mutate(value = (V5.x/sums$RBM20.KO_E8[1] + V5.y/sums$RBM20.KO_mTeSR1.1[1] + V5/sums$RBM20.KO_mTeSR1.3[1])/3*100000000) %>% dplyr::select(c(chr, V2, V3, value))
write.table(bed, "/home/shared_folder/gene_expression_data/bedgraph_RUES2/RBM20.KO_RBM20.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

bed <- read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_E8_TTN.bed", header=FALSE) %>% dplyr::select(-V4) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_mTeSR1.3_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_mTeSR1.4_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_mTeSR1.5_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.PROM.KO_mTeSR1.6_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3"))
bed[is.na(bed)] = 0
bed = bed %>% dplyr::mutate(chr = paste0("chr", V1)) %>% dplyr::mutate(value = (V5.x/sums$TTN.PROM.KO_E8[1] + V5.y/sums$TTN.PROM.KO_mTeSR1.3[1] + V5.x.x/sums$TTN.PROM.KO_mTeSR1.4[1] + V5.y.y/sums$TTN.PROM.KO_mTeSR1.5[1] + V5/sums$TTN.PROM.KO_mTeSR1.6[1])/5*100000000) %>% dplyr::select(c(chr, V2, V3, value))
write.table(bed, "/home/shared_folder/gene_expression_data/bedgraph_RUES2/TTN.PROM.KO_RBM20.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

bed <- read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.KO100_E8_RBM20.bed", header=FALSE) %>% dplyr::select(-V4) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.KO100_mTeSR1.3_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3")) %>%
  dplyr::full_join(read.delim("/home/shared_folder/gene_expression_data/bed_RUES2/TTN.KO115_mTeSR1.3_RBM20.bed", header=FALSE) %>% dplyr::select(-V4), by = c("V1", "V2", "V3"))
bed[is.na(bed)] = 0
bed = bed %>% dplyr::mutate(chr = paste0("chr", V1)) %>% dplyr::mutate(value = (V5.x/sums$TTN.KO100_E8[1] + V5.y/sums$TTN.KO100_mTeSR1.3[1] + V5/sums$TTN.KO115_mTeSR1.3[1])/3*100000000) %>% dplyr::select(c(chr, V2, V3, value))
write.table(bed, "/home/shared_folder/gene_expression_data/bedgraph_RUES2/TTN.KO_RBM20.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

##differential splicing analysis
#functions
get_JCcounts <- function(str){
  val = strsplit(str, ",")[[1]]
  sum = 0
  for(i in 1:length(val)){
    sum = sum + as.numeric(val[i])
  }
  return(sum)
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

A3SS.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/A3SS.MATS.JC.txt"))
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A3SS.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(A3SS.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
A3SS.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(A3SS.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A3SS.WT_vs_RBM20KO = A3SS.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(A3SS.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

A3SS.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNpromKO/A3SS.MATS.JC.txt")
A3SS.WT_vs_TTNpromKO = A3SS.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A3SS.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(A3SS.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
A3SS.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(A3SS.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
A3SS.WT_vs_TTNpromKO = A3SS.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A3SS.WT_vs_TTNpromKO = A3SS.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, IncLevelDifference, FDR))
colnames(A3SS.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

A3SS.WT_vs_TTNKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNKO/A3SS.MATS.JC.txt")
A3SS.WT_vs_TTNKO = A3SS.WT_vs_TTNKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A3SS.WT_vs_TTNKO$JC_count_WT = as.character(lapply(A3SS.WT_vs_TTNKO$IJC_SAMPLE_1, get_JCcounts))
A3SS.WT_vs_TTNKO$JC_count_TTNKO = as.character(lapply(A3SS.WT_vs_TTNKO$IJC_SAMPLE_2, get_JCcounts))
A3SS.WT_vs_TTNKO = A3SS.WT_vs_TTNKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A3SS.WT_vs_TTNKO = A3SS.WT_vs_TTNKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNKO, IncLevelDifference, FDR))
colnames(A3SS.WT_vs_TTNKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNKO", "delta_PSI_TTNKOvsWT", "FDR_TTNKOvsWT")

A3SS_intersect = A3SS.WT_vs_RBM20KO %>% dplyr::inner_join(A3SS.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT")) %>%
  dplyr::inner_join(A3SS.WT_vs_TTNKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT"))
A3SS_intersect = A3SS_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_TTNKO)) %>% dplyr::filter(JC_sum > 16)

A5SS.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/A5SS.MATS.JC.txt"))
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A5SS.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(A5SS.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
A5SS.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(A5SS.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A5SS.WT_vs_RBM20KO = A5SS.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(A5SS.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

A5SS.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNpromKO/A5SS.MATS.JC.txt")
A5SS.WT_vs_TTNpromKO = A5SS.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A5SS.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(A5SS.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
A5SS.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(A5SS.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
A5SS.WT_vs_TTNpromKO = A5SS.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A5SS.WT_vs_TTNpromKO = A5SS.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, IncLevelDifference, FDR))
colnames(A5SS.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

A5SS.WT_vs_TTNKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNKO/A5SS.MATS.JC.txt")
A5SS.WT_vs_TTNKO = A5SS.WT_vs_TTNKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_lE:", longExonStart_0base, "-", longExonEnd, "_sE:", shortES, "-", shortEE, "_fE:", flankingES, "-", flankingEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
A5SS.WT_vs_TTNKO$JC_count_WT = as.character(lapply(A5SS.WT_vs_TTNKO$IJC_SAMPLE_1, get_JCcounts))
A5SS.WT_vs_TTNKO$JC_count_TTNKO = as.character(lapply(A5SS.WT_vs_TTNKO$IJC_SAMPLE_2, get_JCcounts))
A5SS.WT_vs_TTNKO = A5SS.WT_vs_TTNKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
A5SS.WT_vs_TTNKO = A5SS.WT_vs_TTNKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNKO, IncLevelDifference, FDR))
colnames(A5SS.WT_vs_TTNKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNKO", "delta_PSI_TTNKOvsWT", "FDR_TTNKOvsWT")

A5SS_intersect = A5SS.WT_vs_RBM20KO %>% dplyr::inner_join(A5SS.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT")) %>%
  dplyr::inner_join(A5SS.WT_vs_TTNKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT"))
A5SS_intersect = A5SS_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_TTNKO)) %>% dplyr::filter(JC_sum > 16)

RI.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/RI.MATS.JC.txt"))
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_riE:", riExonStart_0base, "-", riExonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
RI.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(RI.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
RI.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(RI.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
RI.WT_vs_RBM20KO = RI.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(RI.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

RI.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNpromKO/RI.MATS.JC.txt")
RI.WT_vs_TTNpromKO = RI.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_riE:", riExonStart_0base, "-", riExonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
RI.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(RI.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
RI.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(RI.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
RI.WT_vs_TTNpromKO = RI.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
RI.WT_vs_TTNpromKO = RI.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, IncLevelDifference, FDR))
colnames(RI.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

RI.WT_vs_TTNKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNKO/RI.MATS.JC.txt")
RI.WT_vs_TTNKO = RI.WT_vs_TTNKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_riE:", riExonStart_0base, "-", riExonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
RI.WT_vs_TTNKO$JC_count_WT = as.character(lapply(RI.WT_vs_TTNKO$IJC_SAMPLE_1, get_JCcounts))
RI.WT_vs_TTNKO$JC_count_TTNKO = as.character(lapply(RI.WT_vs_TTNKO$IJC_SAMPLE_2, get_JCcounts))
RI.WT_vs_TTNKO = RI.WT_vs_TTNKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
RI.WT_vs_TTNKO = RI.WT_vs_TTNKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNKO, IncLevelDifference, FDR))
colnames(RI.WT_vs_TTNKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNKO", "delta_PSI_TTNKOvsWT", "FDR_TTNKOvsWT")

RI_intersect = RI.WT_vs_RBM20KO %>% dplyr::inner_join(RI.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT")) %>%
  dplyr::inner_join(RI.WT_vs_TTNKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT"))
RI_intersect = RI_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_TTNKO)) %>% dplyr::filter(JC_sum > 16)

MXE.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/MXE.MATS.JC.txt"))
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E1:", X1stExonStart_0base, "-", X1stExonEnd, "_E2:", X2ndExonStart_0base, "-", X2ndExonEnd, "_uE:", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
MXE.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(MXE.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
MXE.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(MXE.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
MXE.WT_vs_RBM20KO = MXE.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(MXE.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

MXE.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNpromKO/MXE.MATS.JC.txt")
MXE.WT_vs_TTNpromKO = MXE.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E1:", X1stExonStart_0base, "-", X1stExonEnd, "_E2:", X2ndExonStart_0base, "-", X2ndExonEnd, "_uE:", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
MXE.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(MXE.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
MXE.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(MXE.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
MXE.WT_vs_TTNpromKO = MXE.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
MXE.WT_vs_TTNpromKO = MXE.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, IncLevelDifference, FDR))
colnames(MXE.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

MXE.WT_vs_TTNKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNKO/MXE.MATS.JC.txt")
MXE.WT_vs_TTNKO = MXE.WT_vs_TTNKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E1:", X1stExonStart_0base, "-", X1stExonEnd, "_E2:", X2ndExonStart_0base, "-", X2ndExonEnd, "_uE:", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
MXE.WT_vs_TTNKO$JC_count_WT = as.character(lapply(MXE.WT_vs_TTNKO$IJC_SAMPLE_1, get_JCcounts))
MXE.WT_vs_TTNKO$JC_count_TTNKO = as.character(lapply(MXE.WT_vs_TTNKO$IJC_SAMPLE_2, get_JCcounts))
MXE.WT_vs_TTNKO = MXE.WT_vs_TTNKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
MXE.WT_vs_TTNKO = MXE.WT_vs_TTNKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNKO, IncLevelDifference, FDR))
colnames(MXE.WT_vs_TTNKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNKO", "delta_PSI_TTNKOvsWT", "FDR_TTNKOvsWT")

MXE_intersect = MXE.WT_vs_RBM20KO %>% dplyr::inner_join(MXE.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT")) %>%
  dplyr::inner_join(MXE.WT_vs_TTNKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT"))
MXE_intersect = MXE_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_TTNKO)) %>% dplyr::filter(JC_sum > 16)

SE.WT_vs_RBM20KO <- read.delim(paste0("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_RBM20KO/SE.MATS.JC.txt"))
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E:", exonStart_0base, "-", exonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
SE.WT_vs_RBM20KO$JC_count_WT = as.character(lapply(SE.WT_vs_RBM20KO$IJC_SAMPLE_1, get_JCcounts))
SE.WT_vs_RBM20KO$JC_count_RBM20KO = as.character(lapply(SE.WT_vs_RBM20KO$IJC_SAMPLE_2, get_JCcounts))
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
SE.WT_vs_RBM20KO = SE.WT_vs_RBM20KO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_RBM20KO, IncLevelDifference, FDR))
colnames(SE.WT_vs_RBM20KO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_RBM20KO", "delta_PSI_RBM20KOvsWT", "FDR_RBM20KOvsWT")

SE.WT_vs_TTNpromKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNpromKO/SE.MATS.JC.txt")
SE.WT_vs_TTNpromKO = SE.WT_vs_TTNpromKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E:", exonStart_0base, "-", exonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
SE.WT_vs_TTNpromKO$JC_count_WT = as.character(lapply(SE.WT_vs_TTNpromKO$IJC_SAMPLE_1, get_JCcounts))
SE.WT_vs_TTNpromKO$JC_count_TTNpromKO = as.character(lapply(SE.WT_vs_TTNpromKO$IJC_SAMPLE_2, get_JCcounts))
SE.WT_vs_TTNpromKO = SE.WT_vs_TTNpromKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
SE.WT_vs_TTNpromKO = SE.WT_vs_TTNpromKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNpromKO, IncLevelDifference, FDR))
colnames(SE.WT_vs_TTNpromKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNpromKO", "delta_PSI_TTNpromKOvsWT", "FDR_TTNpromKOvsWT")

SE.WT_vs_TTNKO <- read.delim("/home/shared_folder/splicing_data/RUES2_rMATS_data/WT_vs_TTNKO/SE.MATS.JC.txt")
SE.WT_vs_TTNKO = SE.WT_vs_TTNKO %>% dplyr::mutate(ID = paste0(chr, "_", strand, "_E:", exonStart_0base, "-", exonEnd, "_uE:", upstreamES, "-", upstreamEE, "_dE:", downstreamES, "-", downstreamEE)) %>%
  dplyr::select(c("GeneID", "geneSymbol", "ID", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "FDR", "IncLevelDifference"))
SE.WT_vs_TTNKO$JC_count_WT = as.character(lapply(SE.WT_vs_TTNKO$IJC_SAMPLE_1, get_JCcounts))
SE.WT_vs_TTNKO$JC_count_TTNKO = as.character(lapply(SE.WT_vs_TTNKO$IJC_SAMPLE_2, get_JCcounts))
SE.WT_vs_TTNKO = SE.WT_vs_TTNKO %>% dplyr::mutate(IncLevelDifference = - IncLevelDifference)
SE.WT_vs_TTNKO = SE.WT_vs_TTNKO %>% dplyr::select(c(ID, geneSymbol, GeneID, JC_count_WT, JC_count_TTNKO, IncLevelDifference, FDR))
colnames(SE.WT_vs_TTNKO) = c("variant_ID", "geneSymbol", "geneID", "JC_WT", "JC_TTNKO", "delta_PSI_TTNKOvsWT", "FDR_TTNKOvsWT")

SE_intersect = SE.WT_vs_RBM20KO %>% dplyr::inner_join(SE.WT_vs_TTNpromKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT")) %>%
  dplyr::inner_join(SE.WT_vs_TTNKO, by = c("variant_ID", "geneSymbol", "geneID", "JC_WT"))
SE_intersect = SE_intersect %>% dplyr::mutate(JC_sum = as.numeric(JC_WT) + as.numeric(JC_RBM20KO) + as.numeric(JC_TTNpromKO) + as.numeric(JC_TTNKO)) %>% dplyr::filter(JC_sum > 16)

merge = rbind(A3SS_intersect %>% dplyr::mutate(event = "A3SS"),
              A5SS_intersect %>% dplyr::mutate(event = "A5SS"),
              RI_intersect %>% dplyr::mutate(event = "RI"),
              MXE_intersect %>% dplyr::mutate(event = "MXE"),
              SE_intersect %>% dplyr::mutate(event = "SE"))
merge[is.na(merge)] = ""

transC_lists <- read_excel("/home/shared_folder/chromatin_data/eCLIP_data/transC_lists.xlsx") %>% dplyr::filter(HC_eCLIP == "Yes")
merge = merge %>% dplyr::mutate(eCLIP = ifelse(geneSymbol == "TTN", "TTN", ifelse(geneSymbol %in% transC_lists$Gene, "RBM20_target", "non_target")))

merge = merge %>%
  dplyr::mutate(RBM20KO = ifelse(FDR_RBM20KOvsWT < 0.05 & delta_PSI_RBM20KOvsWT > 0.10, "upregulated", ifelse(FDR_RBM20KOvsWT < 0.05 & delta_PSI_RBM20KOvsWT < -0.10, "downregulated", "NS"))) %>%
  dplyr::mutate(TTNpromKO = ifelse(FDR_TTNpromKOvsWT < 0.05 & delta_PSI_TTNpromKOvsWT > 0.10, "upregulated", ifelse(FDR_TTNpromKOvsWT < 0.05 & delta_PSI_TTNpromKOvsWT < -0.10, "downregulated", "NS"))) %>%
  dplyr::mutate(TTNKO = ifelse(FDR_TTNKOvsWT < 0.05 & delta_PSI_TTNKOvsWT > 0.10, "upregulated", ifelse(FDR_TTNKOvsWT < 0.05 & delta_PSI_TTNKOvsWT < -0.10, "downregulated", "NS")))
merge_sign = merge %>% dplyr::filter(RBM20KO != "NS" | TTNpromKO != "NS" | TTNKO != "NS")
merge_sign = merge_sign %>% dplyr::mutate(gene_name = ifelse(eCLIP == "RBM20_target", geneSymbol, ""))

#splicing events of eCLIP targets correlation - supplementary figure 8 panel C
lin_targets = lm(as.numeric((merge_sign %>% dplyr::filter(eCLIP == "RBM20_target"))$delta_PSI_RBM20KOvsWT) ~ as.numeric((merge_sign %>% dplyr::filter(eCLIP == "RBM20_target"))$delta_PSI_TTNpromKOvsWT))
ggplot(merge_sign %>% dplyr::arrange(desc(eCLIP)), aes(x = delta_PSI_RBM20KOvsWT, y = delta_PSI_TTNpromKOvsWT)) + geom_point(aes(fill = eCLIP), shape = 21, size = 3, alpha = 0.7, color = "black") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") + geom_vline(aes(xintercept = 0), linetype = "dotted") +
  scale_fill_manual(values = c("white", "#1C8356", "#FA0087")) + scale_y_continuous(limits = c(-1, 1)) + scale_x_continuous(limits = c(-1, 1)) +
  geom_abline(aes(slope = lin_targets[["coefficients"]][2], intercept = lin_targets[["coefficients"]][1]), linetype = "dashed", color = "#1C8356", size = 1)
ggsave("/home/shared_folder/output_plots/linear_correlation_deltaPSI_RBM20KO_TTNpromKO.pdf")

lin_targets = lm(as.numeric((merge_sign %>% dplyr::filter(eCLIP == "RBM20_target"))$delta_PSI_RBM20KOvsWT) ~ as.numeric((merge_sign %>% dplyr::filter(eCLIP == "RBM20_target"))$delta_PSI_TTNKOvsWT))
ggplot(merge_sign %>% dplyr::arrange(desc(eCLIP)), aes(x = delta_PSI_RBM20KOvsWT, y = delta_PSI_TTNKOvsWT)) + geom_point(aes(fill = eCLIP), shape = 21, size = 3, alpha = 0.7, color = "black") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") + geom_vline(aes(xintercept = 0), linetype = "dotted") +
  scale_fill_manual(values = c("white", "#1C8356", "#FA0087")) + scale_y_continuous(limits = c(-1, 1)) + scale_x_continuous(limits = c(-1, 1)) +
  geom_abline(aes(slope = lin_targets[["coefficients"]][2], intercept = lin_targets[["coefficients"]][1]), linetype = "dashed", color = "#1C8356", size = 1)
ggsave("/home/shared_folder/output_plots/linear_correlation_deltaPSI_RBM20KO_TTNKO.pdf")


#up and downregulation correlation - supplementary figure 8 panel D
merge_eCLIP = merge_sign %>% dplyr::mutate(sign_RBM20KO = ifelse(delta_PSI_RBM20KOvsWT > 0, "RBM20KO+", "RBM20KO-"),
                                          sign_TTNpromKO = ifelse(delta_PSI_TTNpromKOvsWT > 0, "TTNpromKO+", "TTNpromKO-"),
                                          sign_TTNKO = ifelse(delta_PSI_TTNKOvsWT > 0, "TTNKO+", "TTNKO-"), n = 1) %>%
  dplyr::filter(eCLIP == "RBM20_target")

rand_TTNpromKO = as.data.frame(crossing(sign_RBM20KO = c("RBM20KO+", "RBM20KO-"), sign_TTNpromKO = c("TTNpromKO+", "TTNpromKO-")))
for(i in 1:10000){
  merge_eCLIP$random = sample(merge_eCLIP$sign_TTNpromKO)
  summ = merge_eCLIP %>% dplyr::group_by(sign_RBM20KO, random) %>% dplyr::summarise(x = sum(n))
  colnames(summ) = c("sign_RBM20KO", "sign_TTNpromKO", paste0("random", i))
  rand_TTNpromKO = rand_TTNpromKO %>% dplyr::left_join(summ, by = c("sign_RBM20KO", "sign_TTNpromKO"))
}
melt = melt(rand_TTNpromKO)
summ = merge_eCLIP %>% dplyr::group_by(sign_RBM20KO, sign_TTNpromKO) %>% dplyr::summarise(ref = sum(n))
melt = melt %>% dplyr::left_join(summ, by = c("sign_RBM20KO", "sign_TTNpromKO"))
melt = melt %>% dplyr::mutate(higher = ifelse(value >= ref, 0.0001, 0), lower = ifelse(value <= ref, 0.0001, 0))
summ_TTNpromKO = melt %>% dplyr::group_by(sign_RBM20KO, sign_TTNpromKO) %>% dplyr::summarise(pval_higher = sum(higher), pval_lower = sum(lower), observed = mean(ref), random = mean(value)) %>%
  dplyr::mutate(enrichment = observed/random) %>% dplyr::mutate(sign = ifelse(pval_higher < 0.001 | pval_lower < 0.001, "***", ifelse(pval_higher < 0.01 | pval_lower < 0.01, "**", ifelse(pval_higher < 0.05 | pval_lower < 0.05, "*", ""))))
write.csv(summ_TTNpromKO, "/home/shared_folder/output_plots/pval_random_model_eCLIP_targets_TTNpromKO.csv")
ggplot(melt(summ_TTNpromKO %>% dplyr::mutate(category = paste0(sign_RBM20KO, "\n", sign_TTNpromKO)) %>% dplyr::select(c(category, random, observed, sign)))) + ylab("") + xlab("") +
  geom_col(aes(x = category, y = value, fill = variable), position = "dodge") + geom_text(aes(x = category, y = 36.5, label = sign)) + scale_fill_manual(values = c("grey", "#009E73")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("/home/shared_folder/output_plots/enrichment_eCLIP_target_TTNpromKO.pdf")

rand_TTNKO = as.data.frame(crossing(sign_RBM20KO = c("RBM20KO+", "RBM20KO-"), sign_TTNKO = c("TTNKO+", "TTNKO-")))
for(i in 1:10000){
  merge_eCLIP$random = sample(merge_eCLIP$sign_TTNKO)
  summ = merge_eCLIP %>% dplyr::group_by(sign_RBM20KO, random) %>% dplyr::summarise(x = sum(n))
  colnames(summ) = c("sign_RBM20KO", "sign_TTNKO", paste0("random", i))
  rand_TTNKO = rand_TTNKO %>% dplyr::left_join(summ, by = c("sign_RBM20KO", "sign_TTNKO"))
}
melt = melt(rand_TTNKO)
summ = merge_eCLIP %>% dplyr::group_by(sign_RBM20KO, sign_TTNKO) %>% dplyr::summarise(ref = sum(n))
melt = melt %>% dplyr::left_join(summ, by = c("sign_RBM20KO", "sign_TTNKO"))
melt = melt %>% dplyr::mutate(higher = ifelse(value >= ref, 0.0001, 0), lower = ifelse(value <= ref, 0.0001, 0))
summ_TTNKO = melt %>% dplyr::group_by(sign_RBM20KO, sign_TTNKO) %>% dplyr::summarise(pval_higher = sum(higher), pval_lower = sum(lower), observed = mean(ref), random = mean(value)) %>%
  dplyr::mutate(enrichment = observed/random) %>% dplyr::mutate(sign = ifelse(pval_higher < 0.001 | pval_lower < 0.001, "***", ifelse(pval_higher < 0.01 | pval_lower < 0.01, "**", ifelse(pval_higher < 0.05 | pval_lower < 0.05, "*", ""))))
write.csv(summ_TTNKO, "/home/shared_folder/output_plots/pval_random_model_eCLIP_targets_TTNKO.csv")
ggplot(melt(summ_TTNKO %>% dplyr::mutate(category = paste0(sign_RBM20KO, "\n", sign_TTNKO)) %>% dplyr::select(c(category, random, observed, sign)))) + ylab("") + xlab("") +
  geom_col(aes(x = category, y = value, fill = variable), position = "dodge") + geom_text(aes(x = category, y = 36.5, label = sign)) + scale_fill_manual(values = c("grey", "#009E73")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("/home/shared_folder/output_plots/enrichment_eCLIP_target_TTNKO.pdf")

#wilcoxon test for effect size and sign correlation - supplementary figure 8 panel E
melt_eCLIP = melt(merge_eCLIP %>% dplyr::select(c(variant_ID, sign_RBM20KO, delta_PSI_TTNpromKOvsWT, delta_PSI_TTNKOvsWT)))
melt_eCLIP = melt_eCLIP %>% dplyr::mutate(sign = ifelse(sign_RBM20KO == "RBM20KO+", 1, -1)) %>% dplyr::mutate(product = value*sign)
tab = data.frame(variable = c(), pval = c())
test = wilcox.test(x = (melt_eCLIP %>% dplyr::filter(variable == "delta_PSI_TTNpromKOvsWT"))$product, mu = 0)
tab = rbind(tab, data.frame(variable = c("delta_PSI_TTNpromKOvsWT"), pval = signif(c(test[["p.value"]]), digits = 3)))
test = wilcox.test(x = (melt_eCLIP %>% dplyr::filter(variable == "delta_PSI_TTNKOvsWT"))$product, mu = 0)
tab = rbind(tab, data.frame(variable = c("delta_PSI_TTNKOvsWT"), pval = signif(c(test[["p.value"]]), digits = 3)))
tab = tab %>% dplyr::mutate(pval = ifelse(pval < 0.001, paste0("p = ", pval, " ***"), ifelse(pval < 0.01, paste0("p = ", pval, " **"), ifelse(pval < 0.05, paste0("p = ", pval, " *"), paste0("p = ", pval)))))
melt_eCLIP = melt_eCLIP %>% dplyr::left_join(tab, by = "variable")
melt_eCLIP = melt_eCLIP %>% dplyr::mutate(variable = ifelse(variable == "delta_PSI_TTNpromKOvsWT", "TTN promoter KO vs RBM20 KO", "TTN protein KO vs RBM20 KO"))
ggplot(melt_eCLIP) + geom_violin(aes(x = variable, y = product), fill = "#009E73", alpha = 0.5) + geom_boxplot(aes(x = variable, y = product), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_text(aes(x = variable, y = 1, label = pval)) + xlab("") + ylab("delta PSI concordance\neffect size TTNKO") +
  scale_y_continuous(limits = c(-1, 1))
ggsave("/home/shared_folder/output_plots/wilcoxon_eCLIP_targets.pdf")

melt_sign = melt(merge_sign %>% dplyr::select(c(variant_ID, sign_RBM20KO, delta_PSI_TTNpromKOvsWT, delta_PSI_TTNKOvsWT)))
melt_sign = melt_sign %>% dplyr::mutate(sign = ifelse(sign_RBM20KO == "RBM20KO+", 1, -1)) %>% dplyr::mutate(product = value*sign)
tab = data.frame(variable = c(), pval = c())
test = wilcox.test(x = (melt_sign %>% dplyr::filter(variable == "delta_PSI_TTNpromKOvsWT"))$product, mu = 0)
tab = rbind(tab, data.frame(variable = c("delta_PSI_TTNpromKOvsWT"), pval = signif(c(test[["p.value"]]), digits = 3)))
test = wilcox.test(x = (melt_sign %>% dplyr::filter(variable == "delta_PSI_TTNKOvsWT"))$product, mu = 0)
tab = rbind(tab, data.frame(variable = c("delta_PSI_TTNKOvsWT"), pval = signif(c(test[["p.value"]]), digits = 3)))
tab = tab %>% dplyr::mutate(pval = ifelse(pval < 0.001, paste0("p = ", pval, " ***"), ifelse(pval < 0.01, paste0("p = ", pval, " **"), ifelse(pval < 0.05, paste0("p = ", pval, " *"), paste0("p = ", pval)))))
melt_sign = melt_sign %>% dplyr::left_join(tab, by = "variable")
melt_sign = melt_sign %>% dplyr::mutate(variable = ifelse(variable == "delta_PSI_TTNpromKOvsWT", "TTN promoter KO vs RBM20 KO", "TTN protein KO vs RBM20 KO"))
ggplot(melt_sign) + geom_violin(aes(x = variable, y = product), fill = "#D55E00", alpha = 0.5) + geom_boxplot(aes(x = variable, y = product), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_text(aes(x = variable, y = 1, label = pval)) + xlab("") + ylab("delta PSI concordance\neffect size TTNKO") +
  scale_y_continuous(limits = c(-1, 1))
ggsave("/home/shared_folder/output_plots/wilcoxon_all.pdf")


melt_eCLIP = melt(merge_eCLIP %>% dplyr::select(c(variant_ID, sign_TTNpromKO, sign_TTNKO, delta_PSI_RBM20KOvsWT)))
melt_eCLIP = melt_eCLIP %>% dplyr::mutate(sign_TTNpromKO = ifelse(sign_TTNpromKO == "TTNpromKO+", 1, -1)) %>% dplyr::mutate(product_prom = value*sign_TTNpromKO)
melt_eCLIP = melt_eCLIP %>% dplyr::mutate(sign_TTNKO = ifelse(sign_TTNKO == "TTNKO+", 1, -1)) %>% dplyr::mutate(product_prot = value*sign_TTNKO)
tab = data.frame(variable = c(), pval = c())
test = wilcox.test(x = (melt_eCLIP)$product_prom, mu = 0)
tab = rbind(tab, data.frame(variable = c("product_prom"), pval = signif(c(test[["p.value"]]), digits = 3)))
test = wilcox.test(x = (melt_eCLIP)$product_prot, mu = 0)
tab = rbind(tab, data.frame(variable = c("product_prot"), pval = signif(c(test[["p.value"]]), digits = 3)))
tab = tab %>% dplyr::mutate(pval = ifelse(pval < 0.001, paste0("p = ", pval, " ***"), ifelse(pval < 0.01, paste0("p = ", pval, " **"), ifelse(pval < 0.05, paste0("p = ", pval, " *"), paste0("p = ", pval)))))
melt_eCLIP = melt(melt_eCLIP %>% dplyr::select(c(variant_ID, product_prom, product_prot)))
melt_eCLIP = melt_eCLIP %>% dplyr::left_join(tab, by = "variable")
melt_eCLIP = melt_eCLIP %>% dplyr::mutate(variable = ifelse(variable == "product_prom", "TTN promoter KO vs RBM20 KO", "TTN protein KO vs RBM20 KO"))
ggplot(melt_eCLIP) + geom_violin(aes(x = variable, y = value), fill = "#009E73", alpha = 0.5) + geom_boxplot(aes(x = variable, y = value), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_text(aes(x = variable, y = 1, label = pval)) + xlab("") + ylab("delta PSI concordance\neffect size RBM20KO") +
  scale_y_continuous(limits = c(-1, 1))
ggsave("/home/shared_folder/output_plots/wilcoxon_eCLIP_targets_RBM20_effect.pdf")

melt_sign = melt(merge_sign %>% dplyr::select(c(variant_ID, sign_TTNpromKO, sign_TTNKO, delta_PSI_RBM20KOvsWT)))
melt_sign = melt_sign %>% dplyr::mutate(sign_TTNpromKO = ifelse(sign_TTNpromKO == "TTNpromKO+", 1, -1)) %>% dplyr::mutate(product_prom = value*sign_TTNpromKO)
melt_sign = melt_sign %>% dplyr::mutate(sign_TTNKO = ifelse(sign_TTNKO == "TTNKO+", 1, -1)) %>% dplyr::mutate(product_prot = value*sign_TTNKO)
tab = data.frame(variable = c(), pval = c())
test = wilcox.test(x = (melt_sign)$product_prom, mu = 0)
tab = rbind(tab, data.frame(variable = c("product_prom"), pval = signif(c(test[["p.value"]]), digits = 3)))
test = wilcox.test(x = (melt_sign)$product_prot, mu = 0)
tab = rbind(tab, data.frame(variable = c("product_prot"), pval = signif(c(test[["p.value"]]), digits = 3)))
tab = tab %>% dplyr::mutate(pval = ifelse(pval < 0.001, paste0("p = ", pval, " ***"), ifelse(pval < 0.01, paste0("p = ", pval, " **"), ifelse(pval < 0.05, paste0("p = ", pval, " *"), paste0("p = ", pval)))))
melt_sign = melt(melt_sign %>% dplyr::select(c(variant_ID, product_prom, product_prot)))
melt_sign = melt_sign %>% dplyr::left_join(tab, by = "variable")
melt_sign = melt_sign %>% dplyr::mutate(variable = ifelse(variable == "product_prom", "TTN promoter KO vs RBM20 KO", "TTN protein KO vs RBM20 KO"))
ggplot(melt_sign) + geom_violin(aes(x = variable, y = value), fill = "#D55E00", alpha = 0.5) + geom_boxplot(aes(x = variable, y = value), alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_text(aes(x = variable, y = 1, label = pval)) + xlab("") + ylab("delta PSI concordance\neffect size RBM20KO") +
  scale_y_continuous(limits = c(-1, 1))
ggsave("/home/shared_folder/output_plots/wilcoxon_all_RBM20_effect.pdf")
