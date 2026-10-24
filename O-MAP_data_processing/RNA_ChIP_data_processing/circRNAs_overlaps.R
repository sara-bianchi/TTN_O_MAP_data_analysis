library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rbowtie)
library(crisprBwa)
library(dplyr)
library(readxl)
library(ggplot2)
library(viridis)
library(ggsignif)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 45, hjust = 1)))

dir.create("/home/shared_folder/output")

#alignemnt of probes
probes <- read_excel("/home/shared_folder/probes.xlsx", header = TRUE)
landing_pad = "TTTACATCATCATACATCATCATACATCATCAT"

probes = data.frame(probe = c(probes$`pre-TTN - Odd Probes~SABER1 (O-MAP probes)`, probes$`pre-TTN - Even Probes~SABER1 (O-MAP probes)`))
probes$probe = probes$probe %>% stringr::str_remove(pattern = landing_pad)

genome_dir = "/home/shared_folder/bwa_index/"
dir.create(genome_dir)
bowtie_build(fasta,
              outdir = genome_dir,
              force = TRUE,
              prefix = "hg38")

bwa_index = "/home/shared_folder/bwa_index/hg38"

res = data.frame(query = c(), chr = c(), pos = c(), strand = c(), n_mismatches = c())
for(i in probes$probe){
  res_i = runBwa(i, bwa_index = bwa_index, n_mismatches = 3)
  res = rbind(res, res_i)
}
len = c()
for(i in 1:length(rownames(res))){
  len = append(len, nchar(res$query[i]))
}
res$len = len
res = res %>% dplyr::mutate(end = pos + len)

#importation of GTF file
Homo_sapiens.GRCh38.113 <- read.delim("/home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf", header = FALSE, comment.char="#")
get_geneID <- function(str){
  val = strsplit(str, "; ")[[1]][1]
  val = strsplit(val, " ")[[1]][2]
  return(val)
}
get_transcriptID <- function(str){
  val = strsplit(str, "; ")[[1]][3]
  val = strsplit(val, " ")[[1]][2]
  return(val)
}
get_gene_name <- function(str){
  val = strsplit(str, "; ")[[1]][6]
  val = strsplit(val, " ")[[1]][2]
  return(val)
}
get_transcript_name <- function(str){
  val = strsplit(str, "; ")[[1]][9]
  val = strsplit(val, " ")[[1]][2]
  return(val)
}

gtf = as.data.frame(Homo_sapiens.GRCh38.113)
gtf$geneID = as.character(lapply(gtf$V9, get_geneID))
gtf$transcriptID = as.character(lapply(gtf$V9, get_transcriptID))
gtf$gene_name = as.character(lapply(gtf$V9, get_gene_name))
gtf$transcript_name = as.character(lapply(gtf$V9, get_transcript_name))

gtf_TTN = gtf %>% dplyr::filter(gene_name == "TTN")
gtf_exons = gtf_TTN %>% dplyr::filter(V3 == "exon") %>% dplyr::select(c(V4, V5)) %>% dplyr::distinct()

#bed file introns
gtf_exons = gtf_exons %>% dplyr::arrange(V4)
gtf_introns = data.frame(chr = c(), start = c(), end = c())
for(i in 2:length(rownames(gtf_exons))){
  start = gtf_exons$V5[i - 1]
  end = gtf_exons$V4[i]
  gtf_introns = rbind(gtf_introns, data.frame(chr = c("chr2"), start = c(start), end = c(end)))
}

#import reads aligned on TTN
WT1 <- read.delim("/home/shared_folder/ratio_pulldown_input/WT1.bedGraph", header=FALSE)
WT2 <- read.delim("/home/shared_folder/ratio_pulldown_input/WT2.bedGraph", header=FALSE)
WT3 <- read.delim("/home/shared_folder/ratio_pulldown_input/WT3.bedGraph", header=FALSE)
WT4 <- read.delim("/home/shared_folder/ratio_pulldown_input/WT4.bedGraph", header=FALSE)

KO1 <- read.delim("/home/shared_folder/ratio_pulldown_input/KO1.bedGraph", header=FALSE)
KO2 <- read.delim("/home/shared_folder/ratio_pulldown_input/KO2.bedGraph", header=FALSE)
KO3 <- read.delim("/home/shared_folder/ratio_pulldown_input/KO3.bedGraph", header=FALSE)
KO4 <- read.delim("/home/shared_folder/ratio_pulldown_input/KO4.bedGraph", header=FALSE)

gtf_introns = gtf_introns %>% dplyr::arrange(start)
WT1_introns = c()
WT2_introns = c()
WT3_introns = c()
WT4_introns = c()
KO1_introns = c()
KO2_introns = c()
KO3_introns = c()
KO4_introns = c()
for(i in 1:nrow(gtf_introns)){
  bed = WT1 %>% dplyr::filter(V2 >= gtf_introns$start[i] & V2 < gtf_introns$end[i])
  WT1_introns = append(WT1_introns, sum(bed$V4)/nrow(bed))
  bed = WT2 %>% dplyr::filter(V2 >= gtf_introns$start[i] & V2 < gtf_introns$end[i])
  WT2_introns = append(WT2_introns, sum(bed$V4)/nrow(bed))
  bed = WT3 %>% dplyr::filter(V2 >= gtf_introns$start[i] & V2 < gtf_introns$end[i])
  WT3_introns = append(WT3_introns, sum(bed$V4)/nrow(bed))
  bed = WT4 %>% dplyr::filter(V2 >= gtf_introns$start[i] & V2 < gtf_introns$end[i])
  WT4_introns = append(WT4_introns, sum(bed$V4)/nrow(bed))
  bed = KO1 %>% dplyr::filter(V2 >= gtf_introns$start[i] & V2 < gtf_introns$end[i])
  KO1_introns = append(KO1_introns, sum(bed$V4)/nrow(bed))
  bed = KO2 %>% dplyr::filter(V2 >= gtf_introns$start[i] & V2 < gtf_introns$end[i])
  KO2_introns = append(KO2_introns, sum(bed$V4)/nrow(bed))
  bed = KO3 %>% dplyr::filter(V2 >= gtf_introns$start[i] & V2 < gtf_introns$end[i])
  KO3_introns = append(KO3_introns, sum(bed$V4)/nrow(bed))
  bed = KO4 %>% dplyr::filter(V2 >= gtf_introns$start[i] & V2 < gtf_introns$end[i])
  KO4_introns = append(KO4_introns, sum(bed$V4)/nrow(bed))
}

gtf_introns$WT1 = WT1_introns
gtf_introns$WT2 = WT2_introns
gtf_introns$WT3 = WT3_introns
gtf_introns$WT4 = WT4_introns
gtf_introns$KO1 = KO1_introns
gtf_introns$KO2 = KO2_introns
gtf_introns$KO3 = KO3_introns
gtf_introns$KO4 = KO4_introns

gtf_introns = gtf_introns %>% dplyr::arrange(start)
probes_n = c()
for(i in 1:nrow(gtf_introns)){
  res_i = res %>% dplyr::filter(pos >= gtf_introns$start[i] & end <= gtf_introns$end[i])
  probes_n = append(probes_n, nrow(res_i))
}
gtf_introns$probes_n = probes_n

gtf_introns$start = as.character(gtf_introns$start)
gtf_introns$end = as.character(gtf_introns$end)
gtf_introns = na.omit(gtf_introns)
gtf_introns$probes_n = as.character(gtf_introns$probes_n)
melt_introns = reshape2::melt(gtf_introns)
melt_introns = melt_introns %>% dplyr::mutate(condition = ifelse(variable %in% c("WT1", "WT2", "WT3", "WT4"), "WT", "KO"))

p.val_introns = data.frame(start = c(), p.val = c())
for(i in 1:nrow(gtf_introns)){
  melt_i = melt_introns %>% dplyr::filter(start == gtf_introns$start[i])
  w = wilcox.test(x = (melt_i %>% dplyr::filter(condition == "WT"))$value, y = (melt_i %>% dplyr::filter(condition == "KO"))$value, alternative = "two.sided")
  p.val_introns = rbind(p.val_introns, data.frame(start = c(gtf_introns$start[i]), p.val = c(w[["p.value"]])))
}

gtf_introns = gtf_introns %>% dplyr::arrange(desc(start)) %>% dplyr::mutate(intron_number = 1:nrow(gtf_introns))
gtf_introns = gtf_introns %>% dplyr::left_join(p.val_introns, by = "start")
introns_info = gtf_introns %>% dplyr::select(start, p.val, intron_number)
melt_introns = melt_introns %>% dplyr::left_join(introns_info, by = "start")

introns_targeted = melt_introns %>% dplyr::filter(probes_n != 0)
introns_targeted = introns_targeted %>% dplyr::mutate(introns = paste0("intron number = ", intron_number, "\nprobes = ", probes_n, "\np.val = ", p.val))
introns_targeted = introns_targeted %>% dplyr::arrange(intron_number)
ggplot(introns_targeted, aes(x = condition, y = value)) + geom_boxplot(aes(fill = condition)) + facet_wrap(~introns)
ggsave("/home/shared_folder/output/boxplot_probes_introns.pdf")

write.csv(introns_targeted, "/home/shared_folder/output/introns_targeted.csv")
gtf_introns_targeted = gtf_introns %>% dplyr::filter(start %in% (introns_targeted %>% dplyr::filter(p.val < 0.05))$start)

#circular overlap
WT_RBM20KO_circ_anno_input <- read.csv("/home/shared_folder/circRNAs_input/WT_RBM20KO_circ_anno_input.csv")
WT_RBM20KO_circ_anno_pulldown <- read.csv("/home/shared_folder/circRNAs_pulldown/WT_RBM20KO_circ_anno_pulldown.csv")

WT_RBM20KO_circ_exp_input <- read.csv("/home/shared_folder/circRNAs_input/WT_RBM20KO_circ_exp_input.csv")
WT_RBM20KO_circ_exp_pulldown <- read.csv("/home/shared_folder/circRNAs_pulldown/WT_RBM20KO_circ_exp_pulldown.csv")

WT_RBM20KO_circ_dge_input <- read.csv("/home/shared_folder/circRNAs_input/WT_RBM20KO_circ_dge_input.tsv") %>% dplyr::select(X, logFC)
colnames(WT_RBM20KO_circ_dge_input) = c("circ_id", "logFC")
WT_RBM20KO_circ_dge_pulldown <- read.csv("/home/shared_folder/circRNAs_pulldown/WT_vs_RBM20KO_circ_dge_pulldown.csv") %>% dplyr::select(circ_id, logFC)

WT_RBM20KO_input = WT_RBM20KO_circ_exp_input %>% dplyr::left_join(WT_RBM20KO_circ_anno_input, by = "circ_id") %>% dplyr::left_join(WT_RBM20KO_circ_dge_input, by = "circ_id")
WT_RBM20KO_pulldown = WT_RBM20KO_circ_exp_pulldown %>% dplyr::left_join(WT_RBM20KO_circ_anno_pulldown, by = "circ_id") %>% dplyr::left_join(WT_RBM20KO_circ_dge_pulldown, by = "circ_id")

TTN_circ_input = WT_RBM20KO_input %>% dplyr::filter(gene_name == "TTN") %>% dplyr::mutate(enriched = ifelse(logFC < 0, "WT", "KO"))
TTN_circ_pulldown = WT_RBM20KO_pulldown %>% dplyr::filter(gene_name == "TTN") %>% dplyr::mutate(enriched = ifelse(logFC < 0, "WT", "KO"))

TTN_circ_pulldown$start = as.integer(lapply(TTN_circ_pulldown$circ_id, get_start))
TTN_circ_pulldown$end = as.integer(lapply(TTN_circ_pulldown$circ_id, get_end))

#bed files creation
summ_introns = melt_introns %>% dplyr::filter(condition == "WT") %>% dplyr::group_by(intron_number, condition) %>% dplyr::summarise(mean_WT = mean(value))
gtf_introns = gtf_introns %>% dplyr::left_join(summ_introns %>% dplyr::select(-condition), by = "intron_number")
summ_introns = melt_introns %>% dplyr::filter(condition == "KO") %>% dplyr::group_by(intron_number, condition) %>% dplyr::summarise(mean_KO = mean(value))
gtf_introns = gtf_introns %>% dplyr::left_join(summ_introns %>% dplyr::select(-condition), by = "intron_number")
bed_introns_WT = gtf_introns %>% dplyr::filter(mean_WT > mean_KO & p.val < 0.05) %>% dplyr::select(chr, start, end) %>% dplyr::mutate(chr = "2")
write.table(bed_introns_WT, "/home/shared_folder/output/bed_introns_signWT.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
bed_introns_KO = gtf_introns %>% dplyr::filter(mean_WT < mean_KO & p.val < 0.05) %>% dplyr::select(chr, start, end) %>% dplyr::mutate(chr = "2")
write.table(bed_introns_KO, "/home/shared_folder/output/bed_introns_signKO.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

bed_probes = res %>% dplyr::select(c(chr, pos, end)) %>% dplyr::mutate(chr = "2")
write.table(bed_probes, "/home/shared_folder/output/bed_probes.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

bed_introns_WT = gtf_introns %>% dplyr::filter(mean_WT > mean_KO & p.val < 0.05 & probes_n != 0) %>% dplyr::select(chr, start, end) %>% dplyr::mutate(chr = "2")
write.table(bed_introns_WT, "/home/shared_folder/output/bed_introns_signWT_overlap.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

bed_introns_KO = gtf_introns %>% dplyr::filter(mean_WT < mean_KO & p.val < 0.05 & probes_n != 0) %>% dplyr::select(chr, start, end) %>% dplyr::mutate(chr = "2")
write.table(bed_introns_KO, "/home/shared_folder/output/bed_introns_signKO_overlap.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

circ_WT = TTN_circ_pulldown %>% dplyr::filter(enriched == "WT") %>% dplyr::mutate(chr = "2") %>% dplyr::select(c(chr, start, end))
write.table(circ_WT, "/home/shared_folder/output/circ_WT_pulldown.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
circ_KO = TTN_circ_pulldown %>% dplyr::filter(enriched == "KO") %>% dplyr::mutate(chr = "2") %>% dplyr::select(c(chr, start, end))
write.table(circ_KO, "/home/shared_folder/output/circ_KO_pulldown.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
