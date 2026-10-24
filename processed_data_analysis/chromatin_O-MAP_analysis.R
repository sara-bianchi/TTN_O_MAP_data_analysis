#figure 2
library(dplyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(Rlab)
library(Hmisc)
library(ggplotify)
library(viridis)
library(hicrep)
library(circlize)
library(ComplexHeatmap)
library(PRROC)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(hjust = 1)))


##ChIP O-MAP plots
#virtual 4C normalization
WT_chr2_178806793_178808363 <- read.delim("/home/shared_folder/chromatin_data/virtual4C_data/WT_chr2_178806793_178808363_1M.bedGraph", header=FALSE)
sum = sum(WT_chr2_178806793_178808363$V4)
WT_chr2_178806793_178808363 = WT_chr2_178806793_178808363 %>% dplyr::mutate(V4 = V4/sum*nrow(WT_chr2_178806793_178808363))
WT_chr2_178806793_178808363 = WT_chr2_178806793_178808363 %>% dplyr::mutate(c4 = ".", c5 = V4, c6 = ".", c7 = V2, c8 = V3, comp = ifelse(V4 > 0, "255,0,0", "0,0,255")) %>% dplyr::select(-V4)
write.table(WT_chr2_178806793_178808363, "/home/shared_folder/chromatin_data/virtual4C_data/WT_virtual4C.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

KO_chr2_178806793_178808363 <- read.delim("/home/shared_folder/chromatin_data/virtual4C_data/RBM20KO_chr2_178806793_178808363_1M.bedGraph", header=FALSE)
sum = sum(KO_chr2_178806793_178808363$V4)
KO_chr2_178806793_178808363 = KO_chr2_178806793_178808363 %>% dplyr::mutate(V4 = V4/sum*nrow(KO_chr2_178806793_178808363))
KO_chr2_178806793_178808363 = KO_chr2_178806793_178808363 %>% dplyr::mutate(c4 = ".", c5 = V4, c6 = ".", c7 = V2, c8 = V3, comp = ifelse(V4 > 0, "255,0,0", "0,0,255")) %>% dplyr::select(-V4)
write.table(KO_chr2_178806793_178808363, "/home/shared_folder/chromatin_data/virtual4C_data/RBM20KO_virtual4C.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

#circular plots
GSE234808_WT_Interacting_Domains <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/GSE234808_WT_Interacting_Domains.bed", header=FALSE) %>% dplyr::select(-V4)
GSE234808_KO_Interacting_Domains <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/GSE234808_KO_Interacting_Domains.bed", header=FALSE) %>% dplyr::select(-V4)

cytoband = read.cytoband(species = "hg38")
len = as.data.frame(cytoband[["chr.len"]])
colnames(len) = "len"
len$chr = rownames(len)

WT_OMAP_R1 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/WT_OMAP_R1.bedGraph", header=FALSE)
colnames(WT_OMAP_R1) = c("chr", "start", "end", "WT_R1")
WT_OMAP_R2 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/WT_OMAP_R2.bedGraph", header=FALSE)
colnames(WT_OMAP_R2) = c("chr", "start", "end", "WT_R2")
WT_OMAP = WT_OMAP_R1 %>% dplyr::left_join(WT_OMAP_R2, by = c("chr", "start", "end"))
WT_OMAP = WT_OMAP %>% dplyr::mutate(WT_R1 = as.numeric(WT_R1), WT_R2 = as.numeric(WT_R2))
WT_OMAP = WT_OMAP %>% dplyr::mutate(start = as.integer(start), end = as.integer(end))
WT_OMAP = WT_OMAP %>% dplyr::mutate(value = (WT_R1 + WT_R2)/2) %>% dplyr::filter(!chr %in% c("chr2", "chrY")) %>% dplyr::select(c("chr", "start", "end", "value"))
WT_OMAP = WT_OMAP %>% dplyr::filter(value > 0) %>% dplyr::filter(!chr %in% c("chr2", "chrY"))

rm(WT_OMAP_R1, WT_OMAP_R2)
WT_OMAP = WT_OMAP %>% dplyr::left_join(len, by = "chr")
WT_OMAP = WT_OMAP %>% dplyr::filter(len > end)


KO_OMAP_R1 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/RBM20KO_OMAP_R1.bedGraph", header=FALSE)
colnames(KO_OMAP_R1) = c("chr", "start", "end", "KO_R1")
KO_OMAP_R2 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/RBM20KO_OMAP_R2.bedGraph", header=FALSE)
colnames(KO_OMAP_R2) = c("chr", "start", "end", "KO_R2")
KO_OMAP = KO_OMAP_R1 %>% dplyr::left_join(KO_OMAP_R2, by = c("chr", "start", "end"))
KO_OMAP = KO_OMAP %>% dplyr::mutate(KO_R1 = as.numeric(KO_R1), KO_R2 = as.numeric(KO_R2))
KO_OMAP = KO_OMAP %>% dplyr::mutate(start = as.integer(start), end = as.integer(end))
KO_OMAP = KO_OMAP %>% dplyr::mutate(value = (KO_R1 + KO_R2)/2) %>% dplyr::filter(!chr %in% c("chrY")) %>% dplyr::select(c("chr", "start", "end", "value"))
KO_OMAP = KO_OMAP %>% dplyr::filter(value > 0) %>% dplyr::filter(!chr %in% c("chrY"))

rm(KO_OMAP_R1, KO_OMAP_R2)
KO_OMAP = KO_OMAP %>% dplyr::left_join(len, by = "chr")
KO_OMAP = KO_OMAP %>% dplyr::filter(len > end)


WT_mock_R1 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/WT_mock_R1.bedGraph", header=FALSE)
colnames(WT_mock_R1) = c("chr", "start", "end", "WT_R1")
WT_mock_R2 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/WT_mock_R2.bedGraph", header=FALSE)
colnames(WT_mock_R2) = c("chr", "start", "end", "WT_R2")
WT_mock = WT_mock_R1 %>% dplyr::left_join(WT_mock_R2, by = c("chr", "start", "end"))
WT_mock = WT_mock %>% dplyr::mutate(WT_R1 = as.numeric(WT_R1), WT_R2 = as.numeric(WT_R2))
WT_mock = WT_mock %>% dplyr::mutate(start = as.integer(start), end = as.integer(end))
WT_mock = WT_mock %>% dplyr::mutate(value = (WT_R1 + WT_R2)/2) %>% dplyr::filter(!chr %in% c("chrY")) %>% dplyr::select(c("chr", "start", "end", "value"))
WT_mock = WT_mock %>% dplyr::filter(value > 0) %>% dplyr::filter(!chr %in% c("chrY"))

rm(WT_mock_R1, WT_mock_R2)
WT_mock = WT_mock %>% dplyr::left_join(len, by = "chr")
WT_mock = WT_mock %>% dplyr::filter(len > end)


KO_mock_R1 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/RBM20KO_mock_R1.bedGraph", header=FALSE)
colnames(KO_mock_R1) = c("chr", "start", "end", "KO_R1")
KO_mock_R2 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/RBM20KO_mock_R2.bedGraph", header=FALSE)
colnames(KO_mock_R2) = c("chr", "start", "end", "KO_R2")
KO_mock = KO_mock_R1 %>% dplyr::left_join(KO_mock_R2, by = c("chr", "start", "end"))
KO_mock = KO_mock %>% dplyr::mutate(KO_R1 = as.numeric(KO_R1), KO_R2 = as.numeric(KO_R2))
KO_mock = KO_mock %>% dplyr::mutate(start = as.integer(start), end = as.integer(end))
KO_mock = KO_mock %>% dplyr::mutate(value = (KO_R1 + KO_R2)/2) %>% dplyr::filter(!chr %in% c("chrY")) %>% dplyr::select(c("chr", "start", "end", "value"))
KO_mock = KO_mock %>% dplyr::filter(value > 0) %>% dplyr::filter(!chr %in% c("chrY"))

KO_mock = KO_mock %>% dplyr::left_join(len, by = "chr")
KO_mock = KO_mock %>% dplyr::filter(len > end)

bin_map <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/bin_map_1000000.bed", header=FALSE) %>% dplyr::filter(V1 != "chr2")
colnames(bin_map) = c("chr", "start", "end", "n")
pos = c()
for(i in 1:nrow(bin_map)){
  print(i)
  n_pos = 0
  for(j in 1:nrow(WT_mock)){
    if(bin_map$chr[i] == WT_mock$chr[j] & bin_map$start[i] < WT_mock$start[j] & bin_map$end[i] > WT_mock$end[j]){n_pos = n_pos + 1}
  }
  pos = append(pos, n_pos)
  print(n_pos)
}
bin_map$pos = pos
ggplot(bin_map) + geom_density(aes(x = log2(pos))) + geom_vline(aes(xintercept = log2(100)), color = "red", linetype = "dashed")
WT_pos = bin_map %>% dplyr::filter(pos > 100)
mock_pos = c()
for(i in 1:nrow(WT_OMAP)){
  n_mock = 0
  for(j in 1:nrow(WT_pos)){
    if(WT_pos$chr[j] == WT_OMAP$chr[i] & WT_pos$start[j] < WT_OMAP$start[i] & WT_pos$end[j] > WT_OMAP$end[i]){n_mock = n_mock + 1}
  }
  mock_pos = append(mock_pos, n_mock)
}
WT_OMAP$mock_pos = mock_pos
WT_mock_distribution = bin_map

bin_map <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/bin_map_1000000.bed", header=FALSE)
colnames(bin_map) = c("chr", "start", "end", "n")
pos = c()
for(i in 1:nrow(bin_map)){
  n_pos = 0
  for(j in 1:nrow(KO_mock)){
    if(bin_map$chr[i] == KO_mock$chr[j] & bin_map$start[i] < KO_mock$start[j] & bin_map$end[i] > KO_mock$end[j]){n_pos = n_pos + 1}
  }
  pos = append(pos, n_pos)
}
bin_map$pos = pos
ggplot(bin_map) + geom_density(aes(x = log2(pos))) + geom_vline(aes(xintercept = log2(100)), color = "red", linetype = "dashed")
KO_pos = bin_map %>% dplyr::filter(pos > 100)
mock_pos = c()
for(i in 1:nrow(KO_OMAP)){
  n_mock = 0
  for(j in 1:nrow(KO_pos)){
    if(KO_pos$chr[j] == KO_OMAP$chr[i] & KO_pos$start[j] < KO_OMAP$start[i] & KO_pos$end[j] > KO_OMAP$end[i]){n_mock = n_mock + 1}
  }
  mock_pos = append(mock_pos, n_mock)
}
KO_OMAP$mock_pos = mock_pos
KO_mock_distribution = bin_map

WT_OMAP_filtered = WT_OMAP %>% dplyr::filter(mock_pos == 0)
KO_OMAP_filtered = KO_OMAP %>% dplyr::filter(mock_pos == 0)

WT_interacting_domains = GSE234808_WT_Interacting_Domains %>% dplyr::filter(!V1 %in% c("chrY"))
value2 = WT_interacting_domains %>% dplyr::select(V1)
pdf("/home/shared_folder/output_plots/circular_plot_WT_OMAP_density_mock_normalized.pdf", height = 12, width = 12)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(cytoband_df)
for(i in unique(WT_interacting_domains$V1)){
  value_i = value2 %>% dplyr::filter(V1 == i)
  WT_domains_i = WT_interacting_domains %>% dplyr::filter(V1 == i)
  circos.genomicRect(WT_domains_i, value = value_i, col = "black", border = NA, ytop = -20, ybottom = -30, sector.index = i, lty = 5)
}
circos.genomicDensity(list(KO_OMAP_filtered, WT_OMAP_filtered), col = c("white", "#AA0DFE"), border = FALSE)
dev.off()


KO_interacting_domains = GSE234808_KO_Interacting_Domains %>% dplyr::filter(!V1 %in% c("chrY"))
value2 = KO_interacting_domains %>% dplyr::select(V1)
pdf("/home/shared_folder/output_plots/circular_plot_RBM20KO_OMAP_density_mock_normalized.pdf", height = 12, width = 12)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(cytoband_df)
for(i in unique(KO_interacting_domains$V1)){
  value_i = value2 %>% dplyr::filter(V1 == i)
  KO_domains_i = KO_interacting_domains %>% dplyr::filter(V1 == i)
  circos.genomicRect(KO_domains_i, value = value_i, col = "black", border = NA, ytop = -20, ybottom = -25, sector.index = i, lty = 5)
}
circos.genomicDensity(list(WT_OMAP_filtered, KO_OMAP_filtered), col = c("white", "#AA0DFE"), border = FALSE)
dev.off()


##TIDs FISH analysis
Conglomerate_DNA_Hits_TTN_Distances <- read.csv("/home/shared_folder/chromatin_data/FISH_data/Conglomerate_DNA_Hits_TTN_Distances.csv", sep=";")
cis_interactions = Conglomerate_DNA_Hits_TTN_Distances %>% dplyr::filter(cis_trans == "cis")
cis_interactions = cis_interactions %>% dplyr::mutate(interactor = ifelse(interactor == "ctr", "ctr_chr2", interactor)) %>%
  mutate(genotype = ifelse(genotype == "KO", "RBM20KO", genotype))
ggplot(cis_interactions) + stat_ecdf(aes(x = distance, color = interactor, linetype = genotype), linewidth = 0.5) + 
  scale_color_manual(values = c("#AA0DFE", "black")) + scale_linetype_manual(values = c("dotted", "solid")) +
  ylab("cumulative frequency") + xlab("distance (µm)")
ggsave("/home/shared_folder/output_plots/cumulative_CIDs_FISH.pdf", width = 8, height = 7)

cis_interactions$genotype = factor(cis_interactions$genotype, levels = c("WT", "RBM20KO"))
ggplot(cis_interactions) + stat_ecdf(aes(x = distance, color = interactor), linewidth = 0.5) + 
  scale_color_manual(values = c("#AA0DFE", "black")) + facet_wrap(~genotype) +
  ylab("cumulative frequency") + xlab("distance (µm)") + theme(strip.background = element_blank())
ggsave("/home/shared_folder/output_plots/cumulative_CIDs_FISH_separate.pdf", width = 12, height = 7)

#panel d
trans_interactions = Conglomerate_DNA_Hits_TTN_Distances %>% dplyr::filter(cis_trans == "trans")
trans_interactions = trans_interactions %>% dplyr::mutate(interactor = ifelse(interactor == "ctr", "ctr_chr4", interactor)) %>%
  mutate(genotype = ifelse(genotype == "KO", "RBM20KO", genotype))
ggplot(trans_interactions) + stat_ecdf(aes(x = distance, color = interactor, linetype = genotype), linewidth = 0.5) + 
  scale_color_manual(values = c("black", "#FA0087", "#3283FE", "#FBE426", "#1CBE4F", "#0072B2","#D55E00", "#2ED9FF", "#AA0DFE")) + scale_linetype_manual(values = c("dotted", "solid")) +
  ylab("cumulative frequency") + xlab("distance (µm)")
ggsave("/home/shared_folder/panelG/cumulative_TIDs_FISH.pdf", width = 8, height = 7)

trans_interactions$genotype = factor(trans_interactions$genotype, levels = c("WT", "RBM20KO"))
ggplot(trans_interactions) + stat_ecdf(aes(x = distance, color = interactor), linewidth = 0.5) + 
  scale_color_manual(values = c("black", "#FA0087", "#3283FE", "#FBE426", "#1CBE4F", "#0072B2","#D55E00", "#2ED9FF", "#AA0DFE")) + facet_wrap(~genotype) +
  ylab("cumulative frequency") + xlab("distance (µm)") + theme(strip.background = element_blank())
ggsave("/home/shared_folder/panelG/cumulative_TIDs_FISH_separate.pdf", width = 12, height = 7)


##correlation matrices
compartments_undiff_N1_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_undiff_N1_150000.bedGraph", header=FALSE)
colnames(compartments_undiff_N1_150000) = c("chr", "start", "end", "WT_ESCday00_N1")
compartments_undiff_N2_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_undiff_N2_150000.bedGraph", header=FALSE)
colnames(compartments_undiff_N2_150000) = c("chr", "start", "end", "WT_ESCday00_N2")
compartments_meso_N1_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_meso_N1_150000.bedGraph", header=FALSE)
colnames(compartments_meso_N1_150000) = c("chr", "start", "end", "WT_MESday02_N1")
compartments_meso_N2_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_meso_N2_150000.bedGraph", header=FALSE)
colnames(compartments_meso_N2_150000) = c("chr", "start", "end", "WT_MESday02_N2")
compartments_prog_N1_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_prog_N1_150000.bedGraph", header=FALSE)
colnames(compartments_prog_N1_150000) = c("chr", "start", "end", "WT_CPday05_N1")
compartments_prog_N2_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_prog_N2_150000.bedGraph", header=FALSE)
colnames(compartments_prog_N2_150000) = c("chr", "start", "end", "WT_CPday05_N2")
compartments_cardio_N1_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_cardio_N1_150000.bedGraph", header=FALSE)
colnames(compartments_cardio_N1_150000) = c("chr", "start", "end", "WT_CMday14_N1")
compartments_cardio_N2_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_cardio_N2_150000.bedGraph", header=FALSE)
colnames(compartments_cardio_N2_150000) = c("chr", "start", "end", "WT_CMday14_N2")
compartments_WT_D2_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_WT_D2_150000.bedGraph", header=FALSE)
colnames(compartments_WT_D2_150000) = c("chr", "start", "end", "WT_CMday30_D2")
compartments_WT_D3_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_WT_D3_150000.bedGraph", header=FALSE)
colnames(compartments_WT_D3_150000) = c("chr", "start", "end", "WT_CMday30_D3")
compartments_KO_D2_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_KO_D2_150000.bedGraph", header=FALSE)
colnames(compartments_KO_D2_150000) = c("chr", "start", "end", "KO_CMday30_D2")
compartments_KO_D3_150000 <- read.delim("/home/shared_folder/chromatin_data/compartments_data/compartments_KO_D3_150000.bedGraph", header=FALSE)
colnames(compartments_KO_D3_150000) = c("chr", "start", "end", "KO_CMday30_D3")

compartments_all = compartments_undiff_N1_150000 %>% dplyr::left_join(compartments_undiff_N2_150000, by = c("chr", "start", "end")) %>%
  dplyr::left_join(compartments_meso_N1_150000, by = c("chr", "start", "end")) %>% dplyr::left_join(compartments_meso_N2_150000, by = c("chr", "start", "end")) %>%
  dplyr::left_join(compartments_prog_N1_150000, by = c("chr", "start", "end")) %>% dplyr::left_join(compartments_prog_N2_150000, by = c("chr", "start", "end")) %>%
  dplyr::left_join(compartments_cardio_N1_150000, by = c("chr", "start", "end")) %>% dplyr::left_join(compartments_cardio_N2_150000, by = c("chr", "start", "end")) %>%
  dplyr::left_join(compartments_WT_D2_150000, by = c("chr", "start", "end")) %>% dplyr::left_join(compartments_WT_D3_150000, by = c("chr", "start", "end")) %>%
  dplyr::left_join(compartments_KO_D2_150000, by = c("chr", "start", "end")) %>% dplyr::left_join(compartments_KO_D3_150000, by = c("chr", "start", "end"))

compartments_all = na.omit(compartments_all)
compartments_all = compartments_all %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>% dplyr::select(-c(chr, start, end))
rownames(compartments_all) = compartments_all$bin
compartments_all = compartments_all %>% dplyr::select(-bin)

corr = rcorr(as.matrix(compartments_all))
pearson = as.data.frame(corr[["r"]])
annotation_col = data.frame(name = colnames(compartments_all))
genotype = c()
differentiation = c()
batch = c()
for(i in annotation_col$name){
  genotype = append(genotype, strsplit(i, "_")[[1]][1])
  differentiation = append(differentiation, strsplit(i, "_")[[1]][2])
  batch = append(batch, strsplit(i, "_")[[1]][3])
}
annotation_col$genotype = genotype
annotation_col$differentiation = differentiation
annotation_col$batch = batch
rownames(annotation_col) = annotation_col$name
annotation_col = annotation_col %>% dplyr::select(-name)
as.ggplot(pheatmap(pearson, annotation_row = annotation_col, display_numbers = TRUE, color = viridis::viridis(n = 25, option = "H"), number_color = "black"))
ggsave("/home/shared_folder/output_plots/compartments_pearson.pdf")

#SCC
sample_IDs = colnames(pearson)
tab = data.frame(rn = sample_IDs)
cnames = colnames(tab)
combs = c()
names = c()
for(i in sample_IDs){
  scc_i = c()
  sample_i = cool2matrix(paste0("/home/shared_folder/chromatin_data/HiC_matrices/cool/", i, "_150000_iced.cool"))
  for(k in sample_IDs){
    sample_k = cool2matrix(paste0("/home/shared_folder/chromatin_data/HiC_matrices/cool/", k, "_150000_iced.cool"))
    if(nrow(sample_k) > nrow(sample_i)) {
      score = get.scc(sample_i, sample_k, resol = as.integer(150000), h = 5, lbr = 0, ubr = 5000000)
    } else {
      score = get.scc(sample_k, sample_i, resol = as.integer(150000), h = 5, lbr = 0, ubr = 5000000)
    }
    scc = score[["scc"]]
    scc_i = append(scc_i, scc)
    
    if(i != k){
      ordered = as.vector(sort(c(i, k)))
      n1 = ordered[1]
      n2 = ordered[2]
      comb = paste0(n1, "_", n2)
    } else{comb = paste0(i, "_", k)}
    combs = append(combs, scc)
    names = append(names, comb)
  }
  tab$n = as.numeric(scc_i)
  cnames = append(cnames, i)
  colnames(tab) = cnames
}

rownames(tab) = tab$rn
tab = tab %>% dplyr::select(-rn)

as.ggplot(pheatmap(tab, annotation_row = annotation_col, display_numbers = TRUE, color = viridis::viridis(n = 25, option = "H"), number_color = "black"))
ggsave("/home/shared_folder/output_plots/HiCrep_pearson.pdf")

#WT/KO chr2 correlations
WT_OMAP_R1 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/WT_OMAP_R1_binned100000.bed", header=FALSE)
chr = c()
start = c()
end = c()
value = c()
for(i in WT_OMAP_R1$V1){
  chr = append(chr, strsplit(i, ":")[[1]][1])
  start = append(start, strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1])
  end = append(end, strsplit(strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2], "[|]")[[1]][1])
  value = append(value, strsplit(i, "[|]")[[1]][2])
}
O_MAP_WT_R1 = data.frame(chr = chr, start = start, end = end, WT_OMAP_R1 = value) %>% dplyr::filter(chr == "chr2")

WT_OMAP_R2 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/WT_OMAP_R2_binned100000.bed", header=FALSE)
chr = c()
start = c()
end = c()
value = c()
for(i in WT_OMAP_R2$V1){
  chr = append(chr, strsplit(i, ":")[[1]][1])
  start = append(start, strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1])
  end = append(end, strsplit(strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2], "[|]")[[1]][1])
  value = append(value, strsplit(i, "[|]")[[1]][2])
}
O_MAP_WT_R2 = data.frame(chr = chr, start = start, end = end, WT_OMAP_R2 = value) %>% dplyr::filter(chr == "chr2")

WT_mock_R1 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/WT_mock_R1_binned100000.bed", header=FALSE)
chr = c()
start = c()
end = c()
value = c()
for(i in WT_mock_R1$V1){
  chr = append(chr, strsplit(i, ":")[[1]][1])
  start = append(start, strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1])
  end = append(end, strsplit(strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2], "[|]")[[1]][1])
  value = append(value, strsplit(i, "[|]")[[1]][2])
}
mock_WT_R1 = data.frame(chr = chr, start = start, end = end, WT_mock_R1 = value) %>% dplyr::filter(chr == "chr2")

WT_mock_R2 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/WT_mock_R2_binned100000.bed", header=FALSE)
chr = c()
start = c()
end = c()
value = c()
for(i in WT_mock_R2$V1){
  chr = append(chr, strsplit(i, ":")[[1]][1])
  start = append(start, strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1])
  end = append(end, strsplit(strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2], "[|]")[[1]][1])
  value = append(value, strsplit(i, "[|]")[[1]][2])
}
mock_WT_R2 = data.frame(chr = chr, start = start, end = end, WT_mock_R2 = value) %>% dplyr::filter(chr == "chr2")


RBM20KO_OMAP_R1 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/RBM20KO_OMAP_R1_binned100000.bed", header=FALSE)
chr = c()
start = c()
end = c()
value = c()
for(i in RBM20KO_OMAP_R1$V1){
  chr = append(chr, strsplit(i, ":")[[1]][1])
  start = append(start, strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1])
  end = append(end, strsplit(strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2], "[|]")[[1]][1])
  value = append(value, strsplit(i, "[|]")[[1]][2])
}
O_MAP_RBM20KO_R1 = data.frame(chr = chr, start = start, end = end, RBM20KO_OMAP_R1 = value) %>% dplyr::filter(chr == "chr2")

RBM20KO_OMAP_R2 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/RBM20KO_OMAP_R2_binned100000.bed", header=FALSE)
chr = c()
start = c()
end = c()
value = c()
for(i in RBM20KO_OMAP_R2$V1){
  chr = append(chr, strsplit(i, ":")[[1]][1])
  start = append(start, strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1])
  end = append(end, strsplit(strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2], "[|]")[[1]][1])
  value = append(value, strsplit(i, "[|]")[[1]][2])
}
O_MAP_RBM20KO_R2 = data.frame(chr = chr, start = start, end = end, RBM20KO_OMAP_R2 = value) %>% dplyr::filter(chr == "chr2")

RBM20KO_mock_R1 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/RBM20KO_mock_R1_binned100000.bed", header=FALSE)
chr = c()
start = c()
end = c()
value = c()
for(i in RBM20KO_mock_R1$V1){
  chr = append(chr, strsplit(i, ":")[[1]][1])
  start = append(start, strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1])
  end = append(end, strsplit(strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2], "[|]")[[1]][1])
  value = append(value, strsplit(i, "[|]")[[1]][2])
}
mock_RBM20KO_R1 = data.frame(chr = chr, start = start, end = end, RBM20KO_mock_R1 = value) %>% dplyr::filter(chr == "chr2")

RBM20KO_mock_R2 <- read.delim("/home/shared_folder/pulldown_data/chromatin_data/RBM20KO_mock_R2_binned100000.bed", header=FALSE)
chr = c()
start = c()
end = c()
value = c()
for(i in RBM20KO_mock_R2$V1){
  chr = append(chr, strsplit(i, ":")[[1]][1])
  start = append(start, strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1])
  end = append(end, strsplit(strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2], "[|]")[[1]][1])
  value = append(value, strsplit(i, "[|]")[[1]][2])
}
mock_RBM20KO_R2 = data.frame(chr = chr, start = start, end = end, RBM20KO_mock_R2 = value) %>% dplyr::filter(chr == "chr2")


matrix_WT = O_MAP_WT_R1 %>% dplyr::left_join(O_MAP_WT_R2, by = c("chr", "start", "end")) %>%
  dplyr::left_join(mock_WT_R1, by = c("chr", "start", "end")) %>%
  dplyr::left_join(mock_WT_R2, by = c("chr", "start", "end")) 

matrix_KO = O_MAP_RBM20KO_R1 %>% dplyr::left_join(O_MAP_RBM20KO_R2, by = c("chr", "start", "end")) %>%
  dplyr::left_join(mock_RBM20KO_R1, by = c("chr", "start", "end")) %>%
  dplyr::left_join(mock_RBM20KO_R2, by = c("chr", "start", "end"))

matrix_WT = matrix_WT %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>% dplyr::select(-c(chr, start, end))
matrix_KO = matrix_KO %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>% dplyr::select(-c(chr, start, end))

compartments_WT_N2_100000 <- read.csv("/home/shared_folder/chromatin_data/compartments_data/compartments_WT_N2_100000.bedGraph", header=FALSE, comment.char="#", skip = 0, sep = "\t")
compartments_WT_N2_100000 = compartments_WT_N2_100000[2:nrow(compartments_WT_N2_100000),]
colnames(compartments_WT_N2_100000) = c("chr", "start", "end", "WT_N2_compartments")
compartments_WT_N2_100000 = compartments_WT_N2_100000 %>% dplyr::filter(chr == "chr2") %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>% dplyr::select(-c(chr, start, end))

compartments_WT_N3_100000 <- read.csv("/home/shared_folder/chromatin_data/compartments_data/compartments_WT_N3_100000.bedGraph", header=FALSE, comment.char="#", skip = 0, sep = "\t")
compartments_WT_N3_100000 = compartments_WT_N3_100000[2:nrow(compartments_WT_N3_100000),]
colnames(compartments_WT_N3_100000) = c("chr", "start", "end", "WT_N3_compartments")
compartments_WT_N3_100000 = compartments_WT_N3_100000 %>% dplyr::filter(chr == "chr2") %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>% dplyr::select(-c(chr, start, end))

compartments_KO_N2_100000 <- read.csv("/home/shared_folder/chromatin_data/compartments_data/compartments_KO_N2_100000.bedGraph", header=FALSE, comment.char="#", skip = 0, sep = "\t")
compartments_KO_N2_100000 = compartments_KO_N2_100000[2:nrow(compartments_KO_N2_100000),]
colnames(compartments_KO_N2_100000) = c("chr", "start", "end", "RBM20KO_N2_compartments")
compartments_KO_N2_100000 = compartments_KO_N2_100000 %>% dplyr::filter(chr == "chr2") %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>% dplyr::select(-c(chr, start, end))

compartments_KO_N3_100000 <- read.csv("/home/shared_folder/chromatin_data/compartments_data/compartments_KO_N3_100000.bedGraph", header=FALSE, comment.char="#", skip = 0, sep = "\t")
compartments_KO_N3_100000 = compartments_KO_N3_100000[2:nrow(compartments_KO_N3_100000),]
colnames(compartments_KO_N3_100000) = c("chr", "start", "end", "RBM20KO_N3_compartments")
compartments_KO_N3_100000 = compartments_KO_N3_100000 %>% dplyr::filter(chr == "chr2") %>% dplyr::mutate(bin = paste0(chr, "_", start)) %>% dplyr::select(-c(chr, start, end))

virtual4C_WT_N2 <- read.delim("/home/shared_folder/chromatin_data/virtual4C_data/WT_N2_chr2_178806793_178808363.bedGraph", header=FALSE)
colnames(virtual4C_WT_N2) = c("chr", "start", "end", "WT_N2_virtual4C")
virtual4C_WT_N2 = virtual4C_WT_N2 %>% dplyr::filter(chr == "chr2") %>% dplyr::mutate(bin = paste0(chr, "_", as.integer(start - 1))) %>% dplyr::select(-c(chr, start, end))

virtual4C_WT_N3 <- read.delim("/home/shared_folder/chromatin_data/virtual4C_data/WT_N3_chr2_178806793_178808363.bedGraph", header=FALSE)
colnames(virtual4C_WT_N3) = c("chr", "start", "end", "WT_N3_virtual4C")
virtual4C_WT_N3 = virtual4C_WT_N3 %>% dplyr::filter(chr == "chr2") %>% dplyr::mutate(bin = paste0(chr, "_", as.integer(start - 1))) %>% dplyr::select(-c(chr, start, end))

virtual4C_KO_N2 <- read.delim("/home/shared_folder/chromatin_data/virtual4C_data/KO_N2_chr2_178806793_178808363.bedGraph", header=FALSE)
colnames(virtual4C_KO_N2) = c("chr", "start", "end", "KO_N2_virtual4C")
virtual4C_KO_N2 = virtual4C_KO_N2 %>% dplyr::filter(chr == "chr2") %>% dplyr::mutate(bin = paste0(chr, "_", as.integer(start - 1))) %>% dplyr::select(-c(chr, start, end))

virtual4C_KO_N3 <- read.delim("/home/shared_folder/chromatin_data/virtual4C_data/KO_N3_chr2_178806793_178808363.bedGraph", header=FALSE)
colnames(virtual4C_KO_N3) = c("chr", "start", "end", "KO_N3_virtual4C")
virtual4C_KO_N3 = virtual4C_KO_N3 %>% dplyr::filter(chr == "chr2") %>% dplyr::mutate(bin = paste0(chr, "_", as.integer(start - 1))) %>% dplyr::select(-c(chr, start, end))

matrix_WT = matrix_WT %>% 
  dplyr::left_join(compartments_WT_N2_100000, by = "bin") %>% dplyr::left_join(compartments_WT_N3_100000, by = "bin") %>%
  dplyr::left_join(virtual4C_WT_N2, by = "bin") %>% dplyr::left_join(virtual4C_WT_N3, by = "bin") %>%
  dplyr::select(-bin)

matrix_KO = matrix_KO %>% 
  dplyr::left_join(compartments_KO_N2_100000, by = "bin") %>% dplyr::left_join(compartments_KO_N3_100000, by = "bin") %>%
  dplyr::left_join(virtual4C_KO_N2, by = "bin") %>% dplyr::left_join(virtual4C_KO_N3, by = "bin") %>%
  dplyr::select(-bin)

corr = rcorr(as.matrix(matrix_WT))
pearson_WT = as.data.frame(corr[["r"]])
corr = rcorr(as.matrix(matrix_KO))
pearson_KO = as.data.frame(corr[["r"]])

max = max(c(melt(pearson_WT)$value, melt(pearson_KO)$value))
min = min(c(melt(pearson_WT)$value, melt(pearson_KO)$value))
breaks = seq(min, max, length.out = 100)
colors = viridis::turbo(length(breaks) - 1)

annotation_col = data.frame(name = colnames(matrix_WT))
annotation_col$experiment = c("O-MAP", "O-MAP", "mock", "mock", "HiC_compartments", "HiC_compartments", "TTN_virtual4C", "TTN_virtual4C")
rownames(annotation_col) = annotation_col$name
annotation_col = annotation_col %>% dplyr::select(-name)
as.ggplot(pheatmap(pearson_WT, annotation_row = annotation_col, display_numbers = TRUE, color = colors, number_color = "black", breaks = breaks))
ggsave("/home/shared_folder/output_plots/pearson_WT.pdf")

annotation_col = data.frame(name = colnames(matrix_KO))
annotation_col$experiment = c("O-MAP", "O-MAP", "mock", "mock", "HiC_compartments", "HiC_compartments", "TTN_virtual4C", "TTN_virtual4C")
rownames(annotation_col) = annotation_col$name
annotation_col = annotation_col %>% dplyr::select(-name)
as.ggplot(pheatmap(pearson_KO, annotation_row = annotation_col, display_numbers = TRUE, color = colors, number_color = "black", breaks = breaks))
ggsave("/home/shared_folder/output_plots/pearson_RBM20KO.pdf")


##transC
#circular plots
bin_scores_WT_TTN <- read.csv("/home/shared_folder/chromatin_data/transC_data/bin_scores_WT_TTN.csv") %>% dplyr::mutate(n = 1:2887)
bin_scores_RBM20KO_TTN <- read.csv("/home/shared_folder/chromatin_data/transC_data/bin_scores_RBM20KO_TTN.csv") %>% dplyr::mutate(n = 1:2887)
bin_map_1000000 <- read.delim("/home/shared_folder/chromatin_data/transC_data/bin_map_1000000.bed", header=FALSE)

TIDs = c(1235, 1056, 2063, 2882, 1751, 1813, 18, 1674, 2579, 1749, 1629, 1676, 2496)
clique = c(1951, 1917, 1902, 1889, 1671, 804, 616, 2618, 2620, 2581, 2354, 2370)

viridis_col = viridis(2887, option = "D")
col_fun = colorRamp2(c(1:2887), viridis_col)

cytoband = read.cytoband(species = "hg38")
cytoband_df = cytoband$df %>% dplyr::filter(V1 %in% c(bin_map_1000000$V1))

pdf("/home/shared_folder/output_plots/circular_plot_WT_all.pdf", height = 12, width = 12)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(cytoband_df)
region1 = bin_map_1000000 %>% dplyr::filter(V4 != 427) %>% dplyr::select(c(V1, V2, V3))
region2 = do.call("rbind", replicate(length(rownames(region1)), bin_map_1000000 %>% dplyr::filter(V4 == "427"), simplify = FALSE))
ranks = (bin_scores_WT_TTN %>% dplyr::arrange(X0) %>% dplyr::filter(X0 != 427))$n
col = col_fun(ranks)
circos.genomicLink(region1, region2, lwd = 0.1, col =  col)
draw(Legend(col_fun = col_fun, title = "Trans-C ranking", direction = "horizontal"), x = unit(0.1, "npc"), y = unit(0.1, "npc"))
dev.off()

pdf("/home/shared_folder/output_plots/circular_plot_WT_TIDs.pdf", height = 12, width = 12)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(cytoband_df)
region1 = bin_map_1000000 %>% dplyr::filter(V4 %in% TIDs) %>% dplyr::select(c(V1, V2, V3))
region2 = do.call("rbind", replicate(length(rownames(region1)), bin_map_1000000 %>% dplyr::filter(V4 == "427"), simplify = FALSE))
ranks = (bin_scores_WT_TTN %>% dplyr::arrange(X0) %>% dplyr::filter(X0 %in% TIDs))$n
col = col_fun(ranks)
circos.genomicLink(region1, region2, lwd = 0.1, col =  col)
draw(Legend(col_fun = col_fun, title = "Trans-C ranking", direction = "horizontal"), x = unit(0.1, "npc"), y = unit(0.1, "npc"))
dev.off()

pdf("/home/shared_folder/output_plots/circular_plot_WT_clique.pdf", height = 12, width = 12)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(cytoband_df)
region1 = bin_map_1000000 %>% dplyr::filter(V4 %in% clique) %>% dplyr::select(c(V1, V2, V3))
region2 = do.call("rbind", replicate(length(rownames(region1)), bin_map_1000000 %>% dplyr::filter(V4 == "427"), simplify = FALSE))
ranks = (bin_scores_WT_TTN %>% dplyr::arrange(X0) %>% dplyr::filter(X0 %in% clique))$n
col = col_fun(ranks)
circos.genomicLink(region1, region2, lwd = 0.1, col =  col)
draw(Legend(col_fun = col_fun, title = "Trans-C ranking", direction = "horizontal"), x = unit(0.1, "npc"), y = unit(0.1, "npc"))
dev.off()

pdf("/home/shared_folder/output_plots/circular_plot_RBM20KO_all.pdf", height = 12, width = 12)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(cytoband_df)
region1 = bin_map_1000000 %>% dplyr::filter(V4 != 427) %>% dplyr::select(c(V1, V2, V3))
region2 = do.call("rbind", replicate(length(rownames(region1)), bin_map_1000000 %>% dplyr::filter(V4 == "427"), simplify = FALSE))
ranks = (bin_scores_RBM20KO_TTN %>% dplyr::arrange(X0) %>% dplyr::filter(X0 != 427))$n
col = col_fun(ranks)
circos.genomicLink(region1, region2, lwd = 0.1, col =  col)
draw(Legend(col_fun = col_fun, title = "Trans-C ranking", direction = "horizontal"), x = unit(0.1, "npc"), y = unit(0.1, "npc"))
dev.off()

pdf("/home/shared_folder/output_plots/circular_plot_RBM20KO_TIDs.pdf", height = 12, width = 12)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(cytoband_df)
region1 = bin_map_1000000 %>% dplyr::filter(V4 %in% TIDs) %>% dplyr::select(c(V1, V2, V3))
region2 = do.call("rbind", replicate(length(rownames(region1)), bin_map_1000000 %>% dplyr::filter(V4 == "427"), simplify = FALSE))
ranks = (bin_scores_RBM20KO_TTN %>% dplyr::arrange(X0) %>% dplyr::filter(X0 %in% TIDs))$n
col = col_fun(ranks)
circos.genomicLink(region1, region2, lwd = 0.1, col =  col)
draw(Legend(col_fun = col_fun, title = "Trans-C ranking", direction = "horizontal"), x = unit(0.1, "npc"), y = unit(0.1, "npc"))
dev.off()

pdf("/home/shared_folder/output_plots/circular_plot_RBM20_clique.pdf", height = 12, width = 12)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(cytoband_df)
region1 = bin_map_1000000 %>% dplyr::filter(V4 %in% clique) %>% dplyr::select(c(V1, V2, V3))
region2 = do.call("rbind", replicate(length(rownames(region1)), bin_map_1000000 %>% dplyr::filter(V4 == "427"), simplify = FALSE))
ranks = (bin_scores_RBM20KO_TTN %>% dplyr::arrange(X0) %>% dplyr::filter(X0 %in% clique))$n
col = col_fun(ranks)
circos.genomicLink(region1, region2, lwd = 0.1, col =  col)
draw(Legend(col_fun = col_fun, title = "Trans-C ranking", direction = "horizontal"), x = unit(0.1, "npc"), y = unit(0.1, "npc"))
dev.off()

#AUC
bin_scores_WT = bin_scores_WT_TTN %>% dplyr::filter(! X0 %in% c("427")) %>%
  dplyr::mutate(TIDs = ifelse(X0 %in% TIDs, "TIDs", "")) %>% dplyr::mutate(clique = ifelse(X0 %in% clique, "clique", ""))
bin_scores_WT = bin_scores_WT %>% dplyr::mutate(n = rev(1:length(rownames(bin_scores_WT))))
bin_scores_KO = bin_scores_RBM20KO_TTN %>% dplyr::filter(! X0 %in% c("427")) %>%
  dplyr::mutate(TIDs = ifelse(X0 %in% TIDs, "TIDs", "")) %>% dplyr::mutate(clique = ifelse(X0 %in% clique, "clique", ""))
bin_scores_KO = bin_scores_KO %>% dplyr::mutate(n = rev(1:length(rownames(bin_scores_KO))))

table_random_rankings_WT <- read.csv("/home/shared_folder/chromatin_data/transC_data/table_random_rankings_WT_TTN.csv")
table_random_rankings_WT = as.data.frame(t(table_random_rankings_WT %>% dplyr::select(-X)))

random_WT_TIDs = data.frame(FPR = c(), Sensitivity = c(), n = c())
random_WT_clique = data.frame(FPR = c(), Sensitivity = c(), n = c())
auroc_WT_TIDs = c()
auroc_WT_clique = c()

for(i in 1:1000){
  list_i = data.frame(id = table_random_rankings_WT[, i]) %>%
    dplyr::mutate(TIDs = ifelse(id %in% TIDs, "TIDs", "")) %>% dplyr::mutate(clique = ifelse(id %in% clique, "clique", ""))
  list_i = list_i %>% dplyr::mutate(n = rev(1:length(rownames(list_i))))
  PRROC = roc.curve(scores.class0 = (list_i %>% dplyr::filter(TIDs == "TIDs"))$n,
                    scores.class1 = (list_i %>% dplyr::filter(TIDs != "TIDs"))$n,
                    sorted = TRUE,
                    curve = TRUE)
  curve_i = as.data.frame(PRROC[["curve"]]) %>% dplyr::mutate(n = i)
  colnames(curve_i) = c("FPR", "Sensitivity", "n")
  random_WT_TIDs = rbind(random_WT_TIDs, curve_i)
  auroc_WT_TIDs = append(auroc_WT_TIDs, PRROC[["auc"]])
  
  PRROC = roc.curve(scores.class0 = (list_i %>% dplyr::filter(clique == "clique"))$n,
                    scores.class1 = (list_i %>% dplyr::filter(clique != "clique"))$n,
                    sorted = TRUE,
                    curve = TRUE)
  curve_i = as.data.frame(PRROC[["curve"]]) %>% dplyr::mutate(n = i)
  colnames(curve_i) = c("FPR", "Sensitivity", "n")
  random_WT_clique = rbind(random_WT_clique, curve_i)
  auroc_WT_clique = append(auroc_WT_clique, PRROC[["auc"]])
}
random_WT_TIDs = random_WT_TIDs %>% dplyr::select(c("FPR", "Sensitivity", "n"))
max_WT = random_WT_TIDs %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.95))[1])
min_WT = random_WT_TIDs %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.05))[1])
sum_WT = random_WT_TIDs %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.5))[1])

PRROC_obj_WT_TIDs <- roc.curve(scores.class0 = (bin_scores_WT %>% dplyr::filter(TIDs == "TIDs"))$n,
                          scores.class1 = (bin_scores_WT %>% dplyr::filter(TIDs != "TIDs"))$n,
                          sorted = TRUE,
                          curve = TRUE)
curve_WT_TIDs = as.data.frame(PRROC_obj_WT_TIDs[["curve"]])
colnames(curve_WT_TIDs) = c("FPR", "Sensitivity", "score")
curve_WT_TIDs = curve_WT_TIDs %>% dplyr::mutate(ranking = nrow(curve_WT_TIDs) - score)

empirical_ROC = data.frame(ROC = auroc_WT_TIDs)
n_up = 0
n_down = 0
for(i in auroc_WT_TIDs){
  if(i > PRROC_obj_WT_TIDs[["auc"]]){n_up = n_up + 1} else{n_down = n_down + 1}
}
ggplot() + geom_line(aes(x = curve_WT_TIDs$FPR, y = curve_WT_TIDs$Sensitivity, color = curve_WT_TIDs$ranking), size = 2) + ggtitle(paste0("AUC = ", signif(PRROC_obj_WT_TIDs[["auc"]], digits = 3), "\np-val empirical = ", n_up/1000)) +
  scale_color_viridis(name = "Trans-C ranking", guide = guide_colourbar(reverse = TRUE)) + geom_line(aes(x = sum_WT$FPR, y = sum_WT$Sensitivity)) + geom_line(aes(x = min_WT$FPR, y = min_WT$Sensitivity), linetype = "dashed") +
  geom_line(aes(x = max_WT$FPR, y = max_WT$Sensitivity), linetype = "dashed") + xlab("FPR (false positive rate)") + ylab("sensitivity") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/trans-C_AUC_WT_TIDs.pdf")

random_WT_clique = random_WT_clique %>% dplyr::select(c("FPR", "Sensitivity", "n"))
max_WT = random_WT_clique %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.95))[1])
min_WT = random_WT_clique %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.05))[1])
sum_WT = random_WT_clique %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.5))[1])

PRROC_obj_WT_clique <- roc.curve(scores.class0 = (bin_scores_WT %>% dplyr::filter(clique == "clique"))$n,
                               scores.class1 = (bin_scores_WT %>% dplyr::filter(clique != "clique"))$n,
                               sorted = TRUE,
                               curve = TRUE)
curve_WT_clique = as.data.frame(PRROC_obj_WT_clique[["curve"]])
colnames(curve_WT_clique) = c("FPR", "Sensitivity", "score")
curve_WT_clique = curve_WT_clique %>% dplyr::mutate(ranking = nrow(curve_WT_clique) - score)

empirical_ROC = data.frame(ROC = auroc_WT_clique)
n_up = 0
n_down = 0
for(i in auroc_WT_clique){
  if(i > PRROC_obj_WT_clique[["auc"]]){n_up = n_up + 1} else{n_down = n_down + 1}
}
ggplot() + geom_line(aes(x = curve_WT_clique$FPR, y = curve_WT_clique$Sensitivity, color = curve_WT_clique$ranking), size = 2) + ggtitle(paste0("AUC = ", signif(PRROC_obj_WT_clique[["auc"]], digits = 3), "\np-val empirical = ", n_up/1000)) +
  scale_color_viridis(name = "Trans-C ranking", guide = guide_colourbar(reverse = TRUE)) + geom_line(aes(x = sum_WT$FPR, y = sum_WT$Sensitivity)) + geom_line(aes(x = min_WT$FPR, y = min_WT$Sensitivity), linetype = "dashed") +
  geom_line(aes(x = max_WT$FPR, y = max_WT$Sensitivity), linetype = "dashed") + xlab("FPR (false positive rate)") + ylab("sensitivity") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/trans-C_AUC_WT_clique.pdf")


table_random_rankings_KO <- read.csv("/home/shared_folder/chromatin_data/transC_data/table_random_rankings_KO_TTN.csv")
table_random_rankings_KO = as.data.frame(t(table_random_rankings_KO %>% dplyr::select(-X)))

random_KO_TIDs = data.frame(FPR = c(), Sensitivity = c(), n = c())
random_KO_clique = data.frame(FPR = c(), Sensitivity = c(), n = c())
auroc_KO_TIDs = c()
auroc_KO_clique = c()

for(i in 1:1000){
  list_i = data.frame(id = table_random_rankings_KO[, i]) %>%
    dplyr::mutate(TIDs = ifelse(id %in% TIDs, "TIDs", "")) %>% dplyr::mutate(clique = ifelse(id %in% clique, "clique", ""))
  list_i = list_i %>% dplyr::mutate(n = rev(1:length(rownames(list_i))))
  PRROC = roc.curve(scores.class0 = (list_i %>% dplyr::filter(TIDs == "TIDs"))$n,
                    scores.class1 = (list_i %>% dplyr::filter(TIDs != "TIDs"))$n,
                    sorted = TRUE,
                    curve = TRUE)
  curve_i = as.data.frame(PRROC[["curve"]]) %>% dplyr::mutate(n = i)
  colnames(curve_i) = c("FPR", "Sensitivity", "n")
  random_KO_TIDs = rbind(random_KO_TIDs, curve_i)
  auroc_KO_TIDs = append(auroc_KO_TIDs, PRROC[["auc"]])
  
  PRROC = roc.curve(scores.class0 = (list_i %>% dplyr::filter(clique == "clique"))$n,
                    scores.class1 = (list_i %>% dplyr::filter(clique != "clique"))$n,
                    sorted = TRUE,
                    curve = TRUE)
  curve_i = as.data.frame(PRROC[["curve"]]) %>% dplyr::mutate(n = i)
  colnames(curve_i) = c("FPR", "Sensitivity", "n")
  random_KO_clique = rbind(random_KO_clique, curve_i)
  auroc_KO_clique = append(auroc_KO_clique, PRROC[["auc"]])
}
random_KO_TIDs = random_KO_TIDs %>% dplyr::select(c("FPR", "Sensitivity", "n"))
max_KO = random_KO_TIDs %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.95))[1])
min_KO = random_KO_TIDs %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.05))[1])
sum_KO = random_KO_TIDs %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.5))[1])

PRROC_obj_KO_TIDs <- roc.curve(scores.class0 = (bin_scores_KO %>% dplyr::filter(TIDs == "TIDs"))$n,
                               scores.class1 = (bin_scores_KO %>% dplyr::filter(TIDs != "TIDs"))$n,
                               sorted = TRUE,
                               curve = TRUE)
curve_KO_TIDs = as.data.frame(PRROC_obj_KO_TIDs[["curve"]])
colnames(curve_KO_TIDs) = c("FPR", "Sensitivity", "score")
curve_KO_TIDs = curve_KO_TIDs %>% dplyr::mutate(ranking = nrow(curve_KO_TIDs) - score)

empirical_ROC = data.frame(ROC = auroc_KO_TIDs)
n_up = 0
n_down = 0
for(i in auroc_KO_TIDs){
  if(i > PRROC_obj_KO_TIDs[["auc"]]){n_up = n_up + 1} else{n_down = n_down + 1}
}
ggplot() + geom_line(aes(x = curve_KO_TIDs$FPR, y = curve_KO_TIDs$Sensitivity, color = curve_KO_TIDs$ranking), size = 2) + ggtitle(paste0("AUC = ", signif(PRROC_obj_KO_TIDs[["auc"]], digits = 3), "\np-val empirical = ", n_up/1000)) +
  scale_color_viridis(name = "Trans-C ranking", guide = guide_colourbar(reverse = TRUE)) + geom_line(aes(x = sum_KO$FPR, y = sum_KO$Sensitivity)) + geom_line(aes(x = min_KO$FPR, y = min_KO$Sensitivity), linetype = "dashed") +
  geom_line(aes(x = max_KO$FPR, y = max_KO$Sensitivity), linetype = "dashed") + xlab("FPR (false positive rate)") + ylab("sensitivity") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/trans-C_AUC_RBM20KO_TIDs.pdf")

random_KO_clique = random_KO_clique %>% dplyr::select(c("FPR", "Sensitivity", "n"))
max_KO = random_KO_clique %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.95))[1])
min_KO = random_KO_clique %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.05))[1])
sum_KO = random_KO_clique %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.5))[1])

PRROC_obj_KO_clique <- roc.curve(scores.class0 = (bin_scores_KO %>% dplyr::filter(clique == "clique"))$n,
                                 scores.class1 = (bin_scores_KO %>% dplyr::filter(clique != "clique"))$n,
                                 sorted = TRUE,
                                 curve = TRUE)
curve_KO_clique = as.data.frame(PRROC_obj_KO_clique[["curve"]])
colnames(curve_KO_clique) = c("FPR", "Sensitivity", "score")
curve_KO_clique = curve_KO_clique %>% dplyr::mutate(ranking = nrow(curve_KO_clique) - score)

empirical_ROC = data.frame(ROC = auroc_KO_clique)
n_up = 0
n_down = 0
for(i in auroc_KO_clique){
  if(i > PRROC_obj_KO_clique[["auc"]]){n_up = n_up + 1} else{n_down = n_down + 1}
}
ggplot() + geom_line(aes(x = curve_KO_clique$FPR, y = curve_KO_clique$Sensitivity, color = curve_KO_clique$ranking), size = 2) + ggtitle(paste0("AUC = ", signif(PRROC_obj_KO_clique[["auc"]], digits = 3), "\np-val empirical = ", n_up/1000)) +
  scale_color_viridis(name = "Trans-C ranking", guide = guide_colourbar(reverse = TRUE)) + geom_line(aes(x = sum_KO$FPR, y = sum_KO$Sensitivity)) + geom_line(aes(x = min_KO$FPR, y = min_KO$Sensitivity), linetype = "dashed") +
  geom_line(aes(x = max_KO$FPR, y = max_KO$Sensitivity), linetype = "dashed") + xlab("FPR (false positive rate)") + ylab("sensitivity") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(10,0,10,0)))
ggsave("/home/shared_folder/output_plots/trans-C_AUC_RBM20KO_clique.pdf")

random_KO_eCLIP_HC = random_KO_eCLIP_HC %>% dplyr::select(c("FPR", "Sensitivity", "n"))
max_KO = random_KO_eCLIP_HC %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.95))[1])
min_KO = random_KO_eCLIP_HC %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.05))[1])
sum_KO = random_KO_eCLIP_HC %>% dplyr::group_by(n) %>% dplyr::summarise(FPR = mean(FPR), Sensitivity = quantile(Sensitivity, probs = c(.5))[1])
