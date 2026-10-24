library(dplyr)
library(ggplot2)
library(reshape2)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 45, hjust = 1)))

#input
WT_input1 <- read.delim("/home/shared_folder/output_STAR_input/bed/WT_input1.bed", header=FALSE) %>% dplyr::arrange(V2)
wt_I1 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  wt_i = WT_input1 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(wt_i)) == 1){wt_I1 = rbind(wt_I1, data.frame(position = c(i), val = c(wt_i$V5[1])))}
  else{wt_I1 = rbind(wt_I1, data.frame(position = c(i), val = c(0)))}
}

WT_input2 <- read.delim("/home/shared_folder/output_STAR_input/bed/WT_input2.bed", header=FALSE)
wt_I2 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  wt_i = WT_input2 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(wt_i)) == 1){wt_I2 = rbind(wt_I2, data.frame(position = c(i), val = c(wt_i$V5[1])))}
  else{wt_I2 = rbind(wt_I2, data.frame(position = c(i), val = c(0)))}
}

WT_input3 <- read.delim("/home/shared_folder/output_STAR_input/bed/WT_input3.bed", header=FALSE)
wt_I3 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  wt_i = WT_input3 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(wt_i)) == 1){wt_I3 = rbind(wt_I3, data.frame(position = c(i), val = c(wt_i$V5[1])))}
  else{wt_I3 = rbind(wt_I3, data.frame(position = c(i), val = c(0)))}
}

WT_input4 <- read.delim("/home/shared_folder/output_STAR_input/bed/WT_input4.bed", header=FALSE)
wt_I4 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  wt_i = WT_input4 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(wt_i)) == 1){wt_I4 = rbind(wt_I4, data.frame(position = c(i), val = c(wt_i$V5[1])))}
  else{wt_I4 = rbind(wt_I4, data.frame(position = c(i), val = c(0)))}
}

WT_input_joined = wt_I1 %>% dplyr::left_join(wt_I2, by = "position") %>% dplyr::left_join(wt_I3, by = "position") %>% dplyr::left_join(wt_I4, by = "position")
colnames(WT_input_joined) = c("position", "WT_input1", "WT_input2", "WT_input3", "WT_input4")
  
KO_input1 <- read.delim("/home/shared_folder/output_STAR_input/bed/KO_input1.bed", header=FALSE) %>% dplyr::arrange(V2)
ko_I1 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  ko_i = KO_input1 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(ko_i)) == 1){ko_I1 = rbind(ko_I1, data.frame(position = c(i), val = c(ko_i$V5[1])))}
  else{ko_I1 = rbind(ko_I1, data.frame(position = c(i), val = c(0)))}
}

KO_input2 <- read.delim("/home/shared_folder/output_STAR_input/bed/KO_input2.bed", header=FALSE)
ko_I2 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  ko_i = KO_input2 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(ko_i)) == 1){ko_I2 = rbind(ko_I2, data.frame(position = c(i), val = c(ko_i$V5[1])))}
  else{ko_I2 = rbind(ko_I2, data.frame(position = c(i), val = c(0)))}
}

KO_input3 <- read.delim("/home/shared_folder/output_STAR_input/bed/KO_input3.bed", header=FALSE)
ko_I3 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  ko_i = KO_input3 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(ko_i)) == 1){ko_I3 = rbind(ko_I3, data.frame(position = c(i), val = c(ko_i$V5[1])))}
  else{ko_I3 = rbind(ko_I3, data.frame(position = c(i), val = c(0)))}
}

KO_input4 <- read.delim("/home/shared_folder/output_STAR_input/bed/KO_input4.bed", header=FALSE)
ko_I4 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  ko_i = KO_input4 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(ko_i)) == 1){ko_I4 = rbind(ko_I4, data.frame(position = c(i), val = c(ko_i$V5[1])))}
  else{ko_I4 = rbind(ko_I4, data.frame(position = c(i), val = c(0)))}
}

KO_input_joined = ko_I1 %>% dplyr::left_join(ko_I2, by = "position") %>% dplyr::left_join(ko_I3, by = "position") %>% dplyr::left_join(ko_I4, by = "position")
colnames(KO_input_joined) = c("position", "KO_input1", "KO_input2", "KO_input3", "KO_input4")

#pulldown
WT_pulldown1 <- read.delim("/home/shared_folder/output_STAR_pulldown/bed/WT_pulldown1.bed", header=FALSE) %>% dplyr::arrange(V2)
wt_P1 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  wt_p = WT_pulldown1 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(wt_p)) == 1){wt_P1 = rbind(wt_P1, data.frame(position = c(i), val = c(wt_p$V5[1])))}
  else{wt_P1 = rbind(wt_P1, data.frame(position = c(i), val = c(0)))}
}

WT_pulldown2 <- read.delim("/home/shared_folder/output_STAR_pulldown/bed/WT_pulldown2.bed", header=FALSE) %>% dplyr::arrange(V2)
wt_P2 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  wt_p = WT_pulldown2 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(wt_p)) == 1){wt_P2 = rbind(wt_P2, data.frame(position = c(i), val = c(wt_p$V5[1])))}
  else{wt_P2 = rbind(wt_P2, data.frame(position = c(i), val = c(0)))}
}

WT_pulldown3 <- read.delim("/home/shared_folder/output_STAR_pulldown/bed/WT_pulldown3.bed", header=FALSE) %>% dplyr::arrange(V2)
wt_P3 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  wt_p = WT_pulldown3 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(wt_p)) == 1){wt_P3 = rbind(wt_P3, data.frame(position = c(i), val = c(wt_p$V5[1])))}
  else{wt_P3 = rbind(wt_P3, data.frame(position = c(i), val = c(0)))}
}

WT_pulldown4 <- read.delim("/home/shared_folder/output_STAR_pulldown/bed/WT_pulldown4.bed", header=FALSE) %>% dplyr::arrange(V2)
wt_P4 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  wt_p = WT_pulldown4 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(wt_p)) == 1){wt_P4 = rbind(wt_P4, data.frame(position = c(i), val = c(wt_p$V5[1])))}
  else{wt_P4 = rbind(wt_P4, data.frame(position = c(i), val = c(0)))}
}

WT_pulldown_joined = wt_P1 %>% dplyr::left_join(wt_P2, by = "position") %>% dplyr::left_join(wt_P3, by = "position") %>% dplyr::left_join(wt_P4, by = "position")
colnames(WT_pulldown_joined) = c("position", "WT_pulldown1", "WT_pulldown2", "WT_pulldown3", "WT_pulldown4")

KO_pulldown1 <- read.delim("/home/shared_folder/output_STAR_pulldown/bed/KO_pulldown1.bed", header=FALSE) %>% dplyr::arrange(V2)
ko_P1 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  ko_p = KO_pulldown1 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(ko_p)) == 1){ko_P1 = rbind(ko_P1, data.frame(position = c(i), val = c(ko_p$V5[1])))}
  else{ko_P1 = rbind(ko_P1, data.frame(position = c(i), val = c(0)))}
}

KO_pulldown2 <- read.delim("/home/shared_folder/output_STAR_pulldown/bed/KO_pulldown2.bed", header=FALSE) %>% dplyr::arrange(V2)
ko_P2 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  ko_p = KO_pulldown2 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(ko_p)) == 1){ko_P2 = rbind(ko_P2, data.frame(position = c(i), val = c(ko_p$V5[1])))}
  else{ko_P2 = rbind(ko_P2, data.frame(position = c(i), val = c(0)))}
}

KO_pulldown3 <- read.delim("/home/shared_folder/output_STAR_pulldown/bed/KO_pulldown3.bed", header=FALSE) %>% dplyr::arrange(V2)
ko_P3 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  ko_p = KO_pulldown3 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(ko_p)) == 1){ko_P3 = rbind(ko_P3, data.frame(position = c(i), val = c(ko_p$V5[1])))}
  else{ko_P3 = rbind(ko_P3, data.frame(position = c(i), val = c(0)))}
}

KO_pulldown4 <- read.delim("/home/shared_folder/output_STAR_pulldown/bed/KO_pulldown4.bed", header=FALSE) %>% dplyr::arrange(V2)
ko_P4 = data.frame(position = c(), val = c())
for(i in 178525989:178830802){
  print(i)
  ko_p = KO_pulldown4 %>% dplyr::filter(V2 <= i & V3 >= i)
  if(length(rownames(ko_p)) == 1){ko_P4 = rbind(ko_P4, data.frame(position = c(i), val = c(ko_p$V5[1])))}
  else{ko_P4 = rbind(ko_P4, data.frame(position = c(i), val = c(0)))}
}

KO_pulldown_joined = ko_P1 %>% dplyr::left_join(ko_P2, by = "position") %>% dplyr::left_join(ko_P3, by = "position") %>% dplyr::left_join(ko_P4, by = "position")
colnames(KO_pulldown_joined) = c("position", "KO_pulldown1", "KO_pulldown2", "KO_pulldown3", "KO_pulldown4")

#normalize for library size
rownames(WT_input_joined) = WT_input_joined$position
WT_input_joined = WT_input_joined %>% dplyr::select(-position)
WT_input_joined = as.matrix(WT_input_joined) + 1
WT_input_joined = as.data.frame(WT_input_joined/colSums(WT_input_joined)*nrow(WT_input_joined))
WT_input_joined$position = as.integer(rownames(WT_input_joined))

rownames(WT_pulldown_joined) = WT_pulldown_joined$position
WT_pulldown_joined = WT_pulldown_joined %>% dplyr::select(-position)
WT_pulldown_joined = as.matrix(WT_pulldown_joined) + 1
WT_pulldown_joined = as.data.frame(WT_pulldown_joined/colSums(WT_pulldown_joined)*nrow(WT_input_joined))
WT_pulldown_joined$position = as.integer(rownames(WT_pulldown_joined))

rownames(KO_input_joined) = KO_input_joined$position
KO_input_joined = KO_input_joined %>% dplyr::select(-position)
KO_input_joined = as.matrix(KO_input_joined) + 1
KO_input_joined = as.data.frame(KO_input_joined/colSums(KO_input_joined)*nrow(WT_input_joined))
KO_input_joined$position = as.integer(rownames(KO_input_joined))

rownames(KO_pulldown_joined) = KO_pulldown_joined$position
KO_pulldown_joined = KO_pulldown_joined %>% dplyr::select(-position)
KO_pulldown_joined = as.matrix(KO_pulldown_joined) + 1
KO_pulldown_joined = as.data.frame(KO_pulldown_joined/colSums(KO_pulldown_joined)*nrow(WT_input_joined))
KO_pulldown_joined$position = as.integer(rownames(KO_pulldown_joined))

#ratio over input
dir.create("/home/shared_folder/ratio_pulldown_input")
WT_ratio = data.frame(position = WT_input_joined$position, WT1_logFC = log2((WT_pulldown_joined$WT_pulldown1) / (WT_input_joined$WT_input1)),
                      WT2_logFC = log2((WT_pulldown_joined$WT_pulldown2) / (WT_input_joined$WT_input2)),
                      WT3_logFC = log2((WT_pulldown_joined$WT_pulldown3) / (WT_input_joined$WT_input3)),
                      WT4_logFC = log2((WT_pulldown_joined$WT_pulldown4) / (WT_input_joined$WT_input4)))

KO_ratio = data.frame(position = KO_input_joined$position, KO1_logFC = log2((KO_pulldown_joined$KO_pulldown1) / (KO_input_joined$KO_input1)),
                      KO2_logFC = log2((KO_pulldown_joined$KO_pulldown2) / (KO_input_joined$KO_input2)),
                      KO3_logFC = log2((KO_pulldown_joined$KO_pulldown3) / (KO_input_joined$KO_input3)),
                      KO4_logFC = log2((KO_pulldown_joined$KO_pulldown4) / (KO_input_joined$KO_input4)))

WT_1 = data.frame(chr = rep(c("2"), times = nrow(WT_input_joined)), start = WT_input_joined$position, end = WT_input_joined$position + 1, WT1_logFC = log2((WT_pulldown_joined$WT_pulldown1) / (WT_input_joined$WT_input1)))
write.table(WT_1, "/home/shared_folder/ratio_pulldown_input/WT1.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") 
WT_2 = data.frame(chr = rep(c("2"), times = nrow(WT_input_joined)), start = WT_input_joined$position, end = WT_input_joined$position + 1, WT2_logFC = log2((WT_pulldown_joined$WT_pulldown2) / (WT_input_joined$WT_input2)))
write.table(WT_2, "/home/shared_folder/ratio_pulldown_input/WT2.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
WT_3 = data.frame(chr = rep(c("2"), times = nrow(WT_input_joined)), start = WT_input_joined$position, end = WT_input_joined$position + 1, WT3_logFC = log2((WT_pulldown_joined$WT_pulldown3) / (WT_input_joined$WT_input3)))
write.table(WT_3, "/home/shared_folder/ratio_pulldown_input/WT3.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") 
WT_4 = data.frame(chr = rep(c("2"), times = nrow(WT_input_joined)), start = WT_input_joined$position, end = WT_input_joined$position + 1, WT4_logFC = log2((WT_pulldown_joined$WT_pulldown4) / (WT_input_joined$WT_input4)))
write.table(WT_4, "/home/shared_folder/ratio_pulldown_input/WT4.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") 

KO_1 = data.frame(chr = rep(c("2"), times = nrow(KO_input_joined)), start = KO_input_joined$position, end = KO_input_joined$position + 1, KO1_logFC = log2((KO_pulldown_joined$KO_pulldown1) / (KO_input_joined$KO_input1)))
write.table(KO_1, "/home/shared_folder/ratio_pulldown_input/KO1.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") 
KO_2 = data.frame(chr = rep(c("2"), times = nrow(KO_input_joined)), start = KO_input_joined$position, end = KO_input_joined$position + 1, KO2_logFC = log2((KO_pulldown_joined$KO_pulldown2) / (KO_input_joined$KO_input2)))
write.table(KO_2, "/home/shared_folder/ratio_pulldown_input/KO2.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
KO_3 = data.frame(chr = rep(c("2"), times = nrow(KO_input_joined)), start = KO_input_joined$position, end = KO_input_joined$position + 1, KO3_logFC = log2((KO_pulldown_joined$KO_pulldown3) / (KO_input_joined$KO_input3)))
write.table(KO_3, "/home/shared_folder/ratio_pulldown_input/KO3.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") 
KO_4 = data.frame(chr = rep(c("2"), times = nrow(KO_input_joined)), start = KO_input_joined$position, end = KO_input_joined$position + 1, KO4_logFC = log2((KO_pulldown_joined$KO_pulldown4) / (KO_input_joined$KO_input4)))
write.table(KO_4, "/home/shared_folder/ratio_pulldown_input/KO4.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") 

WT = data.frame(chr = rep(c("2"), times = nrow(WT_input_joined)), start = WT_input_joined$position, end = WT_input_joined$position + 1, WT_logFC = (WT_1$WT1_logFC + WT_2$WT2_logFC + WT_3$WT3_logFC + WT_4$WT4_logFC)/4)
write.table(WT, "/home/shared_folder/ratio_pulldown_input/WT.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") 
KO = data.frame(chr = rep(c("2"), times = nrow(KO_input_joined)), start = WT_input_joined$position, end = WT_input_joined$position + 1, KO_logFC = (KO_1$KO1_logFC + KO_2$KO2_logFC + KO_3$KO3_logFC + KO_4$KO4_logFC)/4)
write.table(KO, "/home/shared_folder/ratio_pulldown_input/KO.bedGraph", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") 



