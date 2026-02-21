#------------------------------------------------------------------------------
#BULK RNA SEQ
#------------------------------------------------------------------------------

#packages
library(tidyr)
library(stringr)
library(grid)
library(viridis)
library(ggthemes)
library(dplyr)
library(tidyverse)
library(edgeR)
library(compareGroups)
library(dbscan)
library(reshape2)
library(scran)
library(fgsea)
library(RNAseqQC)
library(DESeq2)
library(ensembldb)
library(tibble)
library(ggplotify)
library(limma)
library(Glimma)
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggrepel)
library(ggsignif)
library(igraph)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(WGCNA)
library(readr)

#set up environment
theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 45, hjust = 1)))
options(stringsAsFactors = FALSE)
update_geom_defaults("point", aes(size = 4))
set.seed(1234)

pal = c("#3283FE", "#FA0087", "#009E73", "#FBE426", "#56B4E9", "#FEAF16", "#DEA0FD", "#1CBE4F", "#F6222E", "#1CFFCE", "#325A9B", "#AA0DFE","#D55E00", "#2ED9FF", "#f0E442", "#1C8356", "#0072B2", "#CC79A7")

current_time = format(Sys.time(), "%y%m%d")

args = commandArgs(trailingOnly = TRUE)
table = args[1]
output_folder = args[2]
dir.create(output_folder)

alignment = TRUE
mode = "paired"
genome_folder = "/home/shared_folder/genome/"

#table
print("uploading table file")

table = read_excel(table)
if(alignment == TRUE){
  if(mode == "single"){
    if(!"fastq" %in% colnames(table)){
      print("fastq column is missing"); stop()
    } else{
      for(i in table$fastq){
        if(file.exists(i) == FALSE){
          print("incorrect path in fastq column"); stop()
        }
      }
    }
  } else if(mode == "paired"){
    if(!"fastq1" %in% colnames(table) | !"fastq2" %in% colnames(table)){
      print("fastq1 or fastq2 columns are missing"); stop()
    } else{
      for(i in table$fastq1){if(file.exists(i) == FALSE){
        print("incorrect path in fastq1 column"); stop()}
      }
      for(i in table$fastq2){if(file.exists(i) == FALSE){
        print("incorrect path in fastq2 column"); stop()}
      }
    }
  }
}

#alignment
if(alignment == TRUE){
  print("starting alignment")
  
  dir.create(paste0(output_folder, "/counts"));
  
  if(mode == "single"){
    for(i in 1:length(rownames(table))){
      sample_ID = table$sample_ID[i];
      zipped = ifelse(rev(strsplit(table$fastq[i], "\\.")[[1]])[1] == "gz", " --readFilesCommand zcat", "")
      counts_dir = paste0(output_folder, "/counts/", sample_ID);
      dir.create(counts_dir);
      system2("chmod", paste0("777 ", counts_dir));
      if(trimming == FALSE){fastq = table$fastq[i]} else if(trimming == TRUE){
        fastq = paste0(out_dir, "/", sample_ID, "_trimmed.fq.gz")};
      system2("/STAR", paste0("--runThreadN ", CPU, " --genomeDir ", genome_dir, " --readFilesIn ", fastq, " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts", zipped, " --outFileNamePrefix ", counts_dir, "/", sample_ID))}
  } else if(mode == "paired"){
    for(i in 1:length(rownames(table))){sample_ID = table$sample_ID[i];
    zipped1 = ifelse(rev(strsplit(table$fastq1[i], ".")[[1]])[1] == "gz", " --readFilesCommand zcat", "")
    zipped2 = ifelse(rev(strsplit(table$fastq2[i], ".")[[1]])[1] == "gz", " --readFilesCommand zcat", "")
    if(zipped1 == zipped2){zipped = zipped1} else{
      stop("fastq files must be provided in the same zip format for R1 and R2")
    }
    counts_dir = paste0(output_folder, "/counts/", sample_ID);
    if(trimming == FALSE){fastq1 = table$fastq1[i]; fastq2 = table$fastq2[i]} else if(trimming == TRUE){
      fastq1 = paste0(out_dir, "/", sample_ID, "_val_1.fq.gz");
      fastq2 = paste0(out_dir, "/", sample_ID, "_val_2.fq.gz")};
    
    system2("/STAR", paste0("--runThreadN ", CPU, " --genomeDir ", genome_dir, " --readFilesIn ", fastq1, " ", fastq2, " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts", zipped, " --outFileNamePrefix ", counts_dir, "/", sample_ID))}
  }
  print("finishing alignment")
}
