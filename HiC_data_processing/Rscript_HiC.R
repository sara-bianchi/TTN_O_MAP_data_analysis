#------------------------------------------------------------------------------
#HiC pipeline
#------------------------------------------------------------------------------

#packages
library(tidyr)
library(stringr)
library(viridis)
library(ggthemes)
library(dplyr)
library(tidyverse)
library(reshape2)
library(scran)
library(tibble)
library(ggplotify)
library(topGO)
library(ggrepel)
library(readxl)
library(pheatmap)
library(ggplot2)
library(WGCNA)
library(readr)
library(reticulate)
library(Hmisc)
library(factoextra)
library(circlize)
library(Rlab)
library(patchwork)
library(cowplot)
library(hicrep)
library(limma)

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
dir.create("/home/shared_folder/output_HiC")
system2("chmod", paste0("-R 777 /home/shared_folder/output_HiC"))

#upload settings
print("uploading settings")

settings = as.data.frame(read_excel("/home/shared_folder/settings.xlsx", col_names = FALSE))
rownames(settings) = settings$...1
settings = settings %>% dplyr::select(-c("...1", "...2"))
settings = as.data.frame(t(settings))

alignment = TRUE
indexing = TRUE
aligner = "HiC-Pro"
targeted_analysis = TRUE
untargeted_analysis = TRUE

if(untargeted_analysis == TRUE | targeted_analysis == TRUE){analysis = TRUE} else {analysis = FALSE}

if(length(unique(na.omit(settings$genome))) != 1){
  stop("invalid argument for genome")
} else{genome = settings$genome[1]}

if(length(unique(na.omit(settings$chr_sizes))) != 1 | file.exists(settings$chr_sizes[1]) == FALSE){
  stop("invalid argument for chr_sizes")} else{chr_sizes_path = settings$chr_sizes[1]}
chrom_sizes = read.delim(chr_sizes_path, header = FALSE)
colnames(chrom_sizes) = c("chr", "size")

binning = c()
if(alignment == TRUE){
  for(i in unique(na.omit(settings$binning))){if(as.integer(i) == "NA"){
    stop("invalid argument for binning")}
    else{binning = append(binning, as.integer(i))}
  }
}

if(untargeted_analysis == TRUE){
  res_untargeted = c()
  for(i in unique(na.omit(settings$res_untargeted))){if(as.integer(i) == "NA"){
    stop("invalid argument for res_untargeted")}
    else{if(! i %in% binning & alignment == TRUE){
      stop("res_untargeted must be among the binning sizes selected")
    } else{
      res_untargeted = append(res_untargeted, as.integer(i))}
    }
  }
  if(length(unique(na.omit(settings$PC))) != 1 | as.integer(settings$PC[1]) == "NA"){
    stop("invalid argument for PC")
  } else{PC = as.integer(settings$PC[1])}
}

if(targeted_analysis == TRUE){
  if(length(unique(na.omit(settings$res_v4C))) != 0){
    virtual4C = TRUE
    res_v4C = c()
    for(i in unique(na.omit(settings$res_v4C))){if(as.integer(i) == "NA"){
      stop("invalid argument for res_v4C")}
      else {if(! i %in% binning & alignment == TRUE){
        stop("res_v4C must be among the binning sizes selected")
      } else {
        res_v4C = append(res_v4C, as.integer(i))}
      }
    }
    
    viewpoints = c()
    chr = c()
    start = c()
    end = c()
    for(i in unique(na.omit(settings$viewpoint))){
      chr_i = strsplit(i, ":")[[1]][1]
      if(!chr_i %in% chrom_sizes$chr){
        stop("viewpoint not in chr_sizes")
      } else {chr = append(chr, chr_i)
      size_i = (chrom_sizes %>% dplyr::filter(chr == chr_i))$size}
      start_i = strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][1]
      end_i = strsplit(strsplit(i, ":")[[1]][2], "-")[[1]][2]

      if(as.integer(start_i) == "NA" | as.integer(end_i) == "NA"){
        stop("invalid viewpoint")
      } else{
        if(as.integer(start_i) > as.integer(end_i) | as.integer(end_i) >= as.integer(size_i)){
          stop("viewpoint coordinates must be in chr:start-end format")
        } else {
          start = append(start, start_i)
          end = append(end, end_i)
        }
      }
    }
    vp = data.frame(chr, start, end)
  } else {
    virtual4C = FALSE
  }
  
} else {
  virtual4C = FALSE
}

if(length(unique(na.omit(settings$CPU))) != 1 | as.integer(settings$CPU[1]) == "NA"){
  stop("invalid argument for CPU")
} else{CPU = as.integer(settings$CPU[1])}

if(alignment == TRUE){
  if(length(unique(na.omit(settings$RE))) == 0){
    RE = c("Dnase")
  }
}

res = c()
if(untargeted_analysis == TRUE){res = append(res, res_untargeted)}
if(targeted_analysis == TRUE){
  if(virtual4C == TRUE){
    res = append(res, res_v4C)
  }
}

res = unique(res)


#upload table
table = read_excel("/home/shared_folder/table.xlsx")

if(alignment == FALSE){
  if(normalized == FALSE){
    normalization = TRUE
  } else {normalization = FALSE
  }
  
  if(!"bias" %in% colnames(table)){
    print("normalized matrix in non HiC-pro format provided. Loops analysis will be skipped")
  }
} else {normalization = TRUE}

if(alignment == TRUE){
  if(!"fastq1" %in% colnames(table) | !"fastq2" %in% colnames(table)){
    stop("fastq1 or fastq2 columns are missing")
  } else{
    for(i in table$fastq1){if(file.exists(i) == FALSE){
      stop("incorrect path in fastq1 column")}
    }
    for(i in table$fastq2){if(file.exists(i) == FALSE){
      stop("incorrect path in fastq2 column")}
    }
  }
} else {
  if(!"matrix" %in% colnames(table)) {
    stop("matrix column is missing")
  } else{
    formats = c("homer", "hic", "h5", "mcool", "cool", "matrix")
    f = c()
    for(i in table$matrix){
      if(file.exists(i) == FALSE){
        stop("incorrect path in matrix column")}
      format = tail(strsplit(i, "[.]")[[1]], n = 1)
      if(!format %in% formats){
        stop("invalid matrix format; supported formats: .homer, .hic, .h5, .mcool, .cool, .matrix")
      } else{
        f = append(f, format)
      }
    }
    table$format = f
    if("matrix" %in% f){
      for(i in f){if(i != "matrix"){
        stop("different formats are not accepted with hicpro")
      }}
      if(length(res) != 1){
        stop("only one resolution size is possible with hicpro format")
      }
      if(!"bed" %in% colnames(table)) {
        stop(".matrix format requires bed column")
      } else{
        for(i in 1:length(rownames(table))){
          if(file.exists(table$bed[i]) == FALSE){
            stop("incorrect path in bed column")
          }
        }
        if("bias" %in% colnames(table)){
          for(i in 1:length(rownames(table))){
            if(file.exists(table$bias[i]) == FALSE){
              stop("incorrect path in bias column")
            }
          }
        }
      }
    }
  }
}

if(untargeted_analysis == TRUE){
  references = unique(na.omit(settings$reference))
  for(i in references){
    if(! i %in% table$condition){
      stop("reference condition not in table conditions")
    }
  }
}

if(alignment == FALSE){
  if(untargeted_analysis == FALSE & targeted_analysis == TRUE){
    if(virtual4C == FALSE){
      normalization = FALSE
    }
  }
}

#genome indexing
if(normalization == TRUE){
  if(length(unique(na.omit(settings$indexing))) != 1 | as.logical(settings$indexing[1]) == "NA"){
    stop("invalid argument for indexing")
  } else{indexing = as.logical(settings$indexing[1])}

  if(length(unique(na.omit(settings$fasta))) != 1 | file.exists(settings$fasta[1]) == FALSE){
    stop("invalid argument for fasta")} else{fasta = settings$fasta[1]}
  }
  
#indexing
if(indexing == TRUE){
  if(aligner == "HiC-pro"){
    print("starting bowtie index generation")
      
    dir.create(paste0("/home/shared_folder/ooutput_HiC/bowtie_index/"))
    dir.create(paste0("/home/shared_folder/ooutput_HiC/bowtie_index/", genome))
    genome_folder = paste0("/home/shared_folder/ooutput_HiC/bowtie_index/")
    system2("bowtie2-build", paste0(fasta, " ", genome))
    system2("mv", paste0(genome, "* ", genome_folder))
  }
}


#alignment
if(alignment == TRUE){
  
  if(aligner == "HiC-pro"){
    #set up config-hicpro.txt file
    print("starting HiC-pro alignment")
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/HiC_pro/"))
    system2("cp", paste0("/HiC-Pro_3.1.0/config-hicpro.txt", " /home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    system2("chmod", paste0("777 /home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    
    if(RE == "Dnase"){
      system2("sed", paste0("-i 's|GENOME_FRAGMENT = HindIII_resfrag_hg19.bed|GENOME_FRAGMENT =|g' /home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
      system2("sed", paste0("-i 's|LIGATION_SITE = AAGCTAGCTT|LIGATION_SITE =|g' /home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    } else{
      enzymes = ""
      enzymes2 = ""
      for(i in RE){
        enzymes = paste0(enzymes, tolower(i), " ")
        enzymes2 =paste0(enzymes2, i, "_")
      }
      system2("/HiC-Pro_3.1.0/bin/utils/digest_genome.py", paste0("-r ", enzymes, "-o /home/shared_folder/ooutput_HiC/HiC_pro/", genome, "_", enzymes2, "fragments.bed ", fasta))
      system2("sed", paste0("-i 's|LIGATION_SITE = AAGCTAGCTT|LIGATION_SITE =|g' /home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    }
    
    system2("sed", paste0("-i 's|N_CPU = 2|N_CPU = ", as.character(CPU), "|g' ", "/home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    system2("sed", paste0("-i 's|BOWTIE2_IDX_PATH =|BOWTIE2_IDX_PATH = ", genome_folder, "|g' ", "/home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    system2("sed", paste0("-i 's|REFERENCE_GENOME = hg19|REFERENCE_GENOME = ", genome, "|g' ", "/home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    system2("sed", paste0("-i 's|GENOME_SIZE = chrom_hg19.sizes|GENOME_SIZE = ", chr_sizes_path, "|g' ", "/home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    
    bins = ""
    for(i in binning){
      bins = paste0(bins, i, " ")
    }
    
    system2("sed", paste0("-i 's|BIN_SIZE = 20000 40000 150000 500000 1000000|BIN_SIZE = ", bins, "|g' ", "/home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
    
    #set up directories with files
    dir.create(paste0("/home/shared_folder/ooutput_HiC/HiC_pro/input"))
    
    for(i in 1:length(rownames(table))){
      dir.create(paste0("/home/shared_folder/ooutput_HiC/HiC_pro/input/", table$sample_ID[i]))
      system2("mv", paste0(table$fastq1[i], " /home/shared_folder/ooutput_HiC/HiC_pro/input/", table$sample_ID[i], "/", table$sample_ID[i], "_R1.fastq.gz"))
      system2("mv", paste0(table$fastq2[i], " /home/shared_folder/ooutput_HiC/HiC_pro/input/", table$sample_ID[i], "/", table$sample_ID[i], "_R2.fastq.gz"))
    }
    
    #run HiC-pro
    system2("/HiC-Pro_3.1.0/bin/HiC-Pro", paste0("-i /home/shared_folder/ooutput_HiC/HiC_pro/input/ -o /home/shared_folder/ooutput_HiC/HiC_pro/output -c /home/shared_folder/ooutput_HiC/HiC_pro/config-hicpro_", current_time, ".txt"))
  }
}

#ANALYSIS - preprocessing
if(analysis == TRUE){
  
  #file conversion
  
  samples = table$sample_ID
  
  if(alignment == TRUE){if(aligner == "HiC-pro"){
    for(j in res){
      path_mcool = c()
      path_bias = c()
      for(i in 1:length(rownames(table))){
        matrix_i = paste0("/home/shared_folder/ooutput_HiC/HiC_pro/output/hic_results/matrix/", table$sample_ID[i], "/iced/", as.character(j), "/", table$sample_ID[i], "_", as.character(j), "_iced.matrix")
        bed_i = paste0("/home/shared_folder/ooutput_HiC/HiC_pro/output/hic_results/matrix/", table$sample_ID[i], "/raw/", as.character(j), "/", table$sample_ID[i], "_", as.character(j), "_abs.bed")
        bias_i = paste0("/home/shared_folder/ooutput_HiC/HiC_pro/output/hic_results/matrix/", table$sample_ID[i], "/iced/", as.character(j), "/", table$sample_ID[i], "_", as.character(j), "_iced.matrix.biases")
        path_matrix = append(path_matrix, matrix_i)
        path_bed = append(path_bed, bed_i)
        path_bias = append(path_bias, bias_i)
      }
      
      paths = data.frame(samples, path_matrix, path_bed, path_bias)
      colnames(paths) = c("sample_ID", paste0("iced_matrix_", as.character(j)), paste0("bed_", as.character(j)), paste0("bias_", as.character(j)))
      table = table %>% dplyr::left_join(paths, by = "sample_ID")
    }
  }} else {
    
    #hic_pro
    if("bed" %in% colnames(table)){
      j = res[1]
      
      if(normalized == TRUE){
        if("bias" %in% colnames(table)){
          paths = data.frame(samples, table$matrix, table$bed, table$bias)
          colnames(paths) = c("sample_ID", paste0("iced_matrix_", as.character(j)), paste0("bed_", as.character(j)), paste0("bias_", as.character(j)))
          table = table %>% dplyr::left_join(paths, by = "sample_ID")
        } else {
          paths = data.frame(samples, table$matrix, table$bed)
          colnames(paths) = c("sample_ID", paste0("iced_matrix_", as.character(j)), paste0("bed_", as.character(j)))
          table = table %>% dplyr::left_join(paths, by = "sample_ID")
        }
        
      } else {
        
        paths = data.frame(samples, table$matrix, table$bed)
        colnames(paths) = c("sample_ID", paste0("raw_matrix_", as.character(j)), paste0("bed_", as.character(j)))
        table = table %>% dplyr::left_join(paths, by = "sample_ID")
        
        dir.create(paste0("/home/shared_folder/ooutput_HiC/mcool_raw"))
        path_mcool = c()
        for(i in 1:length(rownames(table))){
          matrix_i = table$matrix[i]
          bed_i = table$bed[i]
          out_i = paste0("/home/shared_folder/ooutput_HiC/mcool_raw/", table$sample_ID[i], ".mcool")
          
          system2("hicConvertFormat", paste0("--matrices ", matrix_i, " --bedFileHicpro ", bed_i, " --inputFormat hicpro --outputFormat mcool --outFileName ", out_i))
          path_mcool = append(path_mcool, out_i)
        }
        
        paths = data.frame(samples, path_mcool)
        colnames(paths) = c("sample_ID", "raw_mcool")
        table = table %>% dplyr::left_join(paths, by = "sample_ID")
      }
      
    } else {
      
      is_hic = 0
      for(i in table$format){if(i != "hic") {is_hic = is_hic + 1}}
      
      if(is_hic == 0){
        
        #hic to HiC-pro
        dir.create(paste0("/home/shared_folder/ooutput_HiC/format_conversion"))
        
        out = paste0("/home/shared_folder/ooutput_HiC/format_conversion/")
        paths = table %>% dplyr::select(c("sample_ID", "matrix"))
        path_matrix = c()
        path_bed = c()
        system2("conda", "create -y -n hictk -c conda-forge -c bioconda hictk")
        for(i in 1:length(rownames(paths))){
          conda_run2(
            envname = "hictk",
            cmd_line = paste0("hictk convert ", paths$matrix[i], " ", out, paths$sample_ID[i], ".mcool"),
            echo = TRUE)
          
          path_matrix = append(path_matrix, paste0(out, paths$sample_ID[i], ".mcool"))
        }
        
        if(normalized == TRUE){
          paths = data.frame(samples, path_matrix)
          colnames(paths) = c("sample_ID", "iced_mcool")
          table = table %>% dplyr::left_join(paths, by = "sample_ID")
        } else {
          paths = data.frame(samples, path_matrix)
          colnames(paths) = c("sample_ID", "raw_mcool")
          table = table %>% dplyr::left_join(paths, by = "sample_ID")
        }
        
      } else {
        is_cool = 0
        for(i in table$format){if(!i %in% c("mcool", "cool")) {is_cool = is_cool + 1}}
        
        
        #mcool/cool
        if(is_cool == 0){
          path_mcool = c()
          for(i in table$matrix){path_mcool = append(path_mcool, i)}
          
          if(normalized == TRUE){
            paths = data.frame(samples, path_mcool)
            colnames(paths) = c("sample_ID", "iced_mcool")
            table = table %>% dplyr::left_join(paths, by = "sample_ID")
          } else {
            paths = data.frame(samples, path_mcool)
            colnames(paths) = c("sample_ID", "raw_mcool")
            table = table %>% dplyr::left_join(paths, by = "sample_ID")
          }
          
        } else {
          
          #others (hic/homer)
          
          if(normalized == TRUE){
            dir.create(paste0("/home/shared_folder/ooutput_HiC/mcool_iced"))
            path_mcool = c()
            for(i in 1:length(rownames(table))){
              matrix_i = table$matrix[i]
              out_i = paste0("/home/shared_folder/ooutput_HiC/mcool_iced/", table$sample_ID[i], ".mcool")
              
              system2("hicConvertFormat", paste0("--matrices ", matrix_i, " --inputFormat ", table$format[i], " --outputFormat mcool --outFileName ", out_i))
              path_mcool = append(path_mcool, out_i)
            }
            paths = data.frame(samples, path_mcool)
            colnames(paths) = c("sample_ID", "iced_mcool")
            table = table %>% dplyr::left_join(paths, by = "sample_ID")
            
          } else {
            dir.create(paste0("/home/shared_folder/ooutput_HiC/mcool_raw"))
            path_mcool = c()
            for(i in 1:length(rownames(table))){
              matrix_i = table$matrix[i]
              out_i = paste0("/home/shared_folder/ooutput_HiC/mcool_raw/", table$sample_ID[i], ".mcool")
              
              system2("hicConvertFormat", paste0("--matrices ", matrix_i, " --inputFormat ", table$format[i], " --outputFormat mcool --outFileName ", out_i))
              path_mcool = append(path_mcool, out_i)
            }
            paths = data.frame(samples, path_mcool)
            colnames(paths) = c("sample_ID", "raw_mcool")
            table = table %>% dplyr::left_join(paths, by = "sample_ID")
          }
        }
      }
      
      #mcool/cool to HiC-pro
      dir.create(paste0("/home/shared_folder/ooutput_HiC/format_conversion"))
      for(j in res){
        dir.create(paste0("/home/shared_folder/ooutput_HiC/format_conversion/", as.character(j)))
        out_j = paste0("/home/shared_folder/ooutput_HiC/format_conversion/", as.character(j), "/")
        
        if(normalized == TRUE){
          paths = table %>% dplyr::select(c("sample_ID", "iced_mcool"))
          path_matrix = c()
          path_bed = c()
          colnames(paths) = c("sample_ID", "iced_mcool")
          for(i in 1:length(rownames(paths))){
            system2("python", paste0("/preprocess.py -input cool -file ", paths$iced_mcool[i], " -res ", as.character(j), " -prefix ", table$sample_ID[i], " -genomeFile ", chr_sizes_path))
            system2("mv", paste0(paths$sample_ID[i], "_", as.character(j), ".matrix ", out_j))
            system2("mv", paste0(paths$sample_ID[i], "_", as.character(j), "_abs.bed ", out_j))
            path_matrix = append(path_matrix, paste0(out_j, paths$sample_ID[i], "_", as.character(j), ".matrix"))
            path_bed = append(path_bed, paste0(out_j, paths$sample_ID[i], "_", as.character(j), "_abs.bed"))
          }
          
          paths = data.frame(samples, path_matrix, path_bed)
          colnames(paths) = c("sample_ID", paste0("iced_matrix_", as.character(j)), paste0("bed_", as.character(j)))
          table = table %>% dplyr::left_join(paths, by = "sample_ID")
          
        } else {
          paths = table %>% dplyr::select(c("sample_ID", "raw_mcool"))
          path_matrix = c()
          path_bed = c()
          colnames(paths) = c("sample_ID", "raw_mcool")
          for(i in 1:length(rownames(paths))){
            system2("python", paste0("/preprocess.py -input cool -file ", paths$raw_mcool[i], " -res ", as.character(j), " -prefix ", table$sample_ID[i], " -genomeFile ", chr_sizes_path))
            system2("mv", paste0(paths$sample_ID[i], "_", as.character(j), ".matrix ", out_j))
            system2("mv", paste0(paths$sample_ID[i], "_", as.character(j), "_abs.bed ", out_j))
            path_matrix = append(path_matrix, paste0(out_j, paths$sample_ID[i], "_", as.character(j), ".matrix"))
            path_bed = append(path_bed, paste0(out_j, paths$sample_ID[i], "_", as.character(j), "_abs.bed"))
          }
          
          paths = data.frame(samples, path_matrix, path_bed)
          colnames(paths) = c("sample_ID", paste0("raw_matrix_", as.character(j)), paste0("bed_", as.character(j)))
          table = table %>% dplyr::left_join(paths, by = "sample_ID")
        }
      }
    }
  }
  
  if(normalization == TRUE){
    
    #HiC-pro normalization
    dir.create(paste0("/home/shared_folder/ooutput_HiC/ICE_normalization"))
    
    for(j in res){
      dir.create(paste0("/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j)))
      
      system2("cp", paste0("/HiC-Pro_3.1.0/config-hicpro.txt", " /home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      system2("chmod", paste0("777 /home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      
      system2("sed", paste0("-i 's|GENOME_FRAGMENT = HindIII_resfrag_hg19.bed|GENOME_FRAGMENT =|g' /home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      system2("sed", paste0("-i 's|LIGATION_SITE = AAGCTAGCTT|LIGATION_SITE =|g' /home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      system2("sed", paste0("-i 's|N_CPU = 2|N_CPU = ", as.character(CPU), "|g' ", "/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      system2("sed", paste0("-i 's|BOWTIE2_IDX_PATH =|BOWTIE2_IDX_PATH = ", genome_folder, "|g' ", "/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      system2("sed", paste0("-i 's|REFERENCE_GENOME = hg19|REFERENCE_GENOME = ", genome, "|g' ", "/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      system2("sed", paste0("-i 's|GENOME_SIZE = chrom_hg19.sizes|GENOME_SIZE = ", chr_sizes_path, "|g' ", "/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      system2("sed", paste0("-i 's|BIN_SIZE = 20000 40000 150000 500000 1000000|BIN_SIZE = ", as.character(j), "|g' ", "/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt"))
      
      #set up directories with files
      dir.create(paste0("/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/input"))
      
      paths = table %>% dplyr::select(c("sample_ID", paste0("raw_matrix_", as.character(j))))
      colnames(paths) = c("sample_ID", "raw_matrix")
      for(i in 1:length(rownames(paths))){
        dir.create(paste0("/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/input/", paths$sample_ID[i]))
        system2("cp", paste0(paths$raw_matrix[i], " /home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/input/", paths$sample_ID[i], "/"))
      }
      
      system2("/HiC-Pro_3.1.0/bin/HiC-Pro", paste0("-i /home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/input -o /home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/output -c /home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/config-hicpro_", current_time, ".txt -s ice_norm"))
      
      path_matrix = c()
      path_bias = c()
      for(i in 1:length(rownames(paths))){
        matrix_i = paste0("/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/output/hic_results/matrix/", table$sample_ID[i], "/iced/", as.character(j), "/", table$sample_ID[i], "_", as.character(j), "_iced.matrix")
        bias_i = paste0("/home/shared_folder/ooutput_HiC/ICE_normalization/", as.character(j), "/output/hic_results/matrix/", table$sample_ID[i], "/iced/", as.character(j), "/", table$sample_ID[i], "_", as.character(j), "_iced.matrix.biases")
        path_matrix = append(path_matrix, matrix_i)
        path_bias = append(path_bias, bias_i)
      }
      paths = data.frame(samples, path_matrix, path_bias)
      colnames(paths) = c("sample_ID", paste0("iced_matrix_", as.character(j)), paste0("bias_", as.character(j)))
      table = table %>% dplyr::left_join(paths, by = "sample_ID")
    }
  }
  
  if(normalization == TRUE | "bias" %in% colnames(table) | normalized == TRUE){
    
    #iced mcool conversion
    dir.create(paste0("/home/shared_folder/ooutput_HiC/mcool_iced"))
    
    for(j in res){
      paths = table %>% dplyr::select(c("sample_ID", paste0("iced_matrix_", as.character(j)), paste0("bed_", as.character(j))))
      colnames(paths) = c("sample_ID", "iced_matrix", "bed")
      dir.create(paste0("/home/shared_folder/ooutput_HiC/mcool_iced/", as.character(j)))
      path_mcool = c()
      for(i in 1:length(rownames(paths))){
        matrix_i = paths$iced_matrix[i]
        bed_i = paths$bed[i]
        out_i = paste0("/home/shared_folder/ooutput_HiC/mcool_iced/", as.character(j), "/", paths$sample_ID[i], "_", as.character(j), "_iced.mcool")
        system2("hicConvertFormat", paste0("--matrices ", matrix_i, " --bedFileHicpro ", bed_i, " --inputFormat hicpro --outputFormat mcool --outFileName ", out_i))
        path_mcool = append(path_mcool, out_i)
      }
      
      paths = data.frame(samples, path_mcool)
      colnames(paths) = c("sample_ID", paste0("iced_mcool_", as.character(j)))
      table = table %>% dplyr::left_join(paths, by = "sample_ID")
    }
  }
  
  if(untargeted_analysis == TRUE){
    
    #iced h5 conversion
    dir.create(paste0("/home/shared_folder/ooutput_HiC/h5_iced"))
    dir.create(paste0("/home/shared_folder/ooutput_HiC/cool_iced"))
    
    for(j in res_untargeted){
      paths = table %>% dplyr::select(c("sample_ID", paste0("iced_matrix_", as.character(j)), paste0("bed_", as.character(j))))
      colnames(paths) = c("sample_ID", "iced_matrix", "bed")
      dir.create(paste0("/home/shared_folder/ooutput_HiC/h5_iced/", as.character(j)))
      path_h5 = c()
      for(i in 1:length(rownames(paths))){
        matrix_i = paths$iced_matrix[i]
        bed_i = paths$bed[i]
        out_i = paste0("/home/shared_folder/ooutput_HiC/h5_iced/", as.character(j), "/", paths$sample_ID[i], "_", as.character(j), "_iced.h5")
        system2("hicConvertFormat", paste0("--matrices ", matrix_i, " --bedFileHicpro ", bed_i, " --inputFormat hicpro --outputFormat h5 --outFileName ", out_i))
        path_h5 = append(path_h5, out_i)
      }
      
      paths = data.frame(samples, path_h5)
      colnames(paths) = c("sample_ID", paste0("iced_h5_", as.character(j)))
      table = table %>% dplyr::left_join(paths, by = "sample_ID")
    }
    
    for(j in res_untargeted){
      paths = table %>% dplyr::select(c("sample_ID", paste0("iced_matrix_", as.character(j)), paste0("bed_", as.character(j))))
      colnames(paths) = c("sample_ID", "iced_matrix", "bed")
      dir.create(paste0("/home/shared_folder/ooutput_HiC/cool_iced/", as.character(j)))
      path_cool = c()
      for(i in 1:length(rownames(paths))){
        matrix_i = paths$iced_matrix[i]
        bed_i = paths$bed[i]
        out_i = paste0("/home/shared_folder/ooutput_HiC/cool_iced/", as.character(j), "/", paths$sample_ID[i], "_", as.character(j), "_iced.cool")
        system2("hicConvertFormat", paste0("--matrices ", matrix_i, " --bedFileHicpro ", bed_i, " --inputFormat hicpro --outputFormat cool --outFileName ", out_i))
        path_cool = append(path_cool, out_i)
      }
      
      paths = data.frame(samples, path_cool)
      colnames(paths) = c("sample_ID", paste0("iced_cool_", as.character(j)))
      table = table %>% dplyr::left_join(paths, by = "sample_ID")
    }
  }
  
}


#ANALYSIS - untargeted
if(untargeted_analysis == TRUE){
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/TADs"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/autocorrelation"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/ratio"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/HiCrep"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/dcHiC"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/diff_loops"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/IGV"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/gsea"))
  dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/genomeTracks"))
  
  for(j in res_untargeted){
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/TADs/", as.character(j)))
    TADs_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/TADs/", as.character(j))
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/autocorrelation/", as.character(j)))
    autocorr_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/autocorrelation/", as.character(j))
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/ratio/", as.character(j)))
    ratio_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/ratio/", as.character(j))
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/HiCrep/", as.character(j)))
    HiCrep_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/HiCrep/", as.character(j))
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/dcHiC/", as.character(j)))
    dcHiC_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/dcHiC/", as.character(j))
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/diff_loops/", as.character(j)))
    diff_loops_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/diff_loops/", as.character(j))
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/IGV/", as.character(j)))
    IGV_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/IGV/", as.character(j), "/")
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/gsea/", as.character(j)))
    gsea_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/gsea/", as.character(j), "/")
    
    dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/genomeTracks/", as.character(j)))
    track_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/genomeTracks/", as.character(j), "/")
    
    #TADs
    paths = table %>% dplyr::select(c("sample_ID", paste0("iced_h5_", as.character(j))))
    colnames(paths) = c("sample_ID", "iced_h5")
    for(i in 1:length(rownames(table))){
      prefix_i = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/TADs/", as.character(j), "/", paths$sample_ID[i], "_", as.character(j))
      system2("hicFindTADs", paste0("--matrix ", paths$iced_h5[i], " --correctForMultipleTesting fdr --outPrefix ", prefix_i))
    }
    
    #autocorrelation
    paths = table %>% dplyr::select("sample_ID", paste0("iced_mcool_", as.character(j)))
    colnames(paths) = c("sample_ID", "iced_mcool")
    
    for(c in 1:length(rownames(chrom_sizes))){
      dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/autocorrelation/", as.character(j), "/", chrom_sizes$chr[c]))
    }
    for(i in 1:length(rownames(paths))){
      hic = HiCExperiment::import(paths$iced_mcool[i], format = "mcool", resolution = j)
      gis = InteractionSet::interactions(hic)
      gis$score = HiCExperiment::scores(hic, 'count')
      cnts = gis |> 
        tibble::as_tibble() |> 
        dplyr::relocate(c(seqnames1, seqnames2))
      cnts_dup = cnts |> 
        dplyr::rename(seqnames1 = seqnames2, seqnames2 = seqnames1) |> 
        dplyr::relocate(c(seqnames1, seqnames2))
      cnts = rbind(cnts, cnts_dup)
      res = cnts |> 
        dplyr::group_by(seqnames1, seqnames2) |> 
        dplyr::summarize(n = sum(score, na.rm = TRUE)) |> 
        dplyr::mutate(type = ifelse(seqnames1 == seqnames2, 'cis', 'trans')) |> 
        dplyr::group_by(seqnames1, type) |> 
        dplyr::rename(chr = seqnames1) |>
        dplyr::summarize(n = sum(n, na.rm = TRUE)) |> 
        tidyr::pivot_wider(names_from = type, values_from = n)
      
      if(! "trans" %in% colnames(res)){res$trans = rep(c(0), n = length(rownames(res)))}
      ct = res %>% dplyr::mutate(
        n_total = sum(cis + trans, na.rm = TRUE), 
        cis_pct = cis/n_total, 
        trans_pct = trans/n_total
      )
      
      ct = ct %>% dplyr::select(c("chr", "cis_pct", "trans_pct"))
      ct = melt(ct)
      ggplot(ct, aes(x = chr, y = value*100, fill = variable)) + 
        geom_bar(stat = "identity") +
        labs(x = 'Chromosomes', y = '% of cis  and trans contacts')
      ggsave(paste0(autocorr_path, "/ct_perc_", as.character(j), "_", paths$sample_ID[i], ".pdf"))
      
      for(c in 1:length(rownames(chrom_sizes))){
        chr = paste0(chrom_sizes$chr[c], ":1-", chrom_sizes$size[c])
        hic_i = subsetByOverlaps(hic, GRanges(chr), type = "within")
        
        hic_i = HiContacts::autocorrelate(hic_i)
        p1 = HiContacts::plotMatrix(hic_i,
                                    use.scores = "balanced",
                                    scale = "log10")
        
        p2 = HiContacts::plotMatrix(
          hic_i, 
          use.scores = 'autocorrelated', 
          scale = 'linear', 
          cmap = HiContacts::bgrColors(),
          caption = FALSE
        )
        wrap_plots(p1, p2, ncol = 2, nrow = 1,widths = c(10, 10))
        ggsave(paste0(autocorr_path, "/", chrom_sizes$chr[c], "/autocorrelation_", as.character(j), "_", paths$sample_ID[i], ".pdf"), width = 20, height = 10)
      }
    }
    rm(hic, hic_i)
    
    #ratio
    paths = table %>% dplyr::select("sample_ID", "condition", "batch", paste0("iced_mcool_", as.character(j)))
    colnames(paths) = c("sample_ID", "condition", "batch", "iced_mcool")
    
    for(c in 1:length(rownames(chrom_sizes))){
      dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/ratio/", as.character(j), "/", chrom_sizes$chr[c]))
    }
    
    for(k in references){
      paths_k = paths %>% dplyr::filter(condition == k)
      for(i in unique(paths$condition)){
        if(i != k){
          paths_i = paths %>% dplyr::filter(condition == i)
          for(r in paths_i$batch){
            if(r %in% paths_k$batch){
              paths_i_r = paths_i %>% dplyr::filter(batch == r)
              paths_k_r = paths_k %>% dplyr::filter(batch == r)
              
              hic_i = HiCExperiment::import(paths_i_r$iced_mcool[1], format = "mcool", resolution = j)
              hic_k = HiCExperiment::import(paths_k_r$iced_mcool[1], format = "mcool", resolution = j)
              
              for(c in 1:length(rownames(chrom_sizes))){
                chr = paste0(chrom_sizes$chr[c], ":1-", chrom_sizes$size[c])
                hic_i_c = subsetByOverlaps(hic_i, GRanges(chr), type = "within")
                hic_k_c = subsetByOverlaps(hic_k, GRanges(chr), type = "within")
                
                ratio = HiContacts::divide(hic_k_c, by = hic_i_c, use.scores = "balanced")
                cowplot::plot_grid(
                  HiContacts::plotMatrix(hic_k_c, compare.to = hic_i_c, use.scores = "balanced", scale = "log10"), 
                  HiContacts::plotMatrix(
                    ratio, 
                    use.scores = "balanced.fc", 
                    scale = "log2", 
                    cmap = HiContacts::bwrColors()
                  )
                )
                ggsave(paste0(ratio_path, "/", chrom_sizes$chr[c], "/ratio_", as.character(j), "_", k, "_vs_", i, "_", r, ".pdf"), width = 20, height = 10)
              }
            }
          }
        }
      }
    }
    
    #HiCrep
    paths = table %>% dplyr::select("sample_ID", "condition", "batch", paste0("iced_cool_", as.character(j)))
    colnames(paths) = c("sample_ID", "condition", "batch", "iced_cool")
    rownames(paths) = paths$sample_ID
    
    combinations = c()
    for(i in paths$sample_ID){
      for(k in paths$sample_ID){
        if(i != k){
          ordered = as.vector(sort(c(i, k)))
          n1 = ordered[1]
          n2 = ordered[2]
          combinations = append(combinations, paste0(n1, "_", n2))
        } else{combinations = append(combinations, paste0(i, "_", k))}
      }
    }
    combinations = unique(combinations)
    tab_all = data.frame(combs = combinations)
    
    tab = data.frame(rn = paths$sample_ID)
    cnames = colnames(tab)
    combs = c()
    names = c()
    for(i in 1:length(rownames(paths))){
      scc_i = c()
      sample_i = cool2matrix(paths$iced_cool[i])
      for(k in 1:length(rownames(paths))){
        sample_k = cool2matrix(paths$iced_cool[k])
        if(nrow(sample_k) > nrow(sample_i)) {
          score = get.scc(sample_i, sample_k, resol = as.integer(j), h = 5, lbr = 0, ubr = 5000000)
        } else {
          score = get.scc(sample_k, sample_i, resol = as.integer(j), h = 5, lbr = 0, ubr = 5000000)
        }
        scc = score[["scc"]]
        scc_i = append(scc_i, scc)
        
        if(paths$sample_ID[i] != paths$sample_ID[k]){
          ordered = as.vector(sort(c(paths$sample_ID[i], paths$sample_ID[k])))
          n1 = ordered[1]
          n2 = ordered[2]
          comb = paste0(n1, "_", n2)
        } else{comb = paste0(paths$sample_ID[i], "_", paths$sample_ID[k])}
        combs = append(combs, scc)
        names = append(names, comb)
      }
      tab$n = as.numeric(scc_i)
      cnames = append(cnames, paths$sample_ID[i])
      colnames(tab) = cnames
    }
    
    rownames(tab) = tab$rn
    tab = tab %>% dplyr::select(-rn)
    ann_row = as.data.frame(paths %>% dplyr::select(batch))
    rownames(ann_row) = paths$sample_ID
    ann_col = as.data.frame(paths %>% dplyr::select(condition))
    rownames(ann_col) = paths$sample_ID
    tab = as.data.frame(apply(tab, 2, as.numeric))
    rownames(tab) = colnames(tab)
    if(length(rownames(paths)) > 2){
      as.ggplot(pheatmap(tab, annotation_row = ann_row, annotation_col = ann_col, display_numbers = T, number_format = "%.2f", number_color = "black"))
      ggsave(filename = paste0(HiCrep_path, "/scc_score.pdf"), width = 10, height = 10)
    }
              
    for(c in chrom_sizes$chr){
      tab = data.frame(rn = paths$sample_ID)
      cnames = colnames(tab)
      combs = c()
      names = c()
      for(i in 1:length(rownames(paths))){
        scc_i = c()
        sample_i = cool2matrix(paths$iced_cool[i], chr = c)
        for(k in 1:length(rownames(paths))){
          sample_k = cool2matrix(paths$iced_cool[k], chr = c)
          if(nrow(sample_k) > nrow(sample_i)) {
            score = get.scc(sample_i, sample_k, resol = as.integer(j), h = 5, lbr = 0, ubr = 5000000)
          } else {
            score = get.scc(sample_k, sample_i, resol = as.integer(j), h = 5, lbr = 0, ubr = 5000000)
          }
          scc = score[["scc"]]
          scc_i = append(scc_i, scc)
          
          if(paths$sample_ID[i] != paths$sample_ID[k]){
            ordered = as.vector(sort(c(paths$sample_ID[i], paths$sample_ID[k])))
            n1 = ordered[1]
            n2 = ordered[2]
            comb = paste0(n1, "_", n2)
          } else{comb = paste0(paths$sample_ID[i], "_", paths$sample_ID[k])}
          combs = append(combs, scc)
          names = append(names, comb)
        }
        tab$n = as.numeric(scc_i)
        cnames = append(cnames, paths$sample_ID[i])
        colnames(tab) = cnames
      }
      tab_i = data.frame(comb = names, chr = combs) %>% dplyr::distinct()
      colnames(tab_i) = c("combs", c)
      tab_all = tab_all %>% dplyr::left_join(tab_i, by = "combs")
      
      rownames(tab) = tab$rn
      tab = tab %>% dplyr::select(-rn)
      ann_row = as.data.frame(paths %>% dplyr::select(batch))
      rownames(ann_row) = paths$sample_ID
      ann_col = as.data.frame(paths %>% dplyr::select(condition))
      rownames(ann_col) = paths$sample_ID
      tab = as.data.frame(apply(tab, 2, as.numeric))
      rownames(tab) = colnames(tab)
      if(length(rownames(paths)) > 2){
        as.ggplot(pheatmap(tab, annotation_row = ann_row, annotation_col = ann_col, display_numbers = T, number_format = "%.2f", number_color = "black"))
        ggsave(filename = paste0(HiCrep_path, "/scc_score_", c, ".pdf"), width = 10, height = 10)
      }
    }
    
    melt = melt(tab_all)
    ggplot(melt) + geom_point(aes(x = variable, y = value, color = combs))
    ggsave(filename = paste0(HiCrep_path, "/scc_score_all.pdf"), width = 20, height = 10)
    
    #dcHiC
    if("bias" %in% colnames(table) | paste0("bias_", as.character(j)) %in% colnames(table)){
      paths = table %>% dplyr::select(c("sample_ID", "condition", paste0("bed_", as.character(j)), paste0("iced_matrix_", as.character(j)), paste0("bias_", as.character(j))))
      colnames(paths) = c("sample_ID", "condition", "bed", "iced_matrix", "bias")
      new_beds = c()
      dir.create("biases")
      dir.create("beds")
      for(i in 1:length(rownames(paths))){
        bed_i = read.delim(paths$bed[i], header = FALSE)
        bias_i = read.delim(paths$bias[i], header = FALSE)
        n_bed = length(rownames(bed_i))
        n_bias = length(rownames(bias_i))
        diff = n_bed - n_bias
        V1 = c(rep("NaN", diff))
        extra_bins = as.data.frame(V1)
        bias_i = rbind(bias_i, extra_bins)
        bed_i$bias = bias_i$V1
        colnames(bed_i) = c("chr", "start", "end", "n", "bias")
        bed_i = bed_i %>% dplyr::filter(chr %in% chrom_sizes$chr) %>% dplyr::filter(! chr %in% c("chrY", "chrM"))
        
        new_bias = bed_i %>% dplyr::mutate(center = as.integer(start + as.integer(j) * 1/2)) %>% dplyr::select(c("chr", "center", "bias"))
        biases_path = paste0("biases/", paths$sample_ID[i], ".biases")
        write.table(new_bias, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, na = "nan", file = biases_path)
        system2("gzip", biases_path)
        
        new_bed = bed_i %>% dplyr::select(c("chr", "start", "end", "n"))
        beds_path = paste0("beds/", paths$sample_ID[i], ".bed")
        write.table(new_bed, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, na = "nan", file = beds_path)
        new_beds = append(new_beds, beds_path)
      }
    } else {
      paths = table %>% dplyr::select(c("sample_ID", "condition", paste0("bed_", as.character(j)), paste0("iced_matrix_", as.character(j))))
      colnames(paths) = c("sample_ID", "condition", "bed", "iced_matrix")
      new_beds = c()

      dir.create("beds")
      for(i in 1:length(rownames(paths))){
        bed_i = read.delim(paths$bed[i], header = FALSE)
        colnames(bed_i) = c("chr", "start", "end", "n")
        bed_i = bed_i %>% dplyr::filter(chr %in% chrom_sizes$chr) %>% dplyr::filter(! chr %in% c("chrY", "chrM"))
        
        new_bed = bed_i %>% dplyr::select(c("chr", "start", "end", "n"))
        beds_path = paste0("beds/", paths$sample_ID[i], ".bed")
        write.table(new_bed, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, na = "nan", file = beds_path)
        new_beds = append(new_beds, beds_path)
      }
    }
    
    
    paths$new_bed = new_beds
    
    dcHiC_matrix_all = paths %>% dplyr::select(c("iced_matrix", "new_bed", "sample_ID", "condition"))
    write.table(dcHiC_matrix_all, file = paste0(dcHiC_path, "/matrix_", as.character(j), "_all.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
    
    #calculate PCA 
    conda_run2(
      envname = "dchic",
      cmd_line = paste0("/usr/local/anaconda/envs/dchic/lib/R/bin/Rscript /dcHiC/dchicf.r --file ", dcHiC_path, "/matrix_", as.character(j), "_all.txt --pc ", PC, " --pcatype cis --dirovwt T"),
      echo = TRUE)
    
    #select the best PC
    conda_run2(
      envname = "dchic",
      cmd_line = paste0("/usr/local/anaconda/envs/dchic/lib/R/bin/Rscript /dcHiC/dchicf.r --file ", dcHiC_path, "/matrix_", as.character(j), "_all.txt --pcatype select --genome ", genome, " --dirovwt T"),
      echo = TRUE)
    
    #differential analysis for each reference
    for(k in references){
      for(i in paths$condition){
        if(i != k){
          dcHiC_matrix_i = dcHiC_matrix_all %>% dplyr::filter(condition %in% c(k, i))
          write.table(dcHiC_matrix_i, file = paste0(dcHiC_path, "/matrix_", as.character(j), "_", k, "_", i, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
          path_matrix = paste0(dcHiC_path, "/matrix_", as.character(j), "_", k, "_", i, ".txt")
          
          diffdir = paste0(k, "_vs_", i)
          
          #differential compartments
          conda_run2(
            envname = "dchic",
            cmd_line = paste0("/usr/local/anaconda/envs/dchic/lib/R/bin/Rscript /dcHiC/dchicf.r --file ", path_matrix, " --pcatype analyze --diffdir ", diffdir, " --dirovwt T"),
            echo = TRUE)
          
          #differential subcompartments
          conda_run2(
            envname = "dchic",
            cmd_line = paste0("/usr/local/anaconda/envs/dchic/lib/R/bin/Rscript /dcHiC/dchicf.r --file ", path_matrix, " --pcatype subcomp --diffdir ", diffdir, " --dirovwt T"),
            echo = TRUE)
          
          #enrichment in A compartment
          conda_run2(
            envname = "dchic",
            cmd_line = paste0("/usr/local/anaconda/envs/dchic/lib/R/bin/Rscript /dcHiC/dchicf.r --file ", path_matrix, " --pcatype enrich --diffdir ", diffdir, " --genome ", genome, " --dirovwt T --exclA T --region anchor"),
            echo = TRUE)
          
          if("bias" %in% colnames(table) | paste0("bias_", as.character(j)) %in% colnames(table)){
            #identify fithic loops for each sample
            conda_run2(
              envname = "dchic",
              cmd_line = paste0("/usr/local/anaconda/envs/dchic/lib/R/bin/Rscript /dcHiC/dchicf.r --file ", path_matrix, " --pcatype fithic --diffdir ", diffdir, " --dirovwt T --fithic /fithic/fithic/fithic.py --python python"),
              echo = TRUE)
            
            #differential loops
            conda_run2(
              envname = "dchic",
              cmd_line = paste0("/usr/local/anaconda/envs/dchic/lib/R/bin/Rscript /dcHiC/dchicf.r --file ", path_matrix, " --pcatype dloop --diffdir ", diffdir, " --dirovwt T"),
              echo = TRUE)
          }
          
          #create genome tracks
          conda_run2(
            envname = "dchic",
            cmd_line = paste0("/usr/local/anaconda/envs/dchic/lib/R/bin/Rscript /dcHiC/dchicf.r --file ", path_matrix, " --pcatype viz --diffdir ", diffdir, " --genome ", genome, " --pcgroup pcQnm"),
            echo = TRUE)
          
          #move into selected_path
          system2("mv", paste0("DifferentialResult/", diffdir, " /home/shared_folder/ooutput_HiC/untargeted_analysis/dcHiC/", as.character(j), "/differential_results_", k, "_vs_", i))
          system2("rm", "-r DifferentialResult")
          
          #move to IGV folder
          viz_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/dcHiC/", as.character(j), "/differential_results_", k, "_vs_", i, "/viz/vizIGV_intra/data/")
          
          if("bias" %in% colnames(table) | paste0("bias_", as.character(j)) %in% colnames(table)){
            if(file.exists(paste0(IGV_path, "diff_loops_", as.character(j), "_", k, "_vs_", i, ".bedpe")) == FALSE){
              system2("cp", paste0(viz_path, "differential.intra_compartmentLoops.bedpe ", IGV_path, "diff_loops_", as.character(j), "_", k, "_vs_", i, ".bedpe"))
            }
            
            diff_loops = read.delim(paste0(IGV_path, "diff_loops_", as.character(j), "_", k, "_vs_", i, ".bedpe"), header=FALSE)
            colnames(diff_loops) = c("chr1", "start1", "end1", "chr2", "start2", "end2", "pval")
            diff_loops = diff_loops %>% dplyr::mutate(loop = paste0(chr1, "_", start1, "_", chr2, "_", start2))
            
            FithicResult = read.delim(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/dcHiC/", as.character(j), "/differential_results_", k, "_vs_", i, "/fithic_run/FithicResult.txt"))
            FithicResult$loop = rownames(FithicResult)
            diff_loops = diff_loops %>% dplyr::left_join(FithicResult, by = "loop")
            
            if(file.exists(paste0(IGV_path, "diff_loops_",  as.character(j), "_", k, "_diff_", i, ".bedpe")) == FALSE){
              diff_loops_k = diff_loops %>% dplyr::select(c("chr1", "start1", "end1", "chr2", "start2", "end2", "pval", k))
              colnames(diff_loops_k) = c("chr1", "start1", "end1", "chr2", "start2", "end2", "pval", "k")
              diff_loops_k = diff_loops_k %>% dplyr::filter(k == 1) %>% dplyr::select(-k)
              write.table(diff_loops_k, file = paste0(IGV_path, "diff_loops_", as.character(j), "_", k, "_diff_", i, ".bedpe"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
            }
            
            if(file.exists(paste0(IGV_path, "diff_loops_", as.character(j), "_", i, "_diff_", k, ".bedpe")) == FALSE){
              diff_loops_i = diff_loops %>% dplyr::select(c("chr1", "start1", "end1", "chr2", "start2", "end2", "pval", i))
              colnames(diff_loops_i) = c("chr1", "start1", "end1", "chr2", "start2", "end2", "pval", "i")
              diff_loops_i = diff_loops_i %>% dplyr::filter(i == 1) %>% dplyr::select(-i)
              write.table(diff_loops_i, file = paste0(IGV_path, "diff_loops_", as.character(j), "_", i, "_diff_", k, ".bedpe"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
            }
            if(file.exists(paste0(IGV_path, "diff_loops_log10adjPval_", as.character(j), "_", k, "_vs_", i, ".bedpe")) == FALSE){
              diff_loops = diff_loops %>% dplyr::select(c("chr1", "start1", "end1", "chr2", "start2", "end2", "pval"))
              write.table(diff_loops, file = paste0(IGV_path, "diff_loops_log10adjPval_", as.character(j), "_", k, "_vs_", i, ".bedpe"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
            }
          }
            
          
          if(file.exists(paste0(IGV_path, "diff_compartments_distance_", as.character(j), "_", k, "_vs_", i, ".bedGraph")) == FALSE){
            system2("rm", paste0(viz_path, "differential_compartment.Mahalanobis.bedGraph"))
            system2("gzip", paste0("-dk ", viz_path, "differential_compartment.Mahalanobis.bedGraph.gz"))
            system2("mv", paste0(viz_path, "differential_compartment.Mahalanobis.bedGraph ", IGV_path,  "diff_compartments_distance_", as.character(j), "_", k, "_vs_", i, ".bedGraph"))
          }
          if(file.exists(paste0(IGV_path, "diff_compartments_log10adjPval_", as.character(j), "_", k, "_vs_", i, ".bedGraph")) == FALSE){
            system2("gzip", paste0("-dk ", viz_path, "differential_compartment.log10Padj.bedGraph.gz"))
            system2("mv", paste0(viz_path, "differential_compartment.log10Padj.bedGraph ", IGV_path,  "diff_compartments_log10adjPval_", as.character(j), "_", k, "_vs_", i, ".bedGraph"))
          }
          if(file.exists(paste0(IGV_path, "compartments_", i, "_", as.character(j), ".bedGraph")) == FALSE){
            system2("gzip", paste0("-dk ", viz_path, i, ".PC.bedGraph.gz"))
            system2("mv", paste0(viz_path, i, ".PC.bedGraph ", IGV_path,  "compartments_", i, "_", as.character(j), ".bedGraph"))
          }
          if(file.exists(paste0(IGV_path, "subcompartments_", i, "_", as.character(j), ".bedGraph")) == FALSE){
            system2("gzip", paste0("-dk ", viz_path, "intra_", i, ".subcomp.seg.gz"))
            system2("mv", paste0(viz_path,  "intra_", i, ".subcomp.seg ", IGV_path,  "subcompartments_", i, "_", as.character(j), ".bedGraph"))
            subcomp = read.delim(paste0(IGV_path,  "subcompartments_", i, "_", as.character(j), ".bedGraph")) %>% dplyr::select(-c("Sample", "Num_Probes"))
            system2("rm", paste0(IGV_path, "subcompartments_", i, "_", as.character(j), ".bedGraph"))
            write.table(subcomp, file = paste0(IGV_path, "subcompartments_", i, "_", as.character(j), ".bedGraph"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
          }
          
          if(file.exists(paste0(IGV_path, "compartments_", k, "_", as.character(j), ".bedGraph")) == FALSE){
            system2("gzip", paste0("-dk ", viz_path, k, ".PC.bedGraph.gz"))
            system2("mv", paste0(viz_path, k, ".PC.bedGraph ", IGV_path,  "compartments_", k, "_", as.character(j), ".bedGraph"))
          }
          if(file.exists(paste0(IGV_path, "subcompartments_", k, "_", as.character(j), ".bedGraph")) == FALSE){
            system2("gzip", paste0("-dk ", viz_path, "intra_", k, ".subcomp.seg.gz"))
            system2("mv", paste0(viz_path,  "intra_", k, ".subcomp.seg ", IGV_path,  "subcompartments_", k, "_", as.character(j), ".bedGraph"))
            subcomp = read.delim(paste0(IGV_path,  "subcompartments_", k, "_", as.character(j), ".bedGraph")) %>% dplyr::select(-c("Sample", "Num_Probes"))
            system2("rm", paste0(IGV_path, "subcompartments_", k, "_", as.character(j), ".bedGraph"))
            write.table(subcomp, file = paste0(IGV_path, "subcompartments_", k, "_", as.character(j), ".bedGraph"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
          }
          if(file.exists(paste0(IGV_path, "diff_subcompartments_", k, "_vs_", i, "_", as.character(j), ".bedGraph")) == FALSE){
            subcomp = read.delim(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/dcHiC/", as.character(j), "/differential_results_", k, "_vs_", i, "/fdr_result/intra_sample_group.subcompartments.bedGraph"))
            subcomp = subcomp %>% dplyr::select(c("chr", "start", "end", "padj")) %>% dplyr::mutate(padj = -log10(padj))
            write.table(subcomp, file = paste0(IGV_path, "diff_subcompartments_log10adjPval_", as.character(j), "_", k, "_vs_", i, ".bedGraph"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
          }
          
          samples = (table %>% dplyr::filter(condition %in% c(i, k)))$sample_ID
          viz_path = paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/dcHiC/", as.character(j), "/differential_results_", k, "_vs_", i, "/viz/files/")
          
          for(m in samples){
            if(file.exists(paste0(IGV_path, "compartments_", m, "_", as.character(j), ".bedGraph")) == FALSE){
              system2("cp", paste0(viz_path,  "intra_", m, "_PC.bedGraph ", IGV_path, "compartments_", m, "_", as.character(j), ".bedGraph"))
            }
          }
          
          #gsea
          gsea_file_i = read.delim(paste0(dcHiC_path, "/differential_results_", k, "_vs_", i, "/geneEnrichment/", i, "_geneEnrichment/", i, "_geneEnrichment.anchor.txt"))
          
          if("GeneOntologyMolecularFunction" %in% gsea_file_i$Category){
            mf = gsea_file_i %>% dplyr::filter(Category == "GeneOntologyMolecularFunction")
            ggplot(mf, aes(x = reorder(Name, GenesInTermInQuery), y = GenesInTermInQuery, fill = -log10(QValueFDRBH))) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
              labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO MF genes in A in ", i, " vs ", k))
            ggsave(paste0(gsea_path, "/gsea_", as.character(j), "_MF_", i, "_vs_", k, ".pdf"), width = 15, height = 5)
            
            write.csv(mf, paste0(gsea_path, "/gsea_", as.character(j), "_MF_", i, "_vs_", k, ".csv"))
          }
          
          if("GeneOntologyCellularComponent" %in% gsea_file_i$Category){
            cc = gsea_file_i %>% dplyr::filter(Category == "GeneOntologyCellularComponent")
            ggplot(cc, aes(x = reorder(Name, GenesInTermInQuery), y = GenesInTermInQuery, fill = -log10(QValueFDRBH))) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
              labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO CC genes in A in ", i, " vs ", k))
            ggsave(paste0(gsea_path, "/gsea_", as.character(j), "_CC_", i, "_vs_", k, ".pdf"), width = 15, height = 5)
            
            write.csv(cc, paste0(gsea_path, "/gsea_", as.character(j), "_CC_", i, "_vs_", k, ".csv"))
          }
          
          if("GeneOntologyBiologicalProcess" %in% gsea_file_i$Category){
            bp = gsea_file_i %>% dplyr::filter(Category == "GeneOntologyBiologicalProcess")
            ggplot(bp, aes(x = reorder(Name, GenesInTermInQuery), y = GenesInTermInQuery, fill = -log10(QValueFDRBH))) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
              labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO BP genes in A in ", i, " vs ", k))
            ggsave(paste0(gsea_path, "/gsea_", as.character(j), "_BP_", i, "_vs_", k, ".pdf"), width = 15, height = 5)
            
            write.csv(bp, paste0(gsea_path, "/gsea_", as.character(j), "_BP_", i, "_vs_", k, ".csv"))
          }
          
          
          gsea_file_k = read.delim(paste0(dcHiC_path, "/differential_results_", k, "_vs_", i, "/geneEnrichment/", k, "_geneEnrichment/", k, "_geneEnrichment.anchor.txt"))
          
          if("GeneOntologyMolecularFunction" %in% gsea_file_k$Category){
            mf = gsea_file_k %>% dplyr::filter(Category == "GeneOntologyMolecularFunction")
            ggplot(mf, aes(x = reorder(Name, GenesInTermInQuery), y = GenesInTermInQuery, fill = -log10(QValueFDRBH))) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
              labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO MF genes in A in ", k, " vs ", i))
            ggsave(paste0(gsea_path, "/gsea_", as.character(j), "_MF_", k, "_vs_", i, ".pdf"), width = 15, height = 5)
            
            write.csv(mf, paste0(gsea_path, "/gsea_", as.character(j), "_MF_", k, "_vs_", i, ".csv"))
          }
          
          if("GeneOntologyCellularComponent" %in% gsea_file_k$Category){
            cc = gsea_file_k %>% dplyr::filter(Category == "GeneOntologyCellularComponent")
            ggplot(cc, aes(x = reorder(Name, GenesInTermInQuery), y = GenesInTermInQuery, fill = -log10(QValueFDRBH))) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
              labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO CC genes in A in ", k, " vs ", i))
            ggsave(paste0(gsea_path, "/gsea_", as.character(j), "_CC_", k, "_vs_", i, ".pdf"), width = 15, height = 5)
            
            write.csv(cc, paste0(gsea_path, "/gsea_", as.character(j), "_CC_", k, "_vs_", i, ".csv"))
          }
          
          if("GeneOntologyBiologicalProcess" %in% gsea_file_k$Category){
            bp = gsea_file_k %>% dplyr::filter(Category == "GeneOntologyBiologicalProcess")
            ggplot(bp, aes(x = reorder(Name, GenesInTermInQuery), y = GenesInTermInQuery, fill = -log10(QValueFDRBH))) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
              labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO BP genes in A in ", k, " vs ", i))
            ggsave(paste0(gsea_path, "/gsea_", as.character(j), "_BP_", k, "_vs_", i, ".pdf"), width = 15, height = 5)
            
            write.csv(bp, paste0(gsea_path, "/gsea_", as.character(j), "_BP_", k, "_vs_", i, ".csv"))
          }

          }
        }
      }
    
    #PCA
    if(nrow(table) > 2){
      dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/comp_correlation"))
      PCs = c()
      for(i in 1:length(rownames(table))){
        PC_i = paste0(IGV_path,  "compartments_", table$sample_ID[i], "_", as.character(j), ".bedGraph")
        PCs = append(PCs, PC_i)
      }
      
      table$PC = PCs
      PC_table = read.csv(table$PC[1], header=FALSE, skip = 4, sep = "\t")
      colnames(PC_table) = c("chr", "start", "end", table$sample_ID[1])
      for(i in 2:length(rownames(table))){
        PC_i = read.csv(table$PC[i], header=FALSE, skip = 4, sep = "\t")
        colnames(PC_i) = c("chr", "start", "end", table$sample_ID[i])
        PC_table = PC_table %>% dplyr::left_join(PC_i, by = c("chr", "start", "end"))
      }
      
      PC_matrix = PC_table %>% dplyr::select(-c("chr", "start", "end"))
      PC_matrix <- na.omit(PC_matrix)
      n = length(rownames(table))
      c_names = colnames(PC_matrix) 
      colnames(PC_matrix) = c(1:n)
      PC_matrix = PC_matrix %>% dplyr::mutate(sum = rowSums(PC_matrix)) %>% dplyr::filter(sum != 0) %>% dplyr::filter(sum != n * PC_matrix$'1') %>% dplyr::select(- sum)
      colnames(PC_matrix) = c_names
      PC_matrix_t = t(PC_matrix)
      data.pca <- prcomp(PC_matrix_t, scale = TRUE)
      fviz_pca_ind(data.pca)
      ggsave(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/comp_correlation/PCA_", as.character(j), ".pdf"), width = 12, height = 8)
      
      corr = rcorr(as.matrix(PC_matrix))
      pearson = as.data.frame(corr[["r"]])
      as.ggplot(pheatmap(pearson))
      ggsave(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/comp_correlation/Pearson_", as.character(j), ".pdf"), width = 12, height = 8)
    }
    
    #genomeTrack
    paths = table %>% dplyr::select(c("sample_ID", paste0("iced_h5_", as.character(j))))
    colnames(paths) = c("sample_ID", "iced_h5")
    
    for(c in 1:length(rownames(chrom_sizes))){
      dir.create(paste0("/home/shared_folder/ooutput_HiC/untargeted_analysis/genomeTracks/", as.character(j), "/", chrom_sizes$chr[c]))
    }
    
    for(i in 1:length(rownames(paths))){
      system2("cp", paste0("/individual_samples_track ", track_path, "track_", paths$sample_ID[i]))
      system2("sed", paste0("-i 's|file = file_matrix|file = ", paths$iced_h5[i], "|g' ", track_path, "track_", paths$sample_ID[i]))
      system2("sed", paste0("-i 's|file = file_tads|file = ", TADs_path, "/", paths$sample_ID[i], "_", as.character(j), "_domains.bed|g' ", track_path, "track_", paths$sample_ID[i]))
      system2("sed", paste0("-i 's|file = file_pc|file = ", IGV_path, "compartments_", paths$sample_ID[i], "_", as.character(j), ".bedGraph|g' ", track_path, "track_", paths$sample_ID[i]))
      system2("sed", paste0("-i 's|title = HiC matrix|title = HiC matrix ", paths$sample_ID[i], "|g' ", track_path, "track_", paths$sample_ID[i]))
      system2("sed", paste0("-i 's|title = compartments|title = compartments ", paths$sample_ID[i], "|g' ", track_path, "track_", paths$sample_ID[i]))
      
      for(c in 1:length(rownames(chrom_sizes))){
        system2("sed", paste0("-i 's|depth =|depth = ", chrom_sizes$size[c], "|g' ", track_path, "track_", paths$sample_ID[i]))
        system2("pyGenomeTracks", paste0("--tracks ", track_path, "track_", paths$sample_ID[i], " --region ", chrom_sizes$chr[c], ":1-", chrom_sizes$size[c], " --outFileName ", track_path, chrom_sizes$chr[c], "/track_", paths$sample_ID[i], "_", chrom_sizes$chr[c], ".pdf"))
        system2("sed", paste0("-i 's|depth = ", chrom_sizes$size[c], "|depth =|g' ", track_path, "track_", paths$sample_ID[i]))
      }
    }
    
    for(k in references){
      for(i in unique(table$condition)){
        if(i != k){
          system2("cp", paste0("/comparison_track ", track_path, "track_", k, "_vs_", i))
          system2("sed", paste0("-i 's|file = file_pc_k|file = ", IGV_path, "compartments_", k, "_", as.character(j), ".bedGraph|g' ", track_path, "track_", k, "_vs_", i))
          system2("sed", paste0("-i 's|file = file_pc_i|file = ", IGV_path, "compartments_", i, "_", as.character(j), ".bedGraph|g' ", track_path, "track_", k, "_vs_", i))
          
          system2("sed", paste0("-i 's|file = file_sub_k|file = ", IGV_path, "subcompartments_", k, "_", as.character(j), ".bedGraph|g' ", track_path, "track_", k, "_vs_", i))
          system2("sed", paste0("-i 's|file = file_sub_i|file = ", IGV_path, "subcompartments_", i, "_", as.character(j), ".bedGraph|g' ", track_path, "track_", k, "_vs_", i))
          
          system2("sed", paste0("-i 's|file = distance_file|file = ", IGV_path, "diff_compartments_distance_", as.character(j), "_", k, "_vs_", i, ".bedGraph|g' ", track_path, "track_", k, "_vs_", i))
          system2("sed", paste0("-i 's|file = pval_file_comp|file = ", IGV_path, "diff_compartments_log10adjPval_", as.character(j), "_", k, "_vs_", i, ".bedGraph|g' ", track_path, "track_", k, "_vs_", i))
          system2("sed", paste0("-i 's|file = pval_file_subcomp|file = ", IGV_path, "diff_subcompartments_log10adjPval_", as.character(j), "_", k, "_vs_", i, ".bedGraph|g' ", track_path, "track_", k, "_vs_", i))
          
          system2("sed", paste0("-i 's|title = compartments_k|title = compartments ", k, "|g' ", track_path, "track_", k, "_vs_", i))
          system2("sed", paste0("-i 's|title = compartments_i|title = compartments ", i, "|g' ", track_path, "track_", k, "_vs_", i))
          
          system2("sed", paste0("-i 's|title = subcomp_k|title = subcompartments ", k, "|g' ", track_path, "track_", k, "_vs_", i))
          system2("sed", paste0("-i 's|title = subcomp_i|title = subcompartments ", i, "|g' ", track_path, "track_", k, "_vs_", i))
          
          if("bias" %in% colnames(table) | paste0("bias_", as.character(j)) %in% colnames(table)){
            system2("sed", paste0("-i 's|file = pval_file_loops|file = ", IGV_path, "diff_loops_log10adjPval_", as.character(j), "_", k, "_vs_", i, ".bedpe|g' ", track_path, "track_", k, "_vs_", i))
            system2("sed", paste0("-i 's|file = loops_file_k|file = ", IGV_path, "diff_loops_", as.character(j), "_", k, "_diff_", i, ".bedpe|g' ", track_path, "track_", k, "_vs_", i))
            system2("sed", paste0("-i 's|file = loops_file_i|file = ", IGV_path, "diff_loops_", as.character(j), "_", i, "_diff_", k, ".bedpe|g' ", track_path, "track_", k, "_vs_", i))
            system2("sed", paste0("-i 's|title = loops k|title = ", k, "|g' ", track_path, "track_", k, "_vs_", i))
            system2("sed", paste0("-i 's|title = loops i|title = ", i, "|g' ", track_path, "track_", k, "_vs_", i))
          } else {
            system2("sed", paste0("-i '76, $ d' ", track_path, "track_", k, "_vs_", i))
          }
          
          
          for(c in 1:length(rownames(chrom_sizes))){
            system2("pyGenomeTracks", paste0("--tracks ", track_path, "track_", k, "_vs_", i, " --region ", chrom_sizes$chr[c], ":1-", chrom_sizes$size[c], " --outFileName ", track_path, chrom_sizes$chr[c], "/track_", k, "_vs_", i, "_", chrom_sizes$chr[c], ".pdf"))
          }
        }
      }
    }
  }
}

#ANALYSIS - targeted
if(targeted_analysis == TRUE){
  dir.create(paste0("/home/shared_folder/ooutput_HiC/targeted_analysis"))
  
  #virtual 4C
  if(virtual4C == TRUE){
    dir.create(paste0("/home/shared_folder/ooutput_HiC/targeted_analysis/virtual_4C"))
    
    for(v in 1:length(rownames(vp))){
      v_bed = data.frame(chr = vp$chr[v], start = vp$start[v], end = vp$end[v])
      viewpoint = paste0(vp$chr[v], ":", vp$start[v], ":", vp$end[v])
      viewpoint_f = paste0(vp$chr[v], "_", vp$start[v], "_", vp$end[v])
      write.table(v_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, file = paste0("/home/shared_folder/ooutput_HiC/targeted_analysis/virtual_4C/", viewpoint_f, ".bed"))
    }
    print(colnames(table))
    for(j in res_v4C){
      dir.create(paste0("/home/shared_folder/ooutput_HiC/targeted_analysis/virtual_4C/", as.character(j)))
      vir4C_path = paste0("/home/shared_folder/ooutput_HiC/targeted_analysis/virtual_4C/", as.character(j), "/")
      
      paths = table %>% dplyr::select("sample_ID", paste0("iced_mcool_", as.character(j)))
      colnames(paths) = c("sample_ID", "iced_mcool")
      
      for(i in 1:length(rownames(paths))){
        print(paths$iced_mcool[i])
        hic = HiCExperiment::import(paths$iced_mcool[i], format = "mcool", resolution = j)
        for(v in 1:length(rownames(vp))){
          viewpoint = paste0(vp$chr[v], ":", vp$start[v], "-", vp$end[v])
          viewpoint_f = paste0(vp$chr[v], "_", vp$start[v], "_", vp$end[v])
          v4C = HiContacts::virtual4C(hic, viewpoint = GRanges(viewpoint))
          write.table(v4C %>% as.data.frame() %>% dplyr::mutate(start = as.character(as.integer(start))) %>%
                        dplyr::mutate(end = as.character(as.integer(end))) %>% dplyr::select(c("seqnames", "start", "end", "score")),
                      sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, file = paste0(vir4C_path, paths$sample_ID[i], "_", viewpoint_f, ".bedGraph"))
          system2("cp", paste0("/4C_track ", vir4C_path, "track_", paths$sample_ID[i], "_", viewpoint_f))
          system2("sed", paste0("-i 's|file = file_4C|file = ", vir4C_path, paths$sample_ID[i], "_", viewpoint_f, ".bedGraph|g' ", vir4C_path, "track_", paths$sample_ID[i], "_", viewpoint_f))
          system2("sed", paste0("-i 's|file = file_viewpoint|file = /home/shared_folder/ooutput_HiC/targeted_analysis/virtual_4C/", viewpoint_f, ".bed|g' ", vir4C_path, "track_", paths$sample_ID[i], "_", viewpoint_f))
          
          system2("sed", paste0("-i 's|title = title_4C_auto|title = virtual 4C ", paths$sample_ID[i], "|g' ", vir4C_path, "track_", paths$sample_ID[i], "_", viewpoint_f))
          system2("sed", paste0("-i 's|title = title_4C_zoom|title = virtual 4C ", paths$sample_ID[i], " zoom|g' ", vir4C_path, "track_", paths$sample_ID[i], "_", viewpoint_f))
          
          start_range = as.integer(vp$start[v]) - as.integer(j)*10
          if(start_range < 0){start_range = 0}
          end_range = as.integer(vp$end[v]) + as.integer(j)*10
          chr_v = chrom_sizes %>% dplyr::filter(chr == vp$chr[v])
          if(end_range > as.integer(chr_v$size[1])){end_range = as.integer(chr_v$size[1])}
          range = paste0(vp$chr[v], ":", as.character(start_range), "-", as.character(end_range))
          system2("pyGenomeTracks", paste0("--tracks ", vir4C_path, "track_", paths$sample_ID[i], "_", viewpoint_f, " --region ", range, " --outFileName ", vir4C_path, paths$sample_ID[i], "_", viewpoint_f, ".pdf"))
        }
      }
    }
  }
}
