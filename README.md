# TTN_O_MAP_data_analysis
![Docker Image Version](https://img.shields.io/docker/v/sara-bianchi/TTN_O_MAP_data_analysis?style=flat-square&logo=docker&logoColor=white&color=2496ED)
![Docker Pulls](https://img.shields.io/docker/pulls/sara-bianchi/TTN_O_MAP_data_analysis?style=flat-square&logo=docker&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)

scripts for reproducibility of the bioinforatic analysis of "Nascent transcript O-MAP reveals the molecular architecture of an RBM20-regulated single-locus subnuclear compartment organized around the TTN RNA", Kania, Bianchi et al.

Our analyses include:

1) the alignment and preprocessing of raw data:
- O-MAP data processing: RNA ChIP pulldown and input data for WT and ΔRBM20 hESCs RUES2 derived CMs, where raw reads are available on GEO with the accession number GSE234808
- bulk data processing: bulk RNA seq of hiPSCs WTC11 derived CMs (WT, ΔRBM20, ΔTTNpromoter, ΔTTNpromoter-ΔRBM20), RUES2 KD CMs (Scramble, SAFB KD), and hESCs RUES2 derived CMs (WT, ΔRBM20, ΔTTNpromoter, ΔTTNprotein), where raw data are available on ArrayExprress with the following accession numbers, respecively E-MTAB-17833, E-MTAB-17834, and E-MTAB-17835
- HiC data processing: HiC data processing of RUES2 hESCs including the developmental trajectory of cardiac differentiation (ESCs day 0, MES day 2, CP day 5 and CMs day 14) available on 4DN Data Portal (Experiment Set codes 4DNESGTHHJAC, 4DNES6G787KF, 4DNESSU7RRNR, 4DNESW72J62T), and WT and ΔRBM20 day 30 CMs available on ENA at the following accession number PRJEB104212

2) further analyses and integration of the previous datasets with other published data starting from intermediate preprocessed files obtained in the previous steps and available on Zenodo at the following DOI 10.5281/zenodo.18722155. The repository also contains some settings and envrionment files for the processing of raw data.
It is highly recommanded to run all our scripts using our docker images available on the hedgelab docker repository to ensure reproducibility.

### O-MAP data processing
For this step it is required to use the scripts in the folder O-MAP_data_prcessing and the supporing files on Zenodo at https://zenodo.org/api/records/18722155/draft/files/O-MAP_data_processing.gz/content.
Place the files and the scripts in the respective folders (RNA_ChIP_data_processing and DNA_ChIP_data_processing), make the scripts executable, and download the fastq files for the RNA ChIP pulldown and input in a subfolder RNA_ChIP_data_processing/fastq

#### RNA ChIP data processing
RNA ChIP data processing consits of different steps requiring different docker images (specified below). All the steps can be run on the RNA_ChIP_data_processing folder as following:

1) STAR geome index generation and alignment of input and pulldown data
   
       docker pull hedgelab/bulk_image:image3
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_STAR hedgelab/bulk_image:image3
       docker exec -it RNA_ChIP_STAR /bin/bash
       /home/shared_folder/STAR_alignment.sh
   
2) Conversion of BAM files to .bed tracks on the TTN genes, normalization on the coverage and generation of .bedGraph tracks with the log2FC pulldown/input

       docker pull hedgelab/o_map_data_analysis:image1
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_introns hedgelab/o_map_data_analysis:image1
       docker exec -it RNA_ChIP_introns /bin/bash
       /home/shared_folder/bamToBed_input.sh
       /home/shared_folder/bamToBed_pulldown.sh
       Rscript script_introns_alignment.R
   
4) rMATS analysis of differential splicing events from STAR BAM files of pulldown data and sashimi plots generation

       docker pull hedgelab/rmast:image1
       docker run -itv ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name RNA_ChIP_rMATS hedgelab/rmast:image1 /bin/bash
       /home/shared_folder/splicing_variants_analysis.sh
   
       docker pull hedgelab/hedgelab/sashimiplot:image1
       docker run -d -v $PWD:/home/shared_folder --name RNA_ChIP_sashimi hedgelab/sashimiplot:image1
       docker exec -it RNA_ChIP_sashimi /bin/bash
       /home/shared_folder/script_sashimi.sh
   
6) BWA alignment of input and pulldown data

       docker pull hedgelab/ciriquant:image1
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_CIRIquant hedgelab/ciriquant:image1
       docker exec -it RNA_ChIP_CIRIquant /bin/bash
       /home/shared_folder/BWA_alignment.sh
   
7) CIRI2 detection of circRNAs

       docker pull hedgelab/ciri_full:image3
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_CIRI2 hedgelab/ciri_full:image3
       docker exec -it RNA_ChIP_CIRI2 /bin/bash
       /home/shared_folder/CIRI2.sh
   
8) Quantification of differentially expressed circRNAs with CIRIquant

       docker exec -it RNA_ChIP_CIRIquant /bin/bash
       /home/shared_folder/CIRIquant.sh
   
9) Bowtie index generation, O-MAP probes alignment, and generation of .bed tracks for differentially enriched TTN introns and circRNAs or WT and ΔRBM20

       docker pull hedgelab/crispr:image1
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_tracks hedgelab/crispr:image1
       docker exec -it RNA_ChIP_tracks /bin/bash
       Rscript /home/shared_folder/circRNAs_overlaps.R
    

#### DNA ChIP data processing
DNA ChIP data processing consists in the download of .bw track files from GEO and conversion to the .bedGraph format for circular plots generation and in the binning of chromosome 2 tracks at 100kb to enable comparison with other genomic tracks with Pearson correlation analysis.
The analysis can be run on our docker hedgelab/? and it canbe performed by running the following steps in the terminal from the DNA_ChIP_data_processing folder:

     docker pull hedgelab/o_map_data_analysis:image1
     docker run -d -v ../DNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=DNA_ChIP hedgelab/o_map_data_analysis:image1
     docker exec -it DNA_ChIP /bin/bash
     /home/shared_folder/chromatin_O-MAP_analysis.sh


## Bulk data processing
For this step it is required to use the scripts in the folder bulk_data_prcessing and the supporing files on Zenodo at https://zenodo.org/api/records/18722155/draft/files/bulk_data_processing.gz/content.
Place the files and the scripts in the same folder, make the scripts executable, and download the fastq files for the RUES2 CMs, WTC11 CMs and RUES2 KD CMs in a subfolder bulk_data_processing/fastq

Bulk data processig consits of different steps requiring different docker images (specified below). All the steps can be run on the bulk_data_processing folder as following:
1) STAR geome index generation, trimming, alignment, counts filtering, TPM normalization, PCA, WGCNA modules analysis, differential gene expression, and GSEA analyses
   
       docker pull hedgelab/bulk_image:image3
       docker run -d -v ../bulk_data_processing:/home/shared_folder --privileged=true --name=bulk_DGE hedgelab/bulk_image:image3
       docker exec -it bulk_DGE /bin/bash
       /home/shared_folder/DGE_data_analysis.sh
   
2) For RUES2 CMs onversion of BAM files to .bed tracks on the TTN and RBM20 genes

       docker pull hedgelab/o_map_data_analysis:image1
       docker run -d -v ../bulk_data_processing:/home/shared_folder --privileged=true --name=RUES2_tracks hedgelab/o_map_data_analysis:image1
       docker exec -it RUES2_tracks /bin/bash
       /home/shared_folder/bamToBed_tracks.sh
   
4) rMATS analysis of differential splicing events from STAR BAM files of pulldown data and sashimi plots generation

       docker pull hedgelab/rmast:image1
       docker run -itv ../bulk_data_processing:/home/shared_folder --privileged=true --name bulk_data_rMATS hedgelab/rmast:image1 /bin/bash
       /home/shared_folder/splicing_variants_analysis.sh
   
       docker pull hedgelab/hedgelab/sashimiplot:image1
       docker run -d -v $PWD:/home/shared_folder --name bulk_data_sashimi hedgelab/sashimiplot:image1
       docker exec -it bulk_data_sashimi /bin/bash
       /home/shared_folder/sashimi_plots_RUES2.sh
       /home/shared_folder/sashimi_plots_WTC11.sh


## HiC data processing
For this step it is required to use the scripts in the folder HiC_data_prcessing and the supporing files on Zenodo at https://zenodo.org/api/records/18722155/draft/files/HiC_data_processing.gz/content.
Place the files and the scripts in the same folder, make the scripts executable, and download the fastq files for the RUES2 CMs (developmental tracks and WT-ΔRBM20 day 30 CMs) in a subfolder HiC_data_processing/fastq

The analysis can be run on our docker hedgelab/hic_pipleine:image15 and it canbe performed by running the following steps in the terminal from the HiC_data_processing folder:
1) HiC data processing consists in the Bowtie2 genome index generation, reads alignment with HiC-Pro, replicates similarity analysis with HiCRep, differential compartments and loops analysis with dcHiC and FitHiC, virtual 4C analysis from the TTN promoter with Bioconductor (HiContacts), TADs detection (HiCExplorer), and IGV tracks generation

       docker pull hedgelab/hic_pipeline:image15
       docker run -d -v ../HiC_data_processing:/home/shared_folder --privileged=true --name=HiC_processing hedgelab/hic_pipeline:image15
       docker exec -it HiC_processing /bin/bash
       Rscript /home/shared_folder/Rscript_HiC.R
   
2) Trasn-C analysis from the TTN seed and random background generation

       docker exec -it HiC_processing /bin/bash
       ./ome/shared_folder/script_transC.sh


## Processed data analysis
For this step it is required to use the scripts in the folder processed_data_analysis and the preprocessed files on Zenodo at https://zenodo.org/api/records/18722155/draft/files/processed_data_analysis.gz/content.
Place the files and the scripts in the same folder and make the scripts executable.

Processed data anlysis consits of different steps requiring different docker images (specified below). All the steps can be run on the bulk_data_processing folder as following:
1) pyGenomeTracks generation for visualization of compartments, virtual 4C tracks, O-MAP pulldown CIDs and TIDs, diferentially enriched introns and circRNAs, reads distribution on TTN and RBM20 in RUES2 CMs bulk data:

       docker pull hedgelab/o_map_data_analysis:image1
       docker run -d -v ../processed_data_analysis:/home/shared_folder --privileged=true --name=processed_data_analysis hedgelab/o_map_data_analysis:image1
       docker exec -it processed_data_analysis /bin/bash
       /home/shared_folder/genome_tracks_generation.sh

2) TIDs and CIDs analysis: cumulative distribution of CIDs and TIDs DNA FISH distances, Pearson correlation of genomic tracks on chromosome 2, HiCRep correlation, generation of circular plots for O-MAP TIDs, AUC and p-value calculation of trans-C results over the random background

       docker exec -it processed_data_analysis /bin/bash
       Rscript /home/shared_folder/chromatin_O-MAP_analysis.R
   
3) Validation of orthogonal perturbations of the RBM20 splicing factory in RUES2 CMs from gene expression and splicing data

       docker exec -it processed_data_analysis /bin/bash
       Rscript /home/shared_folder/RBM20KO_TTNpKO_validation.R
   
4) TIDs analysis of compartment changes in WT and ΔRBM20 CMs compared with gene expression of orthognal perturbations of the factory

        docker exec -it processed_data_analysis /bin/bash
        Rscript /home/shared_folder/TIDs_gene_expression_compartments_analysis.R
   
5) Analysis of O-MAP enriched transcripts for the overlap with genes on chromosome 2, RBM20 eCLIP targets and TIDs enclosed genes
   
        docker exec -it processed_data_analysis /bin/bash
        Rscript /home/shared_folder/pulldown_transcripts_analysis.R
   
6) Analyis of splicing variants in pulldown enriched transcripts and analysis of the rescue of the ΔRBM20 misspliced transcripts in the ΔRBM20-ΔTTNpromoter. Analysis of the enrichment of RBP motifs in rescued transcripts and comparison with proteomic pulldown.
   
        docker exec -it processed_data_analysis /bin/bash
        Rscript /home/shared_folder/splicing_analysis.R

7) Enrichment of GO terms in pulldown enriched proteins. Analysis of the number of QKI and SAFB foci per cells in WT, ΔRBM20, ΔTTNpromoter, ΔRBM20-ΔTTNpromoter hESC RUES2 derived CMs. Evaluation of the transcriptional effects of the QKI, SAFB, and ELAVL1 KD compared to the factory perturbations and evaluation of gene expression in TIDs enclosed genes.

        docker exec -it processed_data_analysis /bin/bash
        Rscript /home/shared_folder/proteins_analysis.R

   
For further details on the Bulk RNA seq and HiC processing pipelines refer respectively to Refer to our GitHubs:
- bulk RNAseq pipeline (https://github.com/sara-bianchi/Bulk_RNA_seq_pipeline)
- Hi-C pipeline (https://github.com/sara-bianchi/HiC_pipeline)
