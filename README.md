# TTN_O_MAP_data_analysis
scripts for reproducibility of the bioinforatic analysis of "Nascent trascript O-MAP revels the molecular architecture of a single-locus subnuclear compartment organized around RBM20 and the TTN RNA", Kania, Bianchi et al.

Our analyses include:

1) the alignment and preprocessing of raw data:
- O-MAP data processing: RNA ChIP pulldown and input data for WT and ΔRBM20 hESCs RUES2 derived CMs, where raw reads are available on GEO with the accession number GSE234808
- bulk data processing: bulk RNA seq of hESCs RUES2 derived CMs (WT, ΔRBM20, ΔTTNpromoter, ΔTTNprotein), hiPSCs WTC11 derived CMs (WT, ΔRBM20, ΔTTNpromoter, ΔTTNpromoter-ΔRBM20), RUES2 KD CMs (Scramble, SAFB KD, ELAVL1 KD, QKI KD), where raw data are available on ArrayExprress with the following accession numbers, respecively ? ? ?
- HiC data processing: HiC data processing of RUES2 hESCs including the developmental trajectory of cardiac differentiation (ESCs day 0, MES day 2, CP day 5 and CMs day 14) available on 4DN Data Portal (Experiment Set codes 4DNESGTHHJAC, 4DNES6G787KF, 4DNESSU7RRNR, 4DNESW72J62T), and WT and ΔRBM20 day 30 CMs available on ENA at the following accession number ?

2) further analyses and integration of the previous datasets with other published data starting from intermediate preprocessed files obtained in the previous steps and available on Zenodo at the following DOI 10.5281/zenodo.18722155. The repository also contains some settings and envrionment files for the processing of raw data.
It is highly recommanded to run all our scripts using our docker images available on the hedgelab docker repository to ensure reproducibility.

### O-MAP data processing
For this step it is required to use the scripts in the folder O-MAP_data_prcessing and the supporing files on Zenodo at https://zenodo.org/api/records/18722155/draft/files/O-MAP_data_processing.gz/content.
Place the files and the scripts in the respective folders (RNA_ChIP_data_processing and DNA_ChIP_data_processing), make the scripts executable, and download the fastq files for the RNA ChIP pulldown and input in a subfolder RNA_ChIP_data_processing/fastq

#### RNA ChIP data processing
RNA ChIP data anlysis consits of different steps requiring different docker images (specified below). All the steps can be run on the RNA_ChIP_data_processing folder

1) STAR geome index generation and alignment of input and pulldown data
   
       docker pull hedgelab/bulk_image:image3
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_STAR hedgelab/bulk_image:image3
       docker exec -it RNA_ChIP_STAR /bin/bash
       /home/shared_folder/STAR_alignment.sh
   
2) Conversion of BAM files to .bed tracks on the TTN genes, normalization on the coverage and generation of .bedGraph tracks with the log2FC pulldown/input

       docker exec -it RNA_ChIP_STAR /bin/bash
       ./bamToBed_input.sh
       ./bamToBed_pulldown.sh
       Rscript script_introns_alignment.R
   
3) rMATS analysis of differential splicing events from STAR BAM files of pulldown data and sashimi plots generation

       docker pull hedgelab/rmast:image1
       docker run -itv ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name RNA:ChIP_rMATS hedgelab/rmast:image1 /bin/bash
       /home/shared_folder/splicing_variants_analysis.sh
   
       docker pull hedgelab/hedgelab/sashimiplot:image1
       sudo docker run -d -v $PWD:/home/shared_folder --name RNA:ChIP_sashimi hedgelab/sashimiplot:image1
       /home/shared_folder/script_sashimi.sh
   
4) BWA alignment of input and pulldown data

       docker pull hedgelab/ciriquant:image1
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_CIRIquant hedgelab/ciriquant:image1
       docker exec -it RNA_ChIP_CIRIquant /bin/bash
       /home/shared_folder/BWA_alignment.sh
   
5) CIRI2 detection of circRNAs

       docker pull hedgelab/ciri_full:image3
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_CIRI2 hedgelab/ciri_full:image3
       docker exec -it RNA_ChIP_CIRI2 /bin/bash
       /home/shared_folder/CIRI2.sh
   
6) Quantification of differentially expressed circRNAs with CIRIquant

       docker exec -it RNA_ChIP_CIRIquant /bin/bash
       /home/shared_folder/CIRIquant.sh
   
7) O-MAP probes alignment and generation of .bed tracks for differentially enriched TTN introns and circRNAs or WT and ΔRBM20

       docker pull hedgelab/crispr:image1
       docker run -d -v ../RNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=RNA_ChIP_tracks hedgelab/crispr:image1
       docker exec -it RNA_ChIP_tracks /bin/bash
       Rscript /home/shared_folder/circRNAs_overlaps.R
    

#### DNA ChIP data processing
DNA ChIP data processing consists in the download of .bw track files from GEO and conversion to the .bedGraph format for circular plots generation and in the binning of chromosome 2 tracks at 100kb to enable comparison with other genomic tracks with Pearson correlation analysis.
The analysis can be run on our docker hedgelab/? and it canbe performed by running the following steps in the terminal from the DNA_ChIP_data_processing folder

     docker pull hedgelab/?
     docker run -d -v ../DNA_ChIP_data_processing:/home/shared_folder --privileged=true --name=DNA_ChIP hedgelab/?
     docker exec -it DNA_ChIP /bin/bash
     /home/shared_folder/chromatin_O-MAP_analysis.sh

## bulk data processing

## HiC data processing

## processed data analysis
