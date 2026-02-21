# TTN_O_MAP_data_analysis
scripts for reproducibility of the bioinforatic analysis of "Nascent trascript O-MAP revels the molecular architecture of a single-locus subnuclear compartment organized around RBM20 and the TTN RNA", Kania, Bianchi et al.

Our analyses include:

the alignment and preprocessing of raw data:
- O-MAP data processing: RNA ChIP pulldown and input data for WT and ΔRBM20 hESCs RUES2 derived CMs, where raw reads are available on GEO with the accession number GSE234808
- bulk data processing: bulk RNA seq of hESCs RUES2 derived CMs (WT, ΔRBM20, ΔTTNpromoter, ΔTTNprotein), hiPSCs WTC11 derived CMs (WT, ΔRBM20, ΔTTNpromoter, ΔTTNpromoter-ΔRBM20), RUES2 KD CMs (Scramble, SAFB KD, ELAVL1 KD, QKI KD), where raw data are available on ArrayExprress with the following accession numbers, respecively
- HiC data processing: HiC data processing of RUES2 hESCs including the developmental trajectory of cardiac differentiation (ESCs day 0, MES day 2, CP day 5 and CMs day 14) available on 4DN Data Portal (Experiment Set codes 4DNESGTHHJAC, 4DNES6G787KF, 4DNESSU7RRNR, 4DNESW72J62T), and WT and ΔRBM20 day 30 CMs available on ENA at the following accession number

further analyses and integration of the previous datasets with other published data starting from intermediate preprocessed files obtained in the previous steps and available on Zenodo at the following DOI 10.5281/zenodo.18722155. The repository also contains some settings and envrionment files for the processing of raw data.
It is highly recommanded to run all our scripts using our docker images available on the hedgelab docker repository to ensure reproducibility.

### O-MAP data processing
# RNA ChIP data processing
# DNA ChIP data processing

## bulk data processing

## HiC data processing

## processed data analysis
