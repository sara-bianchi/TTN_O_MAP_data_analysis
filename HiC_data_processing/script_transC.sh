#!/bin/bash

#format conversion to h5
mkdir /home/shared_folder/h5_raw

hicConvertFormat --matrices /home/shared_folder/ooutput_HiC/mcool_raw/WT_D2.mcool --inputFormat cool --outputFormat h5 --outFileName /home/shared_folder/h5_raw/WT_D2.h5
hicConvertFormat --matrices /home/shared_folder/ooutput_HiC/mcool_raw/WT_D3.mcool --inputFormat cool --outputFormat h5 --outFileName /home/shared_folder/h5_raw/WT_D3.h5
hicConvertFormat --matrices /home/shared_folder/ooutput_HiC/mcool_raw/KO_D2.mcool --inputFormat cool --outputFormat h5 --outFileName /home/shared_folder/h5_raw/KO_D2.h5
hicConvertFormat --matrices /home/shared_folder/ooutput_HiC/mcool_raw/KO_D3.mcool --inputFormat cool --outputFormat h5 --outFileName /home/shared_folder/h5_raw/KO_D3.h5

#sum the replicates
hicSumMatrices --matrices home/shared_folder/h5_raw/WT_D2.h5 home/shared_folder/h5_raw/WT_D3.h5 --outFileName /home/shared_folder/h5_raw/WT.h5
hicSumMatrices --matrices home/shared_folder/h5_raw/KO_D2.h5 home/shared_folder/h5_raw/KO_D3.h5 --outFileName /home/shared_folder/h5_raw/KO.h5

#convert in mcool
mkdir /home/shared_folder/transC
mkdir /home/shared_folder/transC/input_mcool
hicConvertFormat --matrices  /home/shared_folder/h5_raw/WT.h5 --inputFormat h5 --outputFormat cool --outFileName  /home/shared_folder/transC/input_mcool/WT.mcool
hicConvertFormat --matrices  /home/shared_folder/h5_raw/KO.h5 --inputFormat h5 --outputFormat cool --outFileName  /home/shared_folder/transC/input_mcool/KO.mcool

#run transC
mkdir /home/shared_folder/transC/output_WT
mkdir /home/shared_folder/transC/output_KO
python /trans-C/code/trans_C.py /home/shared_folder/transC/input_mcool/WT.mcool /home/shared_folder/genome_files/chrom_hg38.sizes 1000000 /home/shared_folder/seeds_TTN.bed /home/shared_folder/transC/output_WT/
python /trans-C/code/trans_C.py /home/shared_folder/transC/input_mcool/KO.mcool /home/shared_folder/genome_files/chrom_hg38.sizes 1000000 /home/shared_folder/seeds_TTN.bed /home/shared_folder/transC/output_KO/

#run random model generation
python /home/shared_folder/python_WT.py
python /home/shared_folder/python_KO.py
