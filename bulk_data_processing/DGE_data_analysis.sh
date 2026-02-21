#!/bin/bash

#genome index creation
mkdir /home/shared_folder/genome
chmod 777 /home/shared_folder/genome
fasta=/home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
gtf=/home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf
/STAR --runThreadN 12 --runMode genomeGenerate --limitGenomeGenerateRAM=150000000000 --genomeDir /home/shared_folder/genome --genomeFastaFiles $fasta --sjdbGTFfile $gtf

#RUES2 DGE analysis
Rscript /home/shared_folder/Rscript_bulk.R /home/shared_folder/settings_RUES2.xlsx /home/shared_folder/table_RUES2.xlsx /home/shared_folder/output_RUES2

#WTC11 DGE analysis
Rscript /home/shared_folder/Rscript_bulk.R /home/shared_folder/settings_WTC11.xlsx /home/shared_folder/table_WTC11.xlsx /home/shared_folder/output_WTC11

#KD DGE analysis
Rscript /home/shared_folder/Rscript_bulk.R /home/shared_folder/settings_KD.xlsx /home/shared_folder/table_KD.xlsx /home/shared_folder/output_KD
