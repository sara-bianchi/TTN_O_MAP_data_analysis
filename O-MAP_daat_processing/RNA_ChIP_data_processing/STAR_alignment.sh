#!/bin/bash

#genome index creation
mkdir /home/shared_folder/genome
chmod 777 /home/shared_folder/genome
fasta=/home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
gtf=/home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf
/STAR --runThreadN 12 --runMode genomeGenerate --limitGenomeGenerateRAM=150000000000 --genomeDir /home/shared_folder/genome --genomeFastaFiles $fasta --sjdbGTFfile $gtf

#input STAR alignment
Rscript /home/shared_folder/Rscript_STAR_alignment.R /home/shared_folder/table_input.xlsx /home/shared_folder/output_STAR_input

#pulldown STAR alignment
Rscript /home/shared_folder/Rscript_STAR_alignment.R /home/shared_folder/table_pulldown.xlsx /home/shared_folder/output_STAR_pulldown
