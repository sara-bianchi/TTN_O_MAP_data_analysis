#!/bin/bash

#input BWA alignment
mkdir BWA_input
/BWA_input/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/KO-Input1-RNA_R1_001.fastq.gz /home/shared_folder/fastq/KO-Input1-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_input/KO_Input1.sam
/BWA_input/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/KO-Input2-RNA_R1_001.fastq.gz /home/shared_folder/fastq/KO-Input2-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_input/KO_Input2.sam
/BWA_input/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/KO-Input3-RNA_R1_001.fastq.gz /home/shared_folder/fastq/KO-Input3-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_input/KO_Input3.sam
/BWA_input/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/KO-Input4-RNA_R1_001.fastq.gz /home/shared_folder/fastq/KO-Input4-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_input/KO_Input4.sam

/BWA_input/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/TTN-Input1-RNA_R1_001.fastq.gz /home/shared_folder/fastq/TTN-Input1-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_input/WT_Input1.sam
/BWA_input/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/TTN-Input2-RNA_R1_001.fastq.gz /home/shared_folder/fastq/TTN-Input2-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_input/WT_Input2.sam
/BWA_input/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/TTN-Input3-RNA_R1_001.fastq.gz /home/shared_folder/fastq/TTN-Input3-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_input/WT_Input3.sam
/BWA_input/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/TTN-Input4-RNA_R1_001.fastq.gz /home/shared_folder/fastq/TTN-Input4-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_input/WT_Input4.sam

#pulldown BWA alignment
mkdir BWA_pulldown
/BWA_pulldown/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/TTN-Pulldown1-RNA_R1_001.fastq.gz /home/shared_folder/fastq/TTN-Pulldown1-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_pulldown/WT_pulldown1.sam
/BWA_pulldown/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/TTN-Pulldown2-RNA_R1_001.fastq.gz /home/shared_folder/fastq/TTN-Pulldown2-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_pulldown/WT_pulldown2.sam
/BWA_pulldown/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/TTN-Pulldown3-RNA_R1_001.fastq.gz /home/shared_folder/fastq/TTN-Pulldown3-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_pulldown/WT_pulldown3.sam
/BWA_pulldown/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/TTN-Pulldown4-RNA_R1_001.fastq.gz /home/shared_folder/fastq/TTN-Pulldown4-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_pulldown/WT_pulldown4.sam

/BWA_pulldown/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/KO-Pulldown1-RNA_R1_001.fastq.gz /home/shared_folder/fastq/KO-Pulldown1-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_pulldown/KO_pulldown1.sam
/BWA_pulldown/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/KO-Pulldown2-RNA_R1_001.fastq.gz /home/shared_folder/fastq/KO-Pulldown2-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_pulldown/KO_pulldown2.sam
/BWA_pulldown/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/KO-Pulldown3-RNA_R1_001.fastq.gz /home/shared_folder/fastq/KO-Pulldown3-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_pulldown/KO_pulldown3.sam
/BWA_pulldown/bwa mem /home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/shared_folder/fastq/KO-Pulldown4-RNA_R1_001.fastq.gz /home/shared_folder/fastq/KO-Pulldown4-RNA_R2_001.fastq.gz > /home/shared_folder/BWA_pulldown/KO_pulldown4.sam
