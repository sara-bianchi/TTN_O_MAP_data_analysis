#!/bin/bash

mkdir /home/shared_folder/output_STAR_input/bed

./bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_input/counts/KO_input1/KO_input1Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_input/counts/KO_input1/KO_input1Aligned.sortedByCoord.out_bam2wig/KO_input1Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_input/bed/KO_input1.bed
./bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_input/counts/KO_input2/KO_input2Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_input/counts/KO_input2/KO_input2Aligned.sortedByCoord.out_bam2wig/KO_input2Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_input/bed/KO_input2.bed
./bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_input/counts/KO_input3/KO_input3Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_input/counts/KO_input3/KO_input3Aligned.sortedByCoord.out_bam2wig/KO_input3Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_input/bed/KO_input3.bed
./bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_input/counts/KO_input4/KO_input4Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_input/counts/KO_input4/KO_input4Aligned.sortedByCoord.out_bam2wig/KO_input4Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_input/bed/KO_input4.bed
./bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_input/counts/WT_input1/WT_input1Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_input/counts/WT_input1/WT_input1Aligned.sortedByCoord.out_bam2wig/WT_input1Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_input/bed/WT_input1.bed
./bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_input/counts/WT_input2/WT_input2Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_input/counts/WT_input2/WT_input2Aligned.sortedByCoord.out_bam2wig/WT_input2Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_input/bed/WT_input2.bed
./bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_input/counts/WT_input3/WT_input3Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_input/counts/WT_input3/WT_input3Aligned.sortedByCoord.out_bam2wig/WT_input3Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_input/bed/WT_input3.bed
./bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_input/counts/WT_input4/WT_input4Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_input/counts/WT_input4/WT_input4Aligned.sortedByCoord.out_bam2wig/WT_input4Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_input/bed/WT_input4.bed
