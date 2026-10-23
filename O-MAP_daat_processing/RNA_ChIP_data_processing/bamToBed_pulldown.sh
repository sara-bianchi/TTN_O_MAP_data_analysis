#!/bin/bash
mkdir /home/shared_folder/output_STAR_pulldown/bed

/bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_pulldown/counts/KO_pulldown1/KO_pulldown1Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_pulldown/counts/KO_pulldown1/KO_pulldown1Aligned.sortedByCoord.out_bam2wig/KO_pulldown1Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_pulldown/bed/KO_pulldown1.bed
/bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_pulldown/counts/KO_pulldown2/KO_pulldown2Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_pulldown/counts/KO_pulldown2/KO_pulldown2Aligned.sortedByCoord.out_bam2wig/KO_pulldown2Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_pulldown/bed/KO_pulldown2.bed
/bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_pulldown/counts/KO_pulldown3/KO_pulldown3Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_pulldown/counts/KO_pulldown3/KO_pulldown3Aligned.sortedByCoord.out_bam2wig/KO_pulldown3Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_pulldown/bed/KO_pulldown3.bed
/bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_pulldown/counts/KO_pulldown4/KO_pulldown4Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_pulldown/counts/KO_pulldown4/KO_pulldown4Aligned.sortedByCoord.out_bam2wig/KO_pulldown4Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_pulldown/bed/KO_pulldown4.bed
/bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_pulldown/counts/WT_pulldown1/WT_pulldown1Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_pulldown/counts/WT_pulldown1/WT_pulldown1Aligned.sortedByCoord.out_bam2wig/WT_pulldown1Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_pulldown/bed/WT_pulldown1.bed
/bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_pulldown/counts/WT_pulldown2/WT_pulldown2Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_pulldown/counts/WT_pulldown2/WT_pulldown2Aligned.sortedByCoord.out_bam2wig/WT_pulldown2Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_pulldown/bed/WT_pulldown2.bed
/bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_pulldown/counts/WT_pulldown3/WT_pulldown3Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_pulldown/counts/WT_pulldown3/WT_pulldown3Aligned.sortedByCoord.out_bam2wig/WT_pulldown3Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_pulldown/bed/WT_pulldown3.bed
/bam2wig/bam2wig -c 2:178525989-178830802 /home/shared_folder/output_STAR_pulldown/counts/WT_pulldown4/WT_pulldown4Aligned.sortedByCoord.out.bam
wig2bed < /home/shared_folder/output_STAR_pulldown/counts/WT_pulldown4/WT_pulldown4Aligned.sortedByCoord.out_bam2wig/WT_pulldown4Aligned.sortedByCoord.out-c2_178525989-178830802.wig --multisplit bar > /home/shared_folder/output_STAR_pulldown/bed/WT_pulldown4.bed
