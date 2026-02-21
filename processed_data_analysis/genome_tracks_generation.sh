#!/bin/bash

#create output folder
mkdir /home/shared_folder/output_plots

#pyGenomeTracks
#figure 2 - panel C
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_OMAP --region chr2:1-242193529 --outFileName /home/shared_folder/output_plots/CIS_pulldown_OMAP.pdf --fontSize 7
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_chromatin --region chr2:1-242193529 --outFileName /home/shared_folder/output_plots/com_virtual4C.pdf --fontSize 7
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_CIDs --region chr2:178000000-179500000 --outFileName /home/shared_folder/output_plots/CID_6362.pdf --fontSize 5 --width 15 --trackLabelFraction 0.1 
#figure 2 - panel D
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_TIDs --region chr9:131000000-138000000 --outFileName /home/shared_folder/output_plots/TID_130_242.pdf --fontSize 5 --width 15 --trackLabelFraction 0.1 
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_TIDs --region chr11:66000000-73000000 --outFileName /home/shared_folder/output_plots/TID_178_75.pdf --fontSize 5 --width 15 --trackLabelFraction 0.1 
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_TIDs --region chr16:84000000-91000000 --outFileName /home/shared_folder/output_plots/TID_245.pdf --fontSize 5 --width 15 --trackLabelFraction 0.1 
#figure 3 - panel D
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_compartments --region chr11:71056000-71815000 --outFileName /home/shared_folder/output_plots/TID_75.pdf --fontSize 7 --width 20
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_compartments --region chr17:78488000-80956000 --outFileName /home/shared_folder/output_plots/TID_141.pdf --fontSize 7 --width 20
#supplementary figure 8 panel B
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_RUES2_TTN --region chr2:178500000-178850000 --outFileName /home/shared_folder/output_plots/TTN_track.pdf --fontSize 7 --width 20
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/track_file_RUES2_RBM20 --region chr10:110600000-110900000 --outFileName /home/shared_folder/output_plots/RBM20_track.pdf --fontSize 7 --width 20
#supplementary figure 13 - panel A
pyGenomeTracks --tracks /home/shared_folder/genome_track_files/comparison_track --region chr2:178500000-178850000 --outFileName /home/shared_folder/output_plots/TTN_probes_overlaps.pdf --fontSize 7 --width 20



