#!/bin/bash

#docker run -itv /30tb/3tb/data/sara/RSEM:/home/shared_folder --name rmast hedgelab/rmast:image1 /bin/bash

#rmats analysis WTC11
python /rmats/rmats.py --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KO.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 250 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO --tmp /temp1 

python /rmats/rmats.py --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_TTNKO.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 250 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_WTC11/WT_vs_TTNKO --tmp /temp2 

python /rmats/rmats.py --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KOTTNKO.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 250 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO --tmp /temp3 

#rmats analysis RUES2
python /rmats/rmats.py --b1 /home/shared_folder/rmats_RUES2/samples_WT.txt --b2 /home/shared_folder/rmats_RUES2/samples_RBM20KO.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 185 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_RUES2/WT_vs_RBM20KO --tmp /temp4

python /rmats/rmats.py --b1 /home/shared_folder/rmats_RUES2/samples_WT.txt --b2 /home/shared_folder/rmats_RUES2/samples_TTNpromKO.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 185 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_RUES2/WT_vs_TTNpromKO --tmp /temp5

python /rmats/rmats.py --b1 /home/shared_folder/rmats_RUES2/samples_WT.txt --b2 /home/shared_folder/rmats_RUES2/samples_TTNKO.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 185 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_RUES2/WT_vs_TTNKO --tmp /temp6

#rmats analysis KD
python /rmats/rmats.py --b1 /home/shared_folder/rmats_KD/samples_SCR.txt --b2 /home/shared_folder/rmats_KD/samples_QKI.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 275 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_KD/SCR_vs_QKI --tmp /temp7

python /rmats/rmats.py --b1 /home/shared_folder/rmats_KD/samples_SCR.txt --b2 /home/shared_folder/rmats_KD/samples_SAFB.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 275 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_KD/SCR_vs_SAFB --tmp /temp8

python /rmats/rmats.py --b1 /home/shared_folder/rmats_KD/samples_SCR.txt --b2 /home/shared_folder/rmats_KD/samples_ELAVL1.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 275 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_KD/SCR_vs_ELAVL1 --tmp /temp9
