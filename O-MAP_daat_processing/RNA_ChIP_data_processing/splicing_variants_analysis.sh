#!/bin/bash

#docker run -itv /30tb/3tb/data/sara/RSEM:/home/shared_folder --name rmast hedgelab/rmast:image1 /bin/bash

#rmats analysis pulldown
python /rmats/rmats.py --b1 /home/shared_folder/rmats_pulldown/samples_WT.txt --b2 /home/shared_folder/rmats_pulldown/samples_RBM20KO.txt --gtf /home/shared_folder/genome_files/Homo_sapiens.GRCh38.113.gtf -t paired --readLength 250 --allow-clipping --variable-read-length --nthread 20 --individual-counts --od /home/shared_folder/rmats_pulldown/WT_vs_RBM20KO --tmp /temp1 

