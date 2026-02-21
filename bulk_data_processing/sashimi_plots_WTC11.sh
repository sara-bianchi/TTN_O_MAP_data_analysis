#!/bin/bash

#docker run -d -v $PWD:/home/shared_folder --name sashimi hedgelab/sashimiplot:image1

mkdir /home/shared_folder/sashimi_WTC11
mkdir /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KO
mkdir /home/shared_folder/sashimi_WTC11/WT_vs_TTNKO
mkdir /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KOTTNKO

#sashimi plots
rmats2sashimiplot --l1 WT --l2 RBM20KO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO/A3SS.MATS_sign.JC.txt --event-type A3SS -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KO/A3SS --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO.gf
rmats2sashimiplot --l1 WT --l2 TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_TTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_TTNKO/A3SS.MATS_sign.JC.txt --event-type A3SS -o /home/shared_folder/sashimi_WTC11/WT_vs_TTNKO/A3SS --group-info /home/shared_folder/rmats_WTC11/WT_vs_TTNKO.gf
rmats2sashimiplot --l1 WT --l2 RBM20KO_TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KOTTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO/A3SS.MATS_sign.JC.txt --event-type A3SS -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KOTTNKO/A3SS --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO.gf

rmats2sashimiplot --l1 WT --l2 RBM20KO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO/A5SS.MATS_sign.JC.txt --event-type A5SS -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KO/A5SS --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO.gf
rmats2sashimiplot --l1 WT --l2 TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_TTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_TTNKO/A5SS.MATS_sign.JC.txt --event-type A5SS -o /home/shared_folder/sashimi_WTC11/WT_vs_TTNKO/A5SS --group-info /home/shared_folder/rmats_WTC11/WT_vs_TTNKO.gf
rmats2sashimiplot --l1 WT --l2 RBM20KO_TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KOTTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO/A5SS.MATS_sign.JC.txt --event-type A5SS -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KOTTNKO/A5SS --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO.gf

rmats2sashimiplot --l1 WT --l2 RBM20KO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO/RI.MATS_sign.JC.txt --event-type RI -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KO/RI --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO.gf
rmats2sashimiplot --l1 WT --l2 TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_TTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_TTNKO/RI.MATS_sign.JC.txt --event-type RI -o /home/shared_folder/sashimi_WTC11/WT_vs_TTNKO/RI --group-info /home/shared_folder/rmats_WTC11/WT_vs_TTNKO.gf
rmats2sashimiplot --l1 WT --l2 RBM20KO_TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KOTTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO/RI.MATS_sign.JC.txt --event-type RI -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KOTTNKO/RI --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO.gf

rmats2sashimiplot --l1 WT --l2 RBM20KO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO/MXE.MATS_sign.JC.txt --event-type MXE -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KO/MXE --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO.gf
rmats2sashimiplot --l1 WT --l2 TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_TTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_TTNKO/MXE.MATS_sign.JC.txt --event-type MXE -o /home/shared_folder/sashimi_WTC11/WT_vs_TTNKO/MXE --group-info /home/shared_folder/rmats_WTC11/WT_vs_TTNKO.gf
rmats2sashimiplot --l1 WT --l2 RBM20KO_TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KOTTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO/MXE.MATS_sign.JC.txt --event-type MXE -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KOTTNKO/MXE --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO.gf

rmats2sashimiplot --l1 WT --l2 RBM20KO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO/SE.MATS_IMMT_sign.JC.txt --event-type SE -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KO/SE --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KO.gf
rmats2sashimiplot --l1 WT --l2 TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_TTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_TTNKO/SE.MATS_IMMT_sign.JC.txt --event-type SE -o /home/shared_folder/sashimi_WTC11/WT_vs_TTNKO/SE --group-info /home/shared_folder/rmats_WTC11/WT_vs_TTNKO.gf
rmats2sashimiplot --l1 WT --l2 RBM20KO_TTNpromKO --b1 /home/shared_folder/rmats_WTC11/samples_WT.txt --b2 /home/shared_folder/rmats_WTC11/samples_RBM20KOTTNKO.txt -e /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO/SE.MATS_IMMT_sign.JC.txt --event-type SE -o /home/shared_folder/sashimi_WTC11/WT_vs_RBM20KOTTNKO/SE --group-info /home/shared_folder/rmats_WTC11/WT_vs_RBM20KOTTNKO.gf

