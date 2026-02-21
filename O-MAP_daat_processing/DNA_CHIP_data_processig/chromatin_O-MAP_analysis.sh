#!/bin/bash

mkdir /home/shared_folder/ChIP_track
wget -P /home/shared_folder/ChIP_track https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234808&format=file&file=GSE234808%5FKO%5FMock%5FRep1%5FL2FC%2Ebw
wget -P /home/shared_folder/ChIP_track https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234808&format=file&file=GSE234808%5FKO%5FMock%5FRep2%5FL2FC%2Ebw
wget -P /home/shared_folder/ChIP_track https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234808&format=file&file=GSE234808%5FKO%5FOMAP%5FRep1%5FL2FC%2Ebw
wget -P /home/shared_folder/ChIP_track https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234808&format=file&file=GSE234808%5FKO%5FOMAP%5FRep2%5FL2FC%2Ebw

wget -P /home/shared_folder/ChIP_track https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234808&format=file&file=GSE234808%5FWT%5FMock%5FRep1%5FL2FC%2Ebw
wget -P /home/shared_folder/ChIP_track https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234808&format=file&file=GSE234808%5FWT%5FMock%5FRep2%5FL2FC%2Ebw
wget -P /home/shared_folder/ChIP_track https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234808&format=file&file=GSE234808%5FWT%5FOMAP%5FRep1%5FL2FC%2Ebw
wget -P /home/shared_folder/ChIP_track https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234808&format=file&file=GSE234808%5FWT%5FOMAP%5FRep2%5FL2FC%2Ebw

#conversion from bigwig to bedgraph for R importation
bigWigToBedGraph /home/shared_folder/chromatin_data/pulldown_data/GSE234808_KO_Mock_Rep1_L2FC.bw /home/shared_folder/chromatin_data/pulldown_data/RBM20KO_mock_R1.bedGraph
bigWigToBedGraph /home/shared_folder/chromatin_data/pulldown_data/GSE234808_KO_Mock_Rep2_L2FC.bw /home/shared_folder/chromatin_data/pulldown_data/RBM20KO_mock_R2.bedGraph
bigWigToBedGraph /home/shared_folder/chromatin_data/pulldown_data/GSE234808_KO_OMAP_Rep1_L2FC.bw /home/shared_folder/chromatin_data/pulldown_data/RBM20KO_OMAP_R1.bedGraph
bigWigToBedGraph /home/shared_folder/chromatin_data/pulldown_data/GSE234808_KO_OMAP_Rep2_L2FC.bw /home/shared_folder/chromatin_data/pulldown_data/RBM20KO_OMAP_R2.bedGraph
bigWigToBedGraph /home/shared_folder/chromatin_data/pulldown_data/GSE234808_WT_Mock_Rep1_L2FC.bw /home/shared_folder/chromatin_data/pulldown_data/WT_mock_R1.bedGraph
bigWigToBedGraph /home/shared_folder/chromatin_data/pulldown_data/GSE234808_WT_Mock_Rep2_L2FC.bw /home/shared_folder/chromatin_data/pulldown_data/WT_mock_R2.bedGraph
bigWigToBedGraph /home/shared_folder/chromatin_data/pulldown_data/GSE234808_WT_OMAP_Rep1_L2FC.bw /home/shared_folder/chromatin_data/pulldown_data/WT_OMAP_R1.bedGraph
bigWigToBedGraph /home/shared_folder/chromatin_data/pulldown_data/GSE234808_WT_OMAP_Rep2_L2FC.bw /home/shared_folder/chromatin_data/pulldown_data/WT_OMAP_R2.bedGraph

#tracks binning for pearson correlation
binSize=100000
cd /home/shared_folder/ChIP_track
sort-bed WT_OMAP_R1.bedGraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > signal.bed
bedops --merge signal.bed | bedops --chop ${binSize} - > bins.bed
bedmap --echo-ref-name --mean bins.bed signal.bed > WT_OMAP_R1_binned100000.bed

sort-bed WT_OMAP_R2.bedGraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > signal.bed
bedops --merge signal.bed | bedops --chop ${binSize} - > bins.bed
bedmap --echo-ref-name --mean bins.bed signal.bed > WT_OMAP_R2_binned100000.bed

sort-bed RBM20KO_OMAP_R1.bedGraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > signal.bed
bedops --merge signal.bed | bedops --chop ${binSize} - > bins.bed
bedmap --echo-ref-name --mean bins.bed signal.bed > RBM20KO_OMAP_R1_binned100000.bed

sort-bed RBM20KO_OMAP_R2.bedGraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > signal.bed
bedops --merge signal.bed | bedops --chop ${binSize} - > bins.bed
bedmap --echo-ref-name --mean bins.bed signal.bed > RBM20KO_OMAP_R2_binned100000.bed

sort-bed WT_mock_R1.bedGraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > signal.bed
bedops --merge signal.bed | bedops --chop ${binSize} - > bins.bed
bedmap --echo-ref-name --mean bins.bed signal.bed > WT_mock_R1_binned100000.bed

sort-bed WT_mock_R2.bedGraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > signal.bed
bedops --merge signal.bed | bedops --chop ${binSize} - > bins.bed
bedmap --echo-ref-name --mean bins.bed signal.bed > WT_mock_R2_binned100000.bed

sort-bed RBM20KO_mock_R1.bedGraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > signal.bed
bedops --merge signal.bed | bedops --chop ${binSize} - > bins.bed
bedmap --echo-ref-name --mean bins.bed signal.bed > RBM20KO_mock_R1_binned100000.bed

sort-bed RBM20KO_mock_R2.bedGraph | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > signal.bed
bedops --merge signal.bed | bedops --chop ${binSize} - > bins.bed
bedmap --echo-ref-name --mean bins.bed signal.bed > RBM20KO_mock_R2_binned100000.bed


