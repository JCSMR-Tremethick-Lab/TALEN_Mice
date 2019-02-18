#!/bin/bash
conda activate deeptools

cd /home/sebastian/Data/Tremethick/TALENs

# extracting separate groups from computeMatrix output prior to plotHeatmap
export matrixLocation="CutRun/deepTools/computeMatrix/scale-region/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/allGenes/"
computeMatrixOperations info -m $matrixLocation/matrix_RPKM.gz
	# WT_01_K36me3_RPKM
	# WT_02_K36me3_RPKM
	# WT_01_PolIIS5_RPKM
	# WT_02_PolIIS5_RPKM
	# KO_01_K36me3_RPKM
	# KO_02_K36me3_RPKM
	# KO_01_PolIIS5_RPKM
	# KO_02_PolIIS5_RPKM
computeMatrixOperations subset -m $matrixLocation/matrix_RPKM.gz --outFileName $matrixLocation/matrix_RPKM_Pol2.gz --samples WT_01_PolIIS5_RPKM WT_02_PolIIS5_RPKM KO_01_PolIIS5_RPKM KO_02_PolIIS5_RPKM &
computeMatrixOperations subset -m $matrixLocation/matrix_RPKM.gz --outFileName $matrixLocation/matrix_RPKM_H3K36me3.gz --samples WT_01_K36me3_RPKM WT_02_K36me3_RPKM KO_01_K36me3_RPKM KO_02_K36me3_RPKM &

export outFile="CutRun/deepTools/plotHeatmap/scale-region/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/allGenes_RPKM_Pol2_kmeans3"
/home/sebastian/miniconda3/envs/deepTools/bin/plotHeatmap --matrixFile $matrixLocation/matrix_RPKM_Pol2.gz \
                                                          --outFileName $outFile.pdf \
                                                          --outFileSortedRegions $outFile.bed \
                                                          --outFileNameMatrix $outFile.tab \
                                                          --kmeans 3 \
                                                          --samplesLabel WT_01_PolIIS5 WT_02_PolIIS5 KO_01_PolIIS5 KO_02_PolIIS5 \
                                                          --labelRotation -90 &

export outFileK36me3="CutRun/deepTools/plotHeatmap/scale-region/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/allGenes_RPKM_H3K36me3_kmeans3"
/home/sebastian/miniconda3/envs/deepTools/bin/plotHeatmap --matrixFile $matrixLocation/matrix_RPKM_H3K36me3.gz\
                                                          --outFileName $outFileK36me3.pdf\
                                                          --outFileSortedRegions $outFileK36me3.bed\
                                                          --outFileNameMatrix $outFileK36me3.tab\
                                                          --kmeans 3 \
                                                          --samplesLabel WT_01_K36me3 WT_02_K36me3 KO_01_K36me3 KO_02_K36me3 \
                                                          --labelRotation -90 &

# reference-point matrices
export matrixLocation="CutRun/deepTools/computeMatrix/reference-point/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"
computeMatrixOperations info -m $matrixLocation/allExons/matrix_RPKM.gz
	# WT_01_K36me3_RPKM
	# WT_02_K36me3_RPKM
	# WT_01_PolIIS5_RPKM
	# WT_02_PolIIS5_RPKM
	# KO_01_K36me3_RPKM
	# KO_02_K36me3_RPKM
	# KO_01_PolIIS5_RPKM
	# KO_02_PolIIS5_RPKM
computeMatrixOperations subset -m $matrixLocation/allExons/matrix_RPKM.gz --outFileName $matrixLocation/allExons/matrix_RPKM_Pol2.gz --samples WT_01_PolIIS5_RPKM WT_02_PolIIS5_RPKM KO_01_PolIIS5_RPKM KO_02_PolIIS5_RPKM &
computeMatrixOperations subset -m $matrixLocation/allExons/matrix_RPKM.gz --outFileName $matrixLocation/allExons/matrix_RPKM_H3K36me3.gz --samples WT_01_K36me3_RPKM WT_02_K36me3_RPKM KO_01_K36me3_RPKM KO_02_K36me3_RPKM &
computeMatrixOperations subset -m $matrixLocation/allIntrons/matrix_RPKM.gz --outFileName $matrixLocation/allIntrons/matrix_RPKM_Pol2.gz --samples WT_01_PolIIS5_RPKM WT_02_PolIIS5_RPKM KO_01_PolIIS5_RPKM KO_02_PolIIS5_RPKM &
computeMatrixOperations subset -m $matrixLocation/allIntrons/matrix_RPKM.gz --outFileName $matrixLocation/allIntrons/matrix_RPKM_H3K36me3.gz --samples WT_01_K36me3_RPKM WT_02_K36me3_RPKM KO_01_K36me3_RPKM KO_02_K36me3_RPKM &

export outFilePol2="CutRun/deepTools/plotHeatmap/reference-point/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/allExons_RPKM_Pol2_kmeans4"
/home/sebastian/miniconda3/envs/deepTools/bin/plotHeatmap --matrixFile $matrixLocation/allExons/matrix_RPKM_Pol2.gz \
                                                          --outFileName $outFilePol2.pdf \
                                                          --outFileSortedRegions $outFilePol2.bed \
                                                          --outFileNameMatrix $outFilePol2.tab \
                                                          --kmeans 4 \
                                                          --samplesLabel WT_01_PolIIS5 WT_02_PolIIS5 KO_01_PolIIS5 KO_02_PolIIS5&

export outFilePol2="CutRun/deepTools/plotHeatmap/reference-point/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/allIntrons_RPKM_Pol2_kmeans4"
/home/sebastian/miniconda3/envs/deepTools/bin/plotHeatmap --matrixFile $matrixLocation/allExons/matrix_RPKM_Pol2.gz \
                                                          --outFileName $outFilePol2.pdf \
                                                          --outFileSortedRegions $outFilePol2.bed \
                                                          --outFileNameMatrix $outFilePol2.tab \
                                                          --kmeans 4 \
                                                          --samplesLabel WT_01_PolIIS5 WT_02_PolIIS5 KO_01_PolIIS5 KO_02_PolIIS5&

export outFileK36me3="CutRun/deepTools/plotHeatmap/scale-region/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/allGenes_RPKM_H3K36me3_kmeans3"
/home/sebastian/miniconda3/envs/deepTools/bin/plotHeatmap --matrixFile $matrixLocation/matrix_RPKM_H3K36me3.gz\
                                                          --outFileName $outFileK36me3.pdf\
                                                          --outFileSortedRegions $outFileK36me3.bed\
                                                          --outFileNameMatrix $outFileK36me3.tab\
                                                          --kmeans 3 \
                                                          --samplesLabel WT_01_K36me3 WT_02_K36me3 KO_01_K36me3 KO_02_K36me3 \
                                                          --labelRotation -90 &
