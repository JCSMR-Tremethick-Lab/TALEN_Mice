export computeMatrixDir="/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/computeMatrix/reference-point/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"

source activate deepTools
cd /home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/computeMatrix/scale-region/GRCm38_ensembl93/publicationPlots

plotHeatmap --matrixFile $computeMatrixDir/allExons/matrix_1xgenome.gz\
            --outFileName NB501086_0221_TSoboleva_JCSMR_CutandRun/allExons_1xgenome.pdf &

plotHeatmap --matrixFile $computeMatrixDir/allExons/matrix_RPKM.gz\
            --outFileName NB501086_0221_TSoboleva_JCSMR_CutandRun/allExons_RPKM.pdf &

plotHeatmap --matrixFile $computeMatrixDir/allExons/matrix_RPGCExact.gz\
            --outFileName NB501086_0221_TSoboleva_JCSMR_CutandRun/allExons_RPGCExact.pdf &

export computeMatrixDir="/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/computeMatrix/scale-region/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"


plotHeatmap --matrixFile $computeMatrixDir/allGenes/matrix_RPGCExact.gz\
            --outFileName NB501086_0221_TSoboleva_JCSMR_CutandRun/allGenes_kmeans4_RPGCExact.pdf\
            --kmeans 4 &

export regionsFileDir="/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/annotationFiles"
export scoreFileDir="/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/bamCoverage/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"
export scoreFileDir2="/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/bamCoverage/GRCm38_ensembl93/180731_NB501086_0217_CutandRun_Tanya/"

computeMatrix scale-regions --regionsFileName $regionsFileDir/splicingCandidates.bed $regionsFileDir/expressedGenes.bed $regionsFileDir/silentGenes.bed \
                            --scoreFileName $scoreFileDir/WT_01_PolIIS5_RPGCExact.bw $scoreFileDir/WT_02_PolIIS5_RPGCExact.bw $scoreFileDir/WT_01_K36me3_RPGCExact.bw $scoreFileDir/WT_02_K36me3_RPGCExact.bw $scoreFileDir/KO_01_PolIIS5_RPGCExact.bw $scoreFileDir/KO_02_PolIIS5_RPGCExact.bw $scoreFileDir/KO_01_K36me3_RPGCExact.bw $scoreFileDir/KO_02_K36me3_RPGCExact.bw \
                            --upstream 2000 \
                            --downstream 2000 \
                            --samplesLabel "PolIIS5 WT1" "PolIIS5 WT2" "PolIIS5 KO1" "PolIIS5 KO2" "H3K36me3 WT1" "H3K36me3 WT2" "H3K36me3 KO1" "H3K36me3 KO2" \
                            --regionBodyLength 5000 \
                            --numberOfProcessors 24 \
                            --outFileName groupedGenes.gz

plotHeatmap --matrixFile groupedGenes.gz \
            --outFileName splicing_expressed_silent.pdf

computeMatrix scale-regions --regionsFileName $regionsFileDir/splicingCandidates.bed $regionsFileDir/expressedGenes.bed \
                            --scoreFileName $scoreFileDir/WT_01_PolIIS5_RPGCExact.bw $scoreFileDir/WT_02_PolIIS5_RPGCExact.bw $scoreFileDir/KO_01_PolIIS5_RPGCExact.bw $scoreFileDir/KO_02_PolIIS5_RPGCExact.bw $scoreFileDir/WT_01_K36me3_RPGCExact.bw $scoreFileDir/WT_02_K36me3_RPGCExact.bw $scoreFileDir/KO_01_K36me3_RPGCExact.bw $scoreFileDir/KO_02_K36me3_RPGCExact.bw \
                            --upstream 2000 \
                            --downstream 2000 \
                            --samplesLabel "PolIIS5 WT1" "PolIIS5 WT2" "PolIIS5 KO1" "PolIIS5 KO2" "H3K36me3 WT1" "H3K36me3 WT2" "H3K36me3 KO1" "H3K36me3 KO2" \
                            --regionBodyLength 5000 \
                            --numberOfProcessors 24 \
                            --outFileName expressedGenes.gz

plotHeatmap --matrixFile expressedGenes.gz \
            --outFileName splicing_expressed.pdf

computeMatrix scale-regions --regionsFileName $regionsFileDir/splicingCandidates.bed \
                            --scoreFileName $scoreFileDir/WT_01_PolIIS5_RPGCExact.bw $scoreFileDir/WT_02_PolIIS5_RPGCExact.bw $scoreFileDir/KO_01_PolIIS5_RPGCExact.bw $scoreFileDir/KO_02_PolIIS5_RPGCExact.bw $scoreFileDir/WT_01_K36me3_RPGCExact.bw $scoreFileDir/WT_02_K36me3_RPGCExact.bw $scoreFileDir/KO_01_K36me3_RPGCExact.bw $scoreFileDir/KO_02_K36me3_RPGCExact.bw $scoreFileDir2/WT_01_H2AL2_3_7_18_RPGCExact.bw $scoreFileDir2/WT_02_H2AL2_24_6_18_RPGCExact.bw $scoreFileDir2/KO_01_H2AL2_3_7_18_RPGCExact.bw $scoreFileDir2/KO_02_H2AL2_24_6_18_RPGCExact.bw \
                            --upstream 2000 \
                            --downstream 2000 \
                            --samplesLabel "PolIIS5 WT1" "PolIIS5 WT2" "PolIIS5 KO1" "PolIIS5 KO2" "H3K36me3 WT1" "H3K36me3 WT2" "H3K36me3 KO1" "H3K36me3 KO2" "H2AL2 WT01" "H2AL2 WT02" "H2AL2 KO01" "H2AL2 KO02"\
                            --regionBodyLength 5000 \
                            --numberOfProcessors 8 \
                            --outFileName diffSplicedGenes_w_H2AL.gz

plotHeatmap --matrixFile diffSplicedGenes_w_H2AL.gz \
            --outFileName spliced_w_H2AL.pdf

plotHeatmap --matrixFile diffSplicedGenes.gz \
            --outFileName spliced.pdf


computeMatrix scale-regions --regionsFileName $regionsFileDir/Lap1Genes.bed $regionsFileDir/RSLap1Genes.bed\
                            --scoreFileName $scoreFileDir/WT_01_PolIIS5_RPGCExact.bw $scoreFileDir/WT_02_PolIIS5_RPGCExact.bw $scoreFileDir/KO_01_PolIIS5_RPGCExact.bw $scoreFileDir/KO_02_PolIIS5_RPGCExact.bw $scoreFileDir/WT_01_K36me3_RPGCExact.bw $scoreFileDir/WT_02_K36me3_RPGCExact.bw $scoreFileDir/KO_01_K36me3_RPGCExact.bw $scoreFileDir/KO_02_K36me3_RPGCExact.bw \
                            --upstream 2000 \
                            --downstream 2000 \
                            --samplesLabel "PolIIS5 WT1" "PolIIS5 WT2" "PolIIS5 KO1" "PolIIS5 KO2" "H3K36me3 WT1" "H3K36me3 WT2" "H3K36me3 KO1" "H3K36me3 KO2" \
                            --regionBodyLength 5000 \
                            --numberOfProcessors 24 \
                            --outFileName Lap1Genes.gz

plotHeatmap --matrixFile Lap1Genes.gz \
            --outFileName Lap1Genes.pdf

computeMatrix scale-regions --regionsFileName $regionsFileDir/Lap1Genes.bed $regionsFileDir/RSLap1Genes.bed\
                            --scoreFileName $scoreFileDir/WT_01_PolIIS5_RPGCExact.bw $scoreFileDir/WT_02_PolIIS5_RPGCExact.bw $scoreFileDir/KO_01_PolIIS5_RPGCExact.bw $scoreFileDir/KO_02_PolIIS5_RPGCExact.bw $scoreFileDir/WT_01_K36me3_RPGCExact.bw $scoreFileDir/WT_02_K36me3_RPGCExact.bw $scoreFileDir/KO_01_K36me3_RPGCExact.bw $scoreFileDir/KO_02_K36me3_RPGCExact.bw $scoreFileDir2/WT_01_H2AL2_3_7_18_RPGCExact.bw $scoreFileDir2/WT_02_H2AL2_24_6_18_RPGCExact.bw $scoreFileDir2/KO_01_H2AL2_3_7_18_RPGCExact.bw $scoreFileDir2/KO_02_H2AL2_24_6_18_RPGCExact.bw \
                            --upstream 2000 \
                            --downstream 2000 \
                            --samplesLabel "PolIIS5 WT1" "PolIIS5 WT2" "PolIIS5 KO1" "PolIIS5 KO2" "H3K36me3 WT1" "H3K36me3 WT2" "H3K36me3 KO1" "H3K36me3 KO2" "H2AL2 WT01" "H2AL2 WT02" "H2AL2 KO01" "H2AL2 KO02"\
                            --regionBodyLength 5000 \
                            --numberOfProcessors 24 \
                            --outFileName Lap1Genes_w_H2AL.gz

plotHeatmap --matrixFile Lap1Genes_w_H2AL.gz \
            --outFileName Lap1Genes_w_H2AL.pdf


computeMatrix scale-regions --regionsFileName $regionsFileDir/A3Genes.bed $regionsFileDir/A5Genes.bed $regionsFileDir/AFGenes.bed $regionsFileDir/ALGenes.bed $regionsFileDir/MXGenes.bed $regionsFileDir/RIGenes.bed $regionsFileDir/SEGenes.bed $regionsFileDir/DTUGenes.bed\
                            --scoreFileName $scoreFileDir/WT_01_PolIIS5_RPGCExact.bw $scoreFileDir/WT_02_PolIIS5_RPGCExact.bw $scoreFileDir/KO_01_PolIIS5_RPGCExact.bw $scoreFileDir/KO_02_PolIIS5_RPGCExact.bw $scoreFileDir/WT_01_K36me3_RPGCExact.bw $scoreFileDir/WT_02_K36me3_RPGCExact.bw $scoreFileDir/KO_01_K36me3_RPGCExact.bw $scoreFileDir/KO_02_K36me3_RPGCExact.bw $scoreFileDir2/WT_01_H2AL2_3_7_18_RPGCExact.bw $scoreFileDir2/WT_02_H2AL2_24_6_18_RPGCExact.bw $scoreFileDir2/KO_01_H2AL2_3_7_18_RPGCExact.bw $scoreFileDir2/KO_02_H2AL2_24_6_18_RPGCExact.bw \
                            --upstream 2000 \
                            --downstream 2000 \
                            --samplesLabel "PolIIS5 WT1" "PolIIS5 WT2" "PolIIS5 KO1" "PolIIS5 KO2" "H3K36me3 WT1" "H3K36me3 WT2" "H3K36me3 KO1" "H3K36me3 KO2" "H2AL2 WT01" "H2AL2 WT02" "H2AL2 KO01" "H2AL2 KO02"\
                            --regionBodyLength 5000 \
                            --numberOfProcessors 8 \
                            --outFileName SUPPAEvents_w_H2AL.gz

plotHeatmap --matrixFile SUPPAEvents_w_H2AL.gz \
            --outFileName SUPPAEvents_w_H2AL.pdf

for i in DTUGenes #A3Genes A5Genes AFGenes ALGenes MXGenes RIGenes SEGenes
do
  computeMatrix scale-regions --regionsFileName $regionsFileDir/${i}.bed\
                              --scoreFileName $scoreFileDir/WT_01_PolIIS5_RPGCExact.bw $scoreFileDir/WT_02_PolIIS5_RPGCExact.bw $scoreFileDir/KO_01_PolIIS5_RPGCExact.bw $scoreFileDir/KO_02_PolIIS5_RPGCExact.bw $scoreFileDir/WT_01_K36me3_RPGCExact.bw $scoreFileDir/WT_02_K36me3_RPGCExact.bw $scoreFileDir/KO_01_K36me3_RPGCExact.bw $scoreFileDir/KO_02_K36me3_RPGCExact.bw $scoreFileDir2/WT_01_H2AL2_3_7_18_RPGCExact.bw $scoreFileDir2/WT_02_H2AL2_24_6_18_RPGCExact.bw $scoreFileDir2/KO_01_H2AL2_3_7_18_RPGCExact.bw $scoreFileDir2/KO_02_H2AL2_24_6_18_RPGCExact.bw \
                              --upstream 2000 \
                              --downstream 2000 \
                              --samplesLabel "PolIIS5 WT1" "PolIIS5 WT2" "PolIIS5 KO1" "PolIIS5 KO2" "H3K36me3 WT1" "H3K36me3 WT2" "H3K36me3 KO1" "H3K36me3 KO2" "H2AL2 WT01" "H2AL2 WT02" "H2AL2 KO01" "H2AL2 KO02"\
                              --regionBodyLength 5000 \
                              --numberOfProcessors 8 \
                              --outFileName ${i}.gz
  plotHeatmap --matrixFile ${i}.gz \
              --outFileName ${i}.pdf
  gdrive upload ${i}.pdf
done
