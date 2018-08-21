export computeMatrixDir="/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/computeMatrix/reference-point/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"

source activate deepTools
cd /home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/plotHeatmap

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
