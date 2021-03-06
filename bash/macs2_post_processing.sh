#!/bin/bash
conda activate deepTools
export macs2dir="/home/sebastian/Data/Tremethick/TALENs/CutRun/macs2/callpeak/broad/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"
export macs2broad="/home/sebastian/Data/Tremethick/TALENs/CutRun/macs2/callpeak/broad/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"

export bw2dir="/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/bamCompare/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"

# analyse WT_01_PolIIS5
cd /home/sebastian/Data/Tremethick/TALENs/CutRun/macs2/post_processing/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5

# coverage over actual Pol2 peaks
computeMatrix reference-point -S $bw2dir/WT_01_PolIIS5_readCount.bw \
                              -R $macs2dir/WT_01_PolIIS5/WT_01_PolIIS5_summits.bed \
                              -p 16 \
                              --smartLabels \
                              --referencePoint center \
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000 \
                              --outFileName WT_01_PolIIS5_macs2_narrow_peak_coverage.gz

plotProfile --matrixFile WT_01_PolIIS5_macs2_narrow_peak_coverage.gz \
            --outFileName WT_01_PolIIS5_macs2_narrow_peak_coverage.pdf \
            --plotType se

# broadPeak -> scale-regions
computeMatrix scale-regions -S $bw2dir/WT_01_PolIIS5_readCount.bw \
                              -R $macs2broad/WT_01_PolIIS5/WT_01_PolIIS5_peaks.broadPeak \
                              -p 16 \
                              --smartLabels \
                              --regionBodyLength 500 \
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000 \
                              --startLabel "Broad peak start" \
                              --endLabel "Broad peak end" \
                              --outFileName WT_01_PolIIS5_macs2_broad_peak_coverage.gz

plotProfile --matrixFile WT_01_PolIIS5_macs2_broad_peak_coverage.gz \
            --outFileName WT_01_PolIIS5_macs2_broad_peak_coverage.pdf \
            --plotType se

# broadPeak -> reference-point
computeMatrix reference-point -S $bw2dir/WT_01_PolIIS5_readCount.bw \
                              -R $macs2broad/WT_01_PolIIS5/WT_01_PolIIS5_peaks.broadPeak \
                              -p 16 \
                              --smartLabels \
                              --referencePoint center \
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000 \
                              --outFileName WT_01_PolIIS5_macs2_broad_peak_coverage_reference_point.gz

plotProfile --matrixFile WT_01_PolIIS5_macs2_broad_peak_coverage_reference_point.gz \
            --outFileName WT_01_PolIIS5_macs2_broad_peak_coverage_reference_point.pdf \
            --plotType se


# check what it looks like at K36me3 peaks
computeMatrix reference-point -S $bw2dir/WT_01_PolIIS5_readCount.bw \
                              -R $macs2dir/WT_01_K36me3/WT_01_K36me3_summits.bed \
                              -p 16 \
                              --smartLabels \
                              --referencePoint center \
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000 \
                              --outFileName WT_01_PolIIS5_macs2_H3K36me3_narrow_peak_coverage.gz

plotProfile --matrixFile WT_01_PolIIS5_macs2_H3K36me3_narrow_peak_coverage.gz \
            --outFileName WT_01_PolIIS5_macs2_H3K36me3_narrow_peak_coverage.pdf \
            --plotType se

#------------------------
# analyse WT_02_PolIIS5
cd /home/sebastian/Data/Tremethick/TALENs/CutRun/macs2/post_processing/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_02_PolIIS5

# narrowPeak -> reference-point
computeMatrix reference-point -S $bw2dir/WT_02_PolIIS5_readCount.bw \
                              -R $macs2dir/WT_02_PolIIS5/WT_02_PolIIS5_summits.bed \
                              -p 16 \
                              --smartLabels \
                              --referencePoint center \
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000 \
                              --outFileName WT_02_PolIIS5_macs2_narrow_peak_coverage.gz

plotProfile --matrixFile WT_02_PolIIS5_macs2_narrow_peak_coverage.gz \
            --outFileName WT_02_PolIIS5_macs2_narrow_peak_coverage.pdf \
            --plotType se

# broadPeak -> scale-regions
computeMatrix scale-regions -S $bw2dir/WT_02_PolIIS5_readCount.bw \
                              -R $macs2broad/WT_02_PolIIS5/WT_02_PolIIS5_peaks.broadPeak \
                              -p 16 \
                              --smartLabels \
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000 \
                              --startLabel "Broad peak start" \
                              --endLabel "Broad peak end" \
                              --outFileName WT_02_PolIIS5_macs2_broad_peak_coverage.gz

plotProfile --matrixFile WT_02_PolIIS5_macs2_broad_peak_coverage.gz \
            --outFileName WT_02_PolIIS5_macs2_broad_peak_coverage.pdf \
            --plotType se
