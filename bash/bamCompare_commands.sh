
cd /home/sebastian/Data/Tremethick/TALENs/CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun
export H3K36reference="WT_01_K36me3.bam"
export PolIIreference="WT_01_PolIIS5.bam"
export outdir="/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/bamCompare/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun"

for i in KO_01_K36me3.bam WT_02_K36me3.bam KO_02_K36me3.bam
do
  sample=$(echo $i|cut -f 1 -d ".");
  bamCompare --scaleFactorsMethod readCount\
             --operation second\
             --bamfile1 $H3K36reference\
             --bamfile2 $i\
             --exactScaling \
             --ignoreForNormalization X MT\
             --binSize 5 \
             --smoothLength 15 \
             --extendReads \
             --centerReads \
             --numberOfProcessors 24\
             --outFileName ${outdir}/${sample}_readCount.bw
done

for i in KO_01_K36me3.bam
do
  sample=$(echo $i|cut -f 1 -d ".");
  bamCompare --scaleFactorsMethod readCount\
             --operation first\
             --bamfile1 $H3K36reference\
             --bamfile2 $i\
             --exactScaling \
             --ignoreForNormalization X MT\
             --binSize 5 \
             --smoothLength 15 \
             --extendReads \
             --centerReads \
             --numberOfProcessors 32\
             --outFileName ${outdir}/WT_01_K36me3_readCount.bw
done


for i in KO_01_PolIIS5.bam KO_02_PolIIS5.bam WT_02_PolIIS5.bam
do
  sample=$(echo $i|cut -f 1 -d ".");
  bamCompare --scaleFactorsMethod readCount\
             --operation second\
             --bamfile1 $PolIIreference\
             --bamfile2 $i\
             --exactScaling \
             --ignoreForNormalization X MT\
             --binSize 5 \
             --smoothLength 15 \
             --extendReads \
             --centerReads \
             --numberOfProcessors 24\
             --outFileName ${outdir}/${sample}_readCount.bw
done

for i in KO_01_PolIIS5.bam
do
  sample=$(echo $i|cut -f 1 -d ".");
  bamCompare --scaleFactorsMethod readCount\
             --operation first\
             --bamfile1 $PolIIreference\
             --bamfile2 $i\
             --exactScaling \
             --ignoreForNormalization X MT\
             --binSize 5 \
             --smoothLength 15 \
             --extendReads \
             --centerReads \
             --numberOfProcessors 24\
             --outFileName ${outdir}/WT_01_PolIIS5_readCount.bw
done
