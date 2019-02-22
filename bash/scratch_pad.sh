java -Djava.io.tmpdir=/home/skurscheid/tmp \
     -Xmx40G \
     -jar /home/skurscheid/Bioinformatics/RNA-SeQC_v1.1.8.jar \
     -t /home/skurscheid/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf \
     -r /home/skurscheid/Data/References/Genomes/Mus_musculus/mm10_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa \
     -strat gc \
     -gc /home/skurscheid/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gc_content.txt \
     -s "10_28_hemi_PFC_ATTCCT|processed_data/GRCm38_ensembl84_ERCC/STAR/full/10_28_hemi_PFC_ATTCCT.aligned.bam|NA" \
     -gatkFlags --num_threads 16 \
     -o processed_data/GRCm38_ensembl84_ERCC/RNASeqQC/10_28_hemi_PFC_ATTCCT


java -Djava.io.tmpdir=/home/skurscheid/tmp \
      -Xmx36G \
      -jar /home/skurscheid/Bioinformatics/picard-tools-1.131/picard.jar CreateSequenceDictionary \
      REFERENCE=/home/skurscheid/Data/References/Genomes/Mus_musculus/mm10_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa \
      OUTPUT=/home/skurscheid/Data/References/Genomes/Mus_musculus/mm10_GRCm38/Mus_musculus.GRCm38.dna.toplevel.dict \
      GENOME_ASSEMBLY=mm10_GRCm38 \
      SPECIES=Mus_musculus


bowtie2 -k 200 --threads 8 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x /home/sebastian/Data/References/Transcriptomes/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.introns/trans_and_introns -1 processed_data/GRCm38_ensembl84/KMA_analysis/experiment/wt/wt1/wt1_R1.fastq.gz -2 processed_data/GRCm38_ensembl84/KMA_analysis/experiment/wt/wt1/wt1_R2.fastq.gz > ./test.bam

RNA-Seq/pizzly/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/WT_18_38/

/home/sebastian/miniconda3/envs/pizzly/bin/pizzly -k 31 --gtf /Data/References/Annotations/Mus_musculus/GRCm38_ensembl93/Mus_musculus.GRCm38.93.gtf --cache /Data/References/Annotations/Mus_musculus/GRCm38_ensembl93/Mus_musculus.GRCm38.93.gtf.cache.txt\
                                                  --fasta /Data/References/Transcriptomes/Mus_musculus/GRCm38_ensembl93/Mus_musculus.GRCm38.cds.all.fa.gz --output WT_18_38 RNA-Seq/kallisto/fusion/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/WT_18_38/fusion.txt

/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --nomodel \
                                                        --extsize 73\
                                                        --shift 37\
                                                        --broad \
                                                        --bdg \
                                                        --nolambda \
                                                        --verbose 3\
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3_35sec.bam\
                                                        --name KO_01_K36me3_35sec_nomodel_broad_nolambda\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3_35sec

sort -k1,1 -k2,2n KO_01_K36me3_35sec_nomodel_broad_control_lambda.bdg > KO_01_K36me3_35sec_nomodel_broad_control_lambda.sorted.bdg
sort -k1,1 -k2,2n KO_01_K36me3_35sec_nomodel_broad_treat_pileup.bdg > KO_01_K36me3_35sec_nomodel_broad_treat_pileup.sorted.bdg

bedGraphToBigWig KO_01_K36me3_35sec_nomodel_broad_control_lambda.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai KO_01_K36me3_35sec_nomodel_broad_control_lambda.bw &
bedGraphToBigWig KO_01_K36me3_35sec_nomodel_broad_treat_pileup.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai KO_01_K36me3_35sec_nomodel_broad_treat_pileup.bw &

sort -k1,1 -k2,2n KO_01_K36me3_35sec_control_lambda.bdg > KO_01_K36me3_35sec_control_lambda.sorted.bdg &
sort -k1,1 -k2,2n KO_01_K36me3_35sec_treat_pileup.bdg > KO_01_K36me3_35sec_treat_pileup.sorted.bdg &

bedGraphToBigWig KO_01_K36me3_35sec_control_lambda.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai KO_01_K36me3_35sec_control_lambda.bw &
bedGraphToBigWig KO_01_K36me3_35sec_treat_pileup.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai KO_01_K36me3_35sec_treat_pileup.bw &


sort -k1,1 -k2,2n KO_01_K36me3_35sec_nomodel_broad_control_lambda.bdg > KO_01_K36me3_35sec_nomodel_broad_control_lambda.sorted.bdg

bedGraphToBigWig KO_01_K36me3_35sec_nomodel_broad_control_lambda.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai KO_01_K36me3_35sec_nomodel_broad_control_lambda.bw &

# try with Pol2 data
# no model, no lambda
/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --nomodel \
                                                        --extsize 85\
                                                        --shift 43\
                                                        --broad \
                                                        --bdg \
                                                        --nolambda \
                                                        --verbose 3\
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5.bam\
                                                        --name WT_01_PolIIS5_nomodel_broad_nolambda\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5

sort -k1,1 -k2,2n CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_nolambda_treat_pileup.bdg > CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_nolambda_treat_pileup.sorted.bdg &
bedGraphToBigWig CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_nolambda_treat_pileup.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_nolambda_treat_pileup.bw &

# no model, with lambda
/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --nomodel \
                                                        --extsize 85\
                                                        --shift 43\
                                                        --broad \
                                                        --bdg \
                                                        --verbose 4\
                                                        --slocal 1200 \
                                                        --llocal 100000 \
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5.bam\
                                                        --name WT_01_PolIIS5_nomodel_broad\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5

sort -k1,1 -k2,2n CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_treat_pileup.bdg > CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_treat_pileup.sorted.bdg &
sort -k1,1 -k2,2n CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_control_lambda.bdg > CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_control_lambda.sorted.bdg

bedGraphToBigWig CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_treat_pileup.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_treat_pileup.bw &
bedGraphToBigWig CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_control_lambda.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_nomodel_broad_control_lambda.bw &


# lambda, model
/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --broad \
                                                        --bdg \
                                                        --verbose 4\
                                                        --slocal 1200 \
                                                        --llocal 100000 \
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5.bam\
                                                        --name WT_01_PolIIS5_broad\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5

sort -k1,1 -k2,2n CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_broad_treat_pileup.bdg > CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_broad_treat_pileup.sorted.bdg
bedGraphToBigWig CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_broad_treat_pileup.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_broad_treat_pileup.bw

#---------- Pol2 KO
# try with Pol2 data
# no model, no lambda
/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --nomodel \
                                                        --extsize 87\
                                                        --shift 45\
                                                        --broad \
                                                        --bdg \
                                                        --nolambda \
                                                        --verbose 3\
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5.bam\
                                                        --name KO_01_PolIIS5_nomodel_broad_nolambda\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5

sort -k1,1 -k2,2n CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_nolambda_treat_pileup.bdg > CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_nolambda_treat_pileup.sorted.bdg
bedGraphToBigWig CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_nolambda_treat_pileup.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_nolambda_treat_pileup.bw

# no model, with lambda
/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --nomodel \
                                                        --extsize 87\
                                                        --shift 45\
                                                        --broad \
                                                        --bdg \
                                                        --verbose 4\
                                                        --slocal 1200 \
                                                        --llocal 100000 \
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5.bam\
                                                        --name KO_01_PolIIS5_nomodel_broad\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5

sort -k1,1 -k2,2n CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_treat_pileup.bdg > CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_treat_pileup.sorted.bdg
sort -k1,1 -k2,2n CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_control_lambda.bdg > CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_control_lambda.sorted.bdg

bedGraphToBigWig CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_treat_pileup.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_treat_pileup.bw
bedGraphToBigWig CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_control_lambda.sorted.bdg /Data/References/Genomes/Mus_musculus/GRCm38_mm10/Mus_musculus.GRCm38.dna.toplevel.fa.fai CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5/KO_01_PolIIS5_nomodel_broad_control_lambda.bw

# model, with lambda
/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --broad \
                                                        --verbose 4\
                                                        --slocal 1200 \
                                                        --llocal 100000 \
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5.bam\
                                                        --name KO_01_PolIIS5_broad\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_PolIIS5
