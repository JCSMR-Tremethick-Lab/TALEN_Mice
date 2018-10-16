#-----------------------------------------
# H3K36me3 Data
#-----------------------------------------

/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --nomodel \
                                                        --extsize 73\
                                                        --shift 37\
                                                        --broad \
                                                        --nolambda \
                                                        --verbose 5\
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3.bam\
                                                        --name KO_01_K36me3_nomodel_broad_nolambda\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3

/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --nomodel \
                                                        --extsize 73\
                                                        --shift 37\
                                                        --broad \
                                                        --slocal 1200 \
                                                        --llocal 100000 \
                                                        --verbose 5\
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3.bam\
                                                        --name KO_01_K36me3_nomodel_broad\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3

/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --broad \
                                                        --verbose 3\
                                                        --slocal 1200 \
                                                        --llocal 100000 \
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3.bam\
                                                        --name KO_01_K36me3_broad\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3

#-----------------------------------------
# H3K36me3 Data - size filtered
#-----------------------------------------
/home/sebastian/miniconda3/envs/py27/bin/macs2 callpeak -f BAMPE\
                                                        --seed 1234\
                                                        --gsize hs\
                                                        --nomodel \
                                                        --bdg \
                                                        --extsize 73\
                                                        --shift 37\
                                                        --broad \
                                                        --nolambda \
                                                        --verbose 5\
                                                        --treatment CutRun/samtools/rmdup/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3_140bp.bam\
                                                        --name KO_01_K36me3_140bp_nomodel_broad_nolambda\
                                                        --outdir CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/KO_01_K36me3
