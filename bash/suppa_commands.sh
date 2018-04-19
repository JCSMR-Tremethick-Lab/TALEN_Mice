#!/usr/bin/env bash

export GTFDir="/home/sebastian/Data/References/Annotations/Homo_sapiens/GRCh38_ensembl84"
export TPMDir="/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_brain/processed_data/GRCm38_ensembl84_cDNA/kallisto"
cd /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_brain/processed_data/GRCm38_ensembl84/suppa

suppa.py generateEvents --input-file $GTFDir/Homo_sapiens.GRCh38.84.gtf\
                        --output-file $GTFDir/Homo_sapiens.GRCh38.84.suppa\
                        --format ioi

suppa.py generateEvents --input-file $GTFDir/Homo_sapiens.GRCh38.84.gtf\
                        --output-file $GTFDir/Homo_sapiens.GRCh38.84.suppa\
                        --format ioe\
                        --event-type SE SS MX RI FL

suppa.py joinFiles --file-extension tpm\
                   --input-files $TPMDir/10_13_wt_OB_AGTCAA/abundance.tsv \
                                 $TPMDir/10_25_wt_OB_AGTTCC/abundance.tsv \
                                 $TPMDir/10_50_wt_OB_ATGTCA/abundance.tsv \
                   --output wt_OB
