#!/usr/bin/env bash

export GTFDir="/Data/References/Annotations/Mus_musculus/GRCm38_ensembl93"
export TPMDir="/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_brain/processed_data/GRCm38_ensembl84_cDNA/kallisto"

if [ -d /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq ];
then
  cd /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq
else
  mkdir -p /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq
  cd /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq
fi


suppa.py generateEvents --input-file $GTFDir/Mus_musculus.GRCm38.93.gtf\
                        --output-file $GTFDir/Mus_musculus.GRCm38.93\
                        --format ioi &

suppa.py generateEvents --input-file $GTFDir/Mus_musculus.GRCm38.93.ERCC.gtf\
                        --output-file $GTFDir/Mus_musculus.GRCm38.93.ERCC.\
                        --format ioi &

suppa.py generateEvents --input-file $GTFDir/Mus_musculus.GRCm38.93.gtf\
                        --output-file $GTFDir/Mus_musculus.GRCm38.93\
                        --format ioe\
                        --event-type SE SS MX RI FL &

suppa.py generateEvents --input-file $GTFDir/Mus_musculus.GRCm38.93.ERCC.gtf\
                        --output-file $GTFDir/Mus_musculus.GRCm38.93.ERCC\
                        --format ioe\
                        --event-type SE SS MX RI FL &



suppa.py joinFiles --file-extension tpm\
                   --input-files 10_35_wt_HIPPO_CGATGT/abundance.tpm \
                                 1_12_wt_HIPPO_ACAGTG/abundance.tpm \
                                 10_42_wt_HIPPO_TGACCA/abundance.tpm \
                   --output HIPPO/WT/abundance

suppa.py joinFiles --file-extension tpm\
                   --input-files 10_16_hemi_HIPPO_CTTGTA/abundance.tpm \
                                 1_21_hemi_HIPPO_GCCAAT/abundance.tpm \
                                 10_28_hemi_HIPPO_CAGATC/abundance.tpm \
                   --output HIPPO/HEMI/abundance

suppa.py joinFiles --file-extension tpm\
                   --input-files 10_14_hemi_OB_GTCCGC/abundance.tpm \
                   10_26_hemi_OB_GTGAAA/abundance.tpm \
                   10_9_hemi_OB_CCGTCC/abundance.tpm \
                   --output OB/HEMI/abundance

suppa.py joinFiles --file-extension tpm\
                   --input-files 10_25_wt_OB_AGTTCC/abundance.tpm \
                   10_50_wt_OB_ATGTCA/abundance.tpm \
                   10_13_wt_OB_AGTCAA/abundance.tpm \
                   --output OB/WT/abundance

suppa.py joinFiles --file-extension tpm\
                   --input-files 1_12_wt_PFC_CGTACG/abundance.tpm \
                   10_35_wt_PFC_GTGGCC/abundance.tpm \
                   10_42_wt_PFC_GTTTCG/abundance.tpm \
                   --output PFC/WT/abundance

suppa.py joinFiles --file-extension tpm\
                   --input-files 10_28_hemi_PFC_ATTCCT/abundance.tpm \
                   10_27_hemi_PFC_ACTGAT/abundance.tpm \
                   1_21_hemi_PFC_GAGTGG/abundance.tpm \
                   --output PFC/HEMI/abundance

# PSI calculation per gene
suppa.py psiPerIsoform -g $GTFDir/Mus_musculus.GRCm38.93.gtf\
                       -e PFC/HEMI/abundance.tpm\
                       -o PFC/HEMI/results/PFC_HEMI\
                       -m INFO &

suppa.py psiPerIsoform -g $GTFDir/Mus_musculus.GRCm38.93.gtf\
                       -e PFC/WT/abundance.tpm\
                       -o PFC/WT/results/PFC_WT\
                       -m INFO &

suppa.py psiPerIsoform -g $GTFDir/Mus_musculus.GRCm38.93.gtf\
                       -e OB/HEMI/abundance.tpm\
                       -o OB/HEMI/results/OB_HEMI\
                       -m INFO &

suppa.py psiPerIsoform -g $GTFDir/Mus_musculus.GRCm38.93.gtf\
                       -e OB/WT/abundance.tpm\
                       -o OB/WT/results/OB_WT\
                       -m INFO &

suppa.py psiPerIsoform -g $GTFDir/Mus_musculus.GRCm38.93.gtf\
                       -e HIPPO/HEMI/abundance.tpm\
                       -o HIPPO/HEMI/results/HIPPO_HEMI\
                       -m INFO &

suppa.py psiPerIsoform -g $GTFDir/Mus_musculus.GRCm38.93.gtf\
                       -e HIPPO/WT/abundance.tpm\
                       -o HIPPO/WT/results/HIPPO_WT\
                       -m INFO &

# PSI per local event
for i in "PFC/HEMI" "PFC/WT" "OB/HEMI" "OB/WT" "HIPPO/WT" "HIPPO/HEMI"
do
  for ioe in $(find $GTFDir -name *.ioe)
  do
    event=$(echo $ioe | cut -f 5 -d "_")
    suppa.py psiPerEvent -i $ioe -e ${i}/abundance.tpm -o ${i}/results/${event} &
  done
done

# differential transcript usage WT - HEMI
for i in "PFC" "OB" "HIPPO"
do
  if [ ! -d ${i}/diff ]
  then
    mkdir -p ${i}/diff
  fi
  ioi=$GTFDir/Mus_musculus.GRCm38.93.ioi
  suppa.py diffSplice --method empirical\
                      --input $ioi\
                      --psi ${i}/WT/results/${i}_WT_isoform.psi ${i}/HEMI/results/${i}_HEMI_isoform.psi\
                      --tpm ${i}/WT/abundance.tpm ${i}/HEMI/abundance.tpm\
                      --area 1000 --lower-bound 0.05 -gc\
                      -o ${i}/diff/results &
done

# differential local events
for i in "PFC" "OB" "HIPPO"
do
  if [ ! -d ${i}/diff ]
  then
    mkdir -p ${i}/diff
  fi
  for ioe in $(find $GTFDir -name *.ioe)
  do
    event=$(echo $ioe | cut -f 5 -d "_")
    suppa.py diffSplice --method empirical\
                        --input $ioe\
                        --psi ${i}/WT/results/${event}.psi ${i}/HEMI/results/${event}.psi\
                        --tpm ${i}/WT/abundance.tpm ${i}/HEMI/abundance.tpm\
                        --area 1000 --lower-bound 0.05 -gc \
                        -o ${i}/diff/${event}
  done
done
