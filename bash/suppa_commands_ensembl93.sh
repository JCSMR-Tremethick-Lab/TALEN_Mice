#!/usr/bin/env bash

suppa.py joinFiles --file-extension tpm\
                   --input-files RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/WT_18_38/abundance.tpm RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/WT_37_39/abundance.tpm RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/WT_46_47/abundance.tpm\
                   --output RNA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/WT/abundance.tpm

suppa.py joinFiles --file-extension tpm\
                   --input-files RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/KO_19_26/abundance.tpm RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/KO_24_25/abundance.tpm RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/KO_44_45/abundance.tpm\
                   --output RNA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/KO/abundance.tpm


export GTFDir="/Data/References/Annotations/Mus_musculus/GRCm38_ensembl93"
export TPMDir="/home/sebastian/Data/Tremethick/TALENs/NA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/"

if [ -d /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq ];
then
  cd /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq
else
  mkdir -p /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq
  cd /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq
fi


suppa.py generateEvents --input-file $GTFDir \
                        --output-file $GTFDir/Mus_musculus.GRCm38.93\
                        --format ioi &

suppa.py generateEvents --input-file $GTFDir\
                        --output-file $GTFDir/Mus_musculus.GRCm38.93.ERCC.\
                        --format ioi &

suppa.py generateEvents --input-file $GTFDir\
                        --output-file $GTFDir/Mus_musculus.GRCm38.93\
                        --format ioe\
                        --event-type SE SS MX RI FL &

suppa.py generateEvents --input-file $GTFDir\
                        --output-file $GTFDir/Mus_musculus.GRCm38.93.ERCC\
                        --format ioe\
                        --event-type SE SS MX RI FL &

# PSI calculation per gene
cd /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq

suppa.py psiPerIsoform -g $GTFDir/Mus_musculus.GRCm38.93.ERCC.gtf\
                       -e WT/abundance.tpm\
                       -o WT/results/WT\
                       -m INFO &

suppa.py psiPerIsoform -g $GTFDir/Mus_musculus.GRCm38.93.ERCC.gtf\
                       -e KO/abundance.tpm\
                       -o KO/results/KO\
                       -m INFO &


# PSI per local event
for i in "WT" "KO"
do
  for ioe in $(find $GTFDir -name *ERCC*.ioe)
  do
    event=$(echo $ioe | cut -f 5 -d "_")
    suppa.py psiPerEvent -i $ioe -e ${i}/abundance.tpm -o ${i}/results/${event} &
  done
done

# differential transcript usage WT - HEMI

ioi=$GTFDir/Mus_musculus.GRCm38.93.ERCC.ioi
suppa.py diffSplice --method empirical\
                    --input $ioi\
                    --psi WT/results/WT_isoform.psi KO/results/KO_isoform.psi\
                    --tpm WT/abundance.tpm KO/abundance.tpm\
                    --area 1000 --lower-bound 0.05 -gc\
                    -o diff/results


# differential local events

for ioe in $(find $GTFDir -name *ERCC*.ioe)
do
  event=$(echo $ioe | cut -f 5 -d "_")
  suppa.py diffSplice --method empirical\
                      --input $ioe\
                      --psi WT/results/${event}.psi KO/results/${event}.psi\
                      --tpm WT/abundance.tpm KO/abundance.tpm\
                      --area 1000 --lower-bound 0.05 -gc \
                      -o diff/${event}
done
