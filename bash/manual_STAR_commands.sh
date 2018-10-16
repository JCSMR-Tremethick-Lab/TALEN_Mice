STAR --runMode alignReads\
     --runThreadN 8\
     --genomeDir /Data/References/Transcriptomes/Mus_musculus/GRCm38_ensembl93/STAR_Index_ERCC\
     --readFilesIn RNA-Seq/trimmed/NB501086_0219_TSoboleva_JCSMR_RNAseq/KO_19_26.end1.fastq.gz RNA-Seq/trimmed/NB501086_0219_TSoboleva_JCSMR_RNAseq/KO_19_26.end2.fastq.gz\
     --readFilesCommand zcat\
     --outTmpDir /home/sebastian/tmp/KO_19_26\
     --outSAMmode Full\
     --outSAMattributes Standard\
     --outSAMtype BAM SortedByCoordinate\
     --alignEndsType EndToEnd\
     --outFileNamePrefix RNA-Seq/STAR/full/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/KO_19_26/ \
     1 > RNA-Seq/STAR/full/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/KO_19_26/log.txt
