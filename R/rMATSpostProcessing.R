# rMATS post processing
# Brain experiment 1
rMATS_dir <- "/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_brain/processed_data/GRCm38_ensembl84/rMATS/rMATS_PFC_wt_vs_hemi"
list.files(rMATS_dir)

RI.MATS.JC <- data.table::fread(file.path(rMATS_dir, "RI.MATS.JC.txt"))
summary(RI.MATS.JC$FDR)
RI.MATS.JCEC <- data.table::fread(file.path(rMATS_dir, "RI.MATS.JCEC.txt"))
A3SS.MATS.JC <- data.table::fread(file.path(rMATS_dir, "A3SS.MATS.JC.txt"))
A5SS.MATS.JC <- data.table::fread(file.path(rMATS_dir, "A5SS.MATS.JC.txt"))
SE.MATS.JC <- data.table::fread(file.path(rMATS_dir, "SE.MATS.JC.txt"))
MXE.MATS.JC <- data.table::fread(file.path(rMATS_dir, "MXE.MATS.JC.txt"))
table(SE.MATS.JC$FDR < 0.1)

table(SE.MATS.JC$FDR < 0.1)
table(RI.MATS.JC$FDR < 0.1)
table(A3SS.MATS.JC$FDR < 0.1)
table(A5SS.MATS.JC$FDR < 0.1)
table(MXE.MATS.JC$FDR < 0.1)

# testes
rMATS_dir <- "/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/processed_data/GRCm38_ensembl84/rMATS_fastq"
list.files(rMATS_dir)
RI.MATS.JC <- data.table::fread(file.path(rMATS_dir, "RI.MATS.JC.txt"))
summary(RI.MATS.JC$FDR)
RI.MATS.JCEC <- data.table::fread(file.path(rMATS_dir, "RI.MATS.JCEC.txt"))
A3SS.MATS.JC <- data.table::fread(file.path(rMATS_dir, "A3SS.MATS.JC.txt"))
A5SS.MATS.JC <- data.table::fread(file.path(rMATS_dir, "A5SS.MATS.JC.txt"))
SE.MATS.JC <- data.table::fread(file.path(rMATS_dir, "SE.MATS.JC.txt"))
MXE.MATS.JC <- data.table::fread(file.path(rMATS_dir, "MXE.MATS.JC.txt"))
table(SE.MATS.JC$FDR < 0.1)

table(SE.MATS.JC$FDR < 0.1)
table(RI.MATS.JC$FDR < 0.1)
table(A3SS.MATS.JC$FDR < 0.1)
table(A5SS.MATS.JC$FDR < 0.1)
table(MXE.MATS.JC$FDR < 0.1)

# quick check of DE genes
d <- "/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_brain_experiment_2/R_analysis/CSV_export"
list.files(d)
deData <- data.table::fread(file.path(d, "sleuth_results_wildtype_naive_vs_fear.csv"))
# 22 DEG naive wt vs ko
# 0 DEG fear wt vs ko <- high variability?
# 2 DEG mutant fear vs naive
# 280 DEG wildtype fear vs naive


                            
                            