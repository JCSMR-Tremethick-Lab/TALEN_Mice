# optimisation of peak calling parameters for macs2 analysis of cut & run data

dataDir <- ("/home/sebastian/Data/Tremethick/TALENs/CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/")
list.files(dataDir)
dir.create(paste(dataDir, "R_Analysis", sep = "/"))
setwd(paste(dataDir, "R_Analysis", sep = "/"))

source(file.path(dataDir, "KO_02_K36me3/KO_02_K36me3_model.r"))

# quick peak at the data to check FE
dT1 <- data.table::fread("/home/sebastian/Data/Tremethick/TALENs/CutRun/macs2/callpeak/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/WT_01_PolIIS5/WT_01_PolIIS5_peaks.xls")
summary(dT1$fold_enrichment)

