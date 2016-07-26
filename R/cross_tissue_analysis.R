# tissue map of transcribed histones in mouse/human,
# with transcribed promoter activity levels
# 1) compile list of all histone genes
# 2) obtain FANTOM5 sequencing data
# 3) obtain matching RNA-Seq data from mouse
# 4) create tissue specific coverage maps for each histone gene

amILocal <- function(machinename = NULL){
  if(is.null(machinename)) stop("Machinename is missing")
  m <- Sys.info()["nodename"]
  mn <- unlist(lapply(strsplit(m, "\\."), function(x) x[1]))
  if (mn == machinename) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# setting working directory and data sources ------------------------------
if (amILocal("JCSMR027564ML")){
  setwd("~/mount/gduserv/Data/Tremethick/TALENs/meta_analysis/R_analysis")
  base_dir <- "~/mount/gduserv/Data/Tremethick/TALENs/meta_analysis/WT"
  load("~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
} else {
  setwd("~/Data/Tremethick/TALENs/meta_analysis/R_analysis")
  base_dir <- "~/Data/Tremethick/TALENs/meta_analysis/WT"
  load("~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
}

# load pre-processed data
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.ob.rda")
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.pfc.rda")
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.hippo.rda")
