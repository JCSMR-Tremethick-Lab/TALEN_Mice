# within tissue analysis of mouse brain RNA-seq data
require(rhdf5)
require(sleuth)
require(biomaRt)
require(tidyr)
require(rtracklayer)


source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

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
# Olfactory Bulb
# 12 DE transcripts @ 10% FDR
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.ob.rda")
rt.ob <- sleuth::sleuth_results(so.ob, "conditionhemi_OB")
rownames(rt.ob) <- rt.ob$target_id
kt.ob <- sleuth::kallisto_table(so.ob)
kt.ob <- tidyr::spread(kt.ob[, c("target_id", "sample", "tpm")], sample, tpm)
rownames(kt.ob) <- kt.ob$target_id
length(which(rt.ob$qval <= 0.1))
rm(so.ob)

# Pre-Frontal Cortex
# 6 DE transcripts @ 10% FDR
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.pfc.rda")
rt.pfc <- sleuth::sleuth_results(so.pfc, "conditionhemi_PFC")
rownames(rt.pfc) <- rt.pfc$target_id
kt.pfc <- sleuth::kallisto_table(so.pfc)
kt.pfc <- tidyr::spread(kt.pfc[, c("target_id", "sample", "tpm")], sample, tpm)
rownames(kt.pfc) <- kt.pfc$target_id
length(which(rt.pfc$qval <= 0.1))
rm(so.pfc)

# HIPPOcampus
# 1181 DE transcripts @ 10% FDR
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.hippo.rda")
rt.hippo <- sleuth::sleuth_results(so.hippo, "conditionhemi_HIPPO")
rownames(rt.hippo) <- rt.hippo$target_id
kt.hippo <- sleuth::kallisto_table(so.hippo)
kt.hippo <- tidyr::spread(kt.hippo[, c("target_id", "sample", "tpm")], sample, tpm)
rownames(kt.hippo) <- kt.hippo$target_id
length(which(rt.hippo$qval <= 0.1))


# look at H2A.Bbd/Lap1 expression -----------------------------------------
# look at H2A.Lap1 genes
toGrep <- c("H2afb3", "H2afb2", "Gm14920")
toGrep <- t2g[grep(paste(toGrep, collapse = "|"), t2g$ext_gene),]$target_id
names(toGrep) <- c("H2afb3", "H2afb2", "Gm14920")
kt.goi <- cbind(kt.pfc[grep(paste(toGrep, collapse = "|"), kt.pfc$target_id),],
                kt.hippo[grep(paste(toGrep, collapse = "|"), kt.hippo$target_id),],
                kt.ob[grep(paste(toGrep, collapse = "|"), kt.ob$target_id),])
kt.goi <- kt.goi[, c(1:7, 9:14, 16:21)]

kt.goi <- gather(kt.goi, sample, value = tpm, c(2:19))
kt.goi$gene <- NA
kt.goi$tissue <- NA

for (i in 1:3){
  kt.goi[which(kt.goi$target_id == toGrep[i]), ]$gene <- names(toGrep)[i]
}

for (i in unique(kt.goi$sample)){
  kt.goi[which(kt.goi$sample == i), ]$tissue <- unlist(lapply(strsplit(i, "_"), function(x) paste(x[3:4], collapse = "_")))
}

