# RUVg analysis
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

# Olfactory Bulb
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.ob.rda")
ktRawCounts.ob <- sleuth::kallisto_table(so.ob, normalized = F)
ktRawCounts.ob <- tidyr::spread(ktRawCounts.ob[, c("target_id", "sample", "est_counts")], sample, est_counts)
rownames(ktRawCounts.ob) <- ktRawCounts.ob$target_id
ktRawCounts.ob <- ktRawCounts.ob[,-1]
rm(so.ob)

# Pre-Frontal Cortex
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.pfc.rda")
ktRawCounts.pfc <- sleuth::kallisto_table(so.pfc, normalized = F)
ktRawCounts.pfc <- tidyr::spread(ktRawCounts.pfc[, c("target_id", "sample", "est_counts")], sample, est_counts)
rownames(ktRawCounts.pfc) <- ktRawCounts.pfc$target_id
ktRawCounts.pfc <- ktRawCounts.pfc[,-1]
rm(so.pfc)

# HIPPOcampus
load("../../NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/so.hippo.rda")
ktRawCounts.hippo <- sleuth::kallisto_table(so.hippo, normalized = F)
ktRawCounts.hippo <- tidyr::spread(ktRawCounts.hippo[, c("target_id", "sample", "est_counts")], sample, est_counts)
rownames(ktRawCounts.hippo) <- ktRawCounts.hippo$target_id
ktRawCounts.hippo <- ktRawCounts.hippo[,-1]
rm(so.hippo)

# merging all tables
ktRawCounts <- cbind(ktRawCounts.ob, ktRawCounts.pfc, ktRawCounts.hippo)


# filtering out transcripts with low read counts (min 2 samples with min 5 reads)
filter <- apply(counts, 1, function(x) length(x[x>5]) >= 2)
table(filter) # removes 50339
filtered <- counts[filter, ]
