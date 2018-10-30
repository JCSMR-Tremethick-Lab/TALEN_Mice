library(data.table)
library(ggplot2)

suppaResults <- lapply(list.files("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/diff", pattern = "dpsi", full.names = T), function (x) fread(x))
names(suppaResults) <- list.files("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/diff", pattern = "dpsi")

suppaResults$results.dpsi.temp.0$ensembl_gene_id <- unlist(lapply(strsplit(suppaResults$results.dpsi.temp.0$Event_id, ";"), function(x) x[1]))
suppaResults$results.dpsi.temp.0[`WT_isoform-KO_isoform_p-val` < 0.1]

suppaResults <- lapply(suppaResults, function(x){
  x$ensembl_gene_id <- unlist(lapply(strsplit(x$Event_id, ";"), function(x) x[1]))
  return(x)
})

l1 <- lapply(suppaResults, function(x){
  table(x[,3] < 0.1)
})

sum(unlist(lapply(l1, function(x) {x["TRUE"]})))

l1 <- lapply(suppaResults, function(x){
  pval <- as.double(unlist(x[,3]))
  el <- length(unique(x[pval < 0.1]$ensembl_gene_id))
})

sum(unlist(l1))

el <- unique(unlist(l1))
length(el)

plotTargets <- suppaResults$results.dpsi.temp.0[`WT_isoform-KO_isoform_p-val` < 0.1]$ensembl_gene_id
plotTargets <- unique(plotTargets)
length(plotTargets)
ggplot(kTtx[plotTargets[12]][!is.na(tpm)], aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(width = 0.1, aes(shape = target_id, colour = target_id)) + 
  theme(legend.position = "right") +
  facet_wrap(~ external_gene_name)


# calculate the mean expression level, across all replicates and c --------
kTmean <- kT[, c("target_id", "tpm")]
kTmean <- kTmean[, mean(tpm), by = list(target_id)]
colnames(kTmean)[2] <- "tpm"
kTmean <- kTmean[!kTmean$tpm  == 0]
kTmean$tpm <- log10(kTmean$tpm + 0.01)
hist(kTmean$tpm)
setkey(kTmean, target_id)

# add mean expression to suppa results ------------------------------------
suppaResults <- lapply(suppaResults, function(i){
  merged <- merge(i, kTmean, by.x = "ensembl_gene_id", by.y = "target_id", all.x = F, all.y = F)
  return(merged)
})

suppaResults$RI.dpsi.temp.0[`RI-RI_p-val` < 0.1][tpm > 1]
suppaResults$MX.dpsi.temp.0[`MX-MX_p-val` < 0.1][tpm > 1]
suppaResults$SE.dpsi.temp.0[`SE-SE_p-val` < 0.1][tpm > 1][order(suppaResults$SE.dpsi.temp.0[`SE-SE_p-val` < 0.1][tpm > 1]$`SE-SE_dPSI`)]
order(suppaResults$SE.dpsi.temp.0[`SE-SE_p-val` < 0.1][tpm > 1]$`SE-SE_dPSI`)

DTUcondition <- suppaResults$results.dpsi.temp.0
DTUcondition <- merge(DTUcondition, kTmean, by.x = "ensembl_gene_id", by.y = "target_id", all.x = F, all.y = F)
DTUcondition <- DTUcondition[!is.na(`WT_isoform-KO_isoform_dPSI`)][!is.na(tpm)]

DTU <- fread("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/diff/results.psivec")
DTU
DTU$ensembl_gene_id <- unlist(lapply(strsplit(DTU$V1, ";"), function(x) x[1]))
DTU$ensembl_transcript_id <- unlist(lapply(strsplit(DTU$V1, ";"), function(x) x[2]))
setkey(kTtx, "transcript_id")
kTtx[DTU$ensembl_transcript_id]
kTtx$transcript_id <-  unlist(lapply(strsplit(kTtx$target_id, "\\."), function(x) x[1]))

DTUpairwise <- data.table(event = DTU$V1, psiPairwise = DTU$WT_isoform_1 - DTU$KO_isoform_1, pair = "pair1")
DTUpairwise <- rbind(data.table(event = DTU$V1, psiPairwise = DTU$WT_isoform_2 - DTU$KO_isoform_2, pair = "pair2"))
DTUpairwise <- rbind(data.table(event = DTU$V1, psiPairwise = DTU$WT_isoform_3 - DTU$KO_isoform_3, pair = "pair3"))
DTUpairwise$target_id <- unlist(lapply(strsplit(DTUpairwise$event, ";"), function(x) x[1]))
DTUpairwise <- merge(DTUpairwise, kTmean, by.x = "target_id", by.y = "target_id", all.x = F, all.y = F)
DTUpairwise <- DTUpairwise[!is.na(psiPairwise)][!is.na(tpm)]

# plotting a la SUPPA2 paper ----------------------------------------------
ggplot(DTUpairwise, aes(x = tpm, y = psiPairwise)) + 
  geom_point(colour = "lightgrey") +
  geom_point(colour = "darkgrey", data = DTUcondition, aes(x = tpm, y = `WT_isoform-KO_isoform_dPSI`)) +
  geom_point(colour = "darkblue", data = DTUcondition[`WT_isoform-KO_isoform_p-val` < 0.1], aes(x = tpm, y = `WT_isoform-KO_isoform_dPSI`)) + 
  xlab("TPM [log10(tpm) + 0.01]") +
  ylab(expression(paste(Delta, "PSI", sep = "")))


# export splicing candidates as BED file for deepTools plotting -----------
DTUcondition[`WT_isoform-KO_isoform_p-val` < 0.1]

lap1Genes <- fread("RS_lap_data copy.txt")
deepToolsUtils::WriteGRangesToBED(grGenes[el], 
                                  out_file = "/home/sebastian/Data/Tremethick/TALENs/Annotation/splicingCandidates.bed")

ensGenes[!el][unique(DTUcondition[`WT_isoform-KO_isoform_p-val` >= 0.1]$ensembl_gene_id)]

deepToolsUtils::WriteGRangesToBED(grGenes[ensGenes[!el][unique(DTUcondition[`WT_isoform-KO_isoform_p-val` >= 0.1]$ensembl_gene_id)]$ensembl_gene_id],
                                  out_file = "/home/sebastian/Data/Tremethick/TALENs/Annotation/expressedGenes.bed")

deepToolsUtils::WriteGRangesToBED(grGenes[ensGenes[!unique(DTUcondition$ensembl_gene_id)]$ensembl_gene_id],
                                  out_file = "/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/annotationFiles/silentGenes.bed")

  
deepToolsUtils::WriteGRangesToBED(grGenes[lap1Genes$V1],
                                  out_file = "/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/annotationFiles/Lap1Genes.bed")

rsLap1Genes <- readLines("RS_Lap1_genes.txt")
rsLap1Genes <- ensGenes[(ensGenes$external_gene_name %in% rsLap1Genes)]$ensembl_gene_id
deepToolsUtils::WriteGRangesToBED(grGenes[rsLap1Genes],
                                  out_file = "/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/annotationFiles/RSLap1Genes.bed")
  
lapply(suppaResults, function(x){
  pval <- as.double(unlist(x[,4]))
  el <- unique(x[pval < 0.1]$ensembl_gene_id)
  event <- unlist(strsplit(names(x)[4], "\\-"))[1]
  deepToolsUtils::WriteGRangesToBED(grGenes[el],
                                    out_file = paste("/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/annotationFiles/", event ,"Genes.bed", sep = ""))
  write.csv(x = x[el], file = paste(event, "_suppa_results.csv", sep = ""))
})

lapply(suppaResults, function(x){
  pval <- as.double(unlist(x[,4]))
  el <- unique(x[pval < 0.1]$ensembl_gene_id)
  event <- unlist(strsplit(names(x)[4], "\\-"))[1]
  intersect(el, rsLap1Genes)
})

lapply(suppaResults, function(x){
  pval <- as.double(unlist(x[,4]))
  el <- unique(x[pval < 0.1]$ensembl_gene_id)
  event <- unlist(strsplit(names(x)[4], "\\-"))[1]
  write.csv(x = x[intersect(el, lap1Genes$V1)], file = paste(event, "intersect_with_lap1_genes_suppa_results.csv", sep = ""))
})

