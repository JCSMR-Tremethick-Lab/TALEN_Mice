library(data.table)

suppaResults <- lapply(list.files("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/diff", pattern = "dpsi", full.names = T), function (x) fread(x))
names(suppaResults) <- list.files("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/suppa/pooled/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq/diff", pattern = "dpsi")

suppaResults$results.dpsi.temp.0$ensembl_gene_id <- unlist(lapply(strsplit(suppaResults$results.dpsi.temp.0$Event_id, ";"), function(x) x[1]))
suppaResults$results.dpsi.temp.0[`WT_isoform-KO_isoform_p-val` < 0.05]

lapply(suppaResults, function(x){
  x[ensembl_gene_id == "ENSMUSG00000000581"]
})

suppaResults <- lapply(suppaResults, function(x){
  x$ensembl_gene_id <- unlist(lapply(strsplit(x$Event_id, ";"), function(x) x[1]))
  return(x)
})

plotTargets <- suppaResults$results.dpsi.temp.0[`WT_isoform-KO_isoform_p-val` < 0.1]$ensembl_gene_id
plotTargets <- unique(plotTargets)
length(plotTargets)
ggplot(kTtx[plotTargets[12]][!is.na(tpm)], aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(width = 0.1, aes(shape = target_id, colour = target_id)) + 
  theme(legend.position = "right") +
  facet_wrap(~ external_gene_name)

kTtx[plotData[1:3]$ensembl_gene_id]
