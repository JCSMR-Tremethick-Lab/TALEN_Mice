require("DEXSeq")
require(biomaRt)
require(parallel)
#
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

# defining global variables
ensemblHost <- "mar2016.archive.ensembl.org"

if (amILocal("JCSMR027564ML")){
 BPPARAM <- MulticoreParam(workers = 7)
  setwd("~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/DEXSeq/count/"
  if (file.exists("~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.DEXSeq.gtf")) {
    flattenedfile <- "~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.DEXSeq.gtf"
  } else {
    stop("GTF file missing!")
  }
} else {
  BPPARAM <- MulticoreParam(workers = 16)
  setwd("~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/DEXSeq/count/"
  flattenedfile <- "~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.DEXSeq.gtf"
}

# biomaRt connection
mouse <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
attribs <- biomaRt::listAttributes(mouse)

# annotation data
if (!file.exists("ens84Exons.rda")){
  ens84Exons <- getBM(attributes = c("ensembl_gene_id",
                                      "ensembl_exon_id",
                                      "chromosome_name",
                                      "exon_chrom_start",
                                      "exon_chrom_end",
                                      "is_constitutive",
                                      "rank"),
                      mart = mouse)
  ens84Exons <- ens84Exons[order(ens84Exons$ensembl_gene_id),]
  save(ens84Exons, file = "ens84Exons.rda")
} else {
  load("ens84Exons.rda")
}

# first have to integrate the size factor estimation from RUVg in order to use ERCC spike ins
# effectively W can be added to the GLM...

files <- list.files(path = dataPath, full.names = T)
names(files) <- unlist(lapply(strsplit(list.files(path = dataPath, full.names = F), "\\."), function(x) x[1]))
tissues <- c("OB", "PFC", "HIPPO")
filesList <- list(OB = list(countfiles = c(files[grep("wt_OB", files)], files[grep("hemi_OB", files)]), 
                            sampleData = data.frame(condition = c(rep("wt", 3), rep("hemi", 3)), 
                                                    row.names = names(c(files[grep("wt_OB", files)], files[grep("hemi_OB", files)])))),
                  PFC = list(countfiles = c(files[grep("wt_PFC", files)], files[grep("hemi_PFC", files)]),
                             sampleData = data.frame(condition = c(rep("wt", 3), rep("hemi", 3)),
                                                     row.names = names(c(files[grep("wt_PFC", files)], files[grep("hemi_PFC", files)])))),
                  HIPPO = list(countfiles = c(files[grep("wt_HIPPO", files)], files[grep("hemi_HIPPO", files)]),
                               sampleData = data.frame(condition = c(rep("wt", 3), rep("hemi", 3)),
                                                       row.names = names(c(files[grep("wt_PFC", files)], files[grep("hemi_PFC", files)])))))


if (!file.exists("DEXSeqDataSetList.rda")){
  DEXSeqDataSetList <- lapply(names(filesList), function(x){
    DEXSeqDataSet <- DEXSeq::DEXSeqDataSetFromHTSeq(countfiles = as.vector(filesList[[x]][["countfiles"]]),
                                                    sampleData = filesList[[x]][["sampleData"]],
                                                    design = ~ sample + exon + condition:exon,
                                                    flattenedfile = flattenedfile)
    DEXSeqDataSet <- DEXSeq::estimateSizeFactors(DEXSeqDataSet)
    DEXSeqDataSet <- DEXSeq::estimateDispersions(DEXSeqDataSet, BPPARAM = BPPARAM)
    DEXSeqDataSet <- DEXSeq::fitDispersionFunction(DEXSeqDataSet)
    DEXSeqDataSet <- DEXSeq::testForDEU(DEXSeqDataSet, BPPARAM = BPPARAM)
    DEXSeqDataSet <- DEXSeq::estimatelog2FoldChanges(DEXSeqDataSet)
    results <- DEXSeq::DEUresultTable(DEXSeqDataSet)
    return(list(DEXSeqDataSet = DEXSeqDataSet,
                results = results))
  })
  names(DEXSeqDataSetList) <- names(filesList)
  save(DEXSeqDataSetList, file = "DEXSeqDataSetList.rda")
} else {
  load("DEXSeqDataSetList.rda")
}

