library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(biomaRt)
library(BS)

# set working & data directory --------------------------------------------
setwd("/home/sebastian/Data/Tremethick/TALENs/R_analysis")
dataDir <- "/home/sebastian/Data/Tremethick/TALENs/CutRun/deepTools/bamCompare/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"


# import coverage data ----------------------------------------------------
covData <- rtracklayer::import(file.path(dataDir, "KO_01_K36me3_readCount.bw"), format = "BigWig")
covData <- sortSeqlevels(covData)
ensGenes$strand <- c("+", "-")[match(ensGenes$strand, c(1, -1))]
grEnsGenes <- GenomicRanges::makeGRangesFromDataFrame(ensGenes, 
                                                      start.field = "start_position", 
                                                      end.field = "end_position", 
                                                      keep.extra.columns = T)

grEnsGenes <- sortSeqlevels(grEnsGenes)

commonLevels <- intersect(seqlevels(grEnsGenes), seqlevels(covData))
grEnsGenes <- keepSeqlevels(grEnsGenes, commonLevels, pruning.mode = "coarse")
covData <- keepSeqlevels(covData, commonLevels, pruning.mode = "coarse")
seqinfo(grEnsGenes) <- seqinfo(covData)

subsetByOverlaps(covData, grEnsGenes)
