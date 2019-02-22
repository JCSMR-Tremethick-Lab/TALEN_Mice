library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)

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


# inspect peak calling data -----------------------------------------------
macs2broad <- "/home/sebastian/Data/Tremethick/TALENs/CutRun/macs2/callpeak/broad/GRCm38_ensembl93/NB501086_0221_TSoboleva_JCSMR_CutandRun/"

broadPeakFiles <- list.files(macs2broad, recursive = T, pattern = "broadPeak", full.names = T)
names(broadPeakFiles) <- list.files(macs2broad, recursive = T, pattern = "broadPeak", full.names = F)
names(broadPeakFiles) <- gsub(".broadPeak", "", unlist(lapply(strsplit(names(broadPeakFiles), "/"), function(x) x[2])))
broadPeakFiles <- broadPeakFiles[grep("35sec", broadPeakFiles, invert = T)]

broadPeaks <- lapply(broadPeakFiles, function(x) data.table::fread(x))
broadPeaksGRL <- lapply(broadPeaks, function(x) {
  GenomicRanges::makeGRangesFromDataFrame(x,
                                          keep.extra.columns = T, 
                                          seqnames.field = "V1", 
                                          start.field = "V2", 
                                          end.field = "V3")
})
broadPeaksGRL <- GenomicRanges::GenomicRangesList(broadPeaksGRL)
names(broadPeaksGRL) <- names(broadPeakFiles)
sapply(broadPeaksGRL, function(x) summary(width(x)))
