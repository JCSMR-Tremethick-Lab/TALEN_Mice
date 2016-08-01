require("DEXSeq")
require(biomaRt)
require(parallel)
require(rtracklayer)
#
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

# defining global variables
ensemblHost <- "mar2016.archive.ensembl.org"

if (amILocal("JCSMR027564ML")){
  mount <- system("mount", intern = T)
  if (length(grep("gduserv", mount)) == 0) {system("sshfs skurscheid@gduserv.anu.edu.au: ~/mount/gduserv/")}
  BPPARAM <- MulticoreParam(workers = 7)
  setwd("~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/DEXSeq/count/"
  # biomaRt connection
  mouse <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
  attribs <- biomaRt::listAttributes(mouse)
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

if (!file.exists("ensGenes.rda")){
  mouse <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
  attribs <- biomaRt::listAttributes(mouse)
  ensGenes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                            "external_gene_name",
                                            "chromosome_name",
                                            "start_position",
                                            "end_position",
                                            "strand",
                                            "band",
                                            "description",
                                            "percentage_gc_content",
                                            "gene_biotype"),
                             mart = mouse)
  save(ensGenes, file = "ensGenes.rda")
  } else {
    load("ensGenes.rda")
    rownames(ensGenes) <- ensGenes$ensembl_gene_id
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
    DEXSeqDataSet <- DEXSeq::testForDEU(DEXSeqDataSet, BPPARAM = BPPARAM)
    DEXSeqDataSet <- DEXSeq::estimateExonFoldChanges(DEXSeqDataSet, BPPARAM = BPPARAM)
    results <- DEXSeq::DEXSeqResults(DEXSeqDataSet)
    return(list(DEXSeqDataSet = DEXSeqDataSet,
                results = results))
  })
  names(DEXSeqDataSetList) <- names(filesList)
  save(DEXSeqDataSetList, file = "DEXSeqDataSetList.rda")
} else {
  load("DEXSeqDataSetList.rda")
}

# prepare data.frame with additinal gene annotation to include in report
GTF <- import(flattenedfile)
names(GTF) <- paste(GTF$gene_id, GTF$exonic_part_number, sep = ":")
dfGTF <- data.frame(type = GTF$type, ensembl_gene_id = GTF$gene_id, row.names = names(GTF))
dfGTF <- dfGTF[dfGTF$type == "exonic_part",]
dfGTF <- merge(dfGTF, ensGenes, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T, sort = F)

lapply(names(DEXSeqDataSetList), function(x){
  if (table(DEXSeqDataSetList[[x]][["results"]]$padj < 0.1)[2] < 21){
    DEXSeqHTML(DEXSeqDataSetList[[x]][["results"]], 
               FDR = 0.1, 
               color = c("red", "blue"),
               path = x,
               extraCols = ensGenes)
  } else {
    DEXSeqHTML(DEXSeqDataSetList[[x]][["results"]], 
               FDR = 0.1, 
               color = c("red", "blue"),
               path = x,
               BPPARAM = BPPARAM) 
  }
})


# some data exploration
res <- DEXSeqDataSetList[["OB"]][["results"]]
res <- res[order(res$padj),]
table(res$padj < 0.1)
plotDEXSeq(res, "ENSMUSG00000027361", legend = T, cex.axis = 1.2, cex = 1.3, lwd = 2, expression = F, norCounts = T, displayTranscripts = T)
