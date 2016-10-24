# STAR/HTSeq analysis
# within tissue analysis of mouse brain RNA-seq data
require(EDASeq)
require(biomaRt)
require(tidyr)
require(rtracklayer)
require(RColorBrewer)
require(RUVSeq)
require(NOISeq)
require(deepToolsUtils)

# load external functions
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

# global variables
ensemblHost <- "mar2016.archive.ensembl.org"
colors <- RColorBrewer::brewer.pal(3, "Set2")

# setting working directory and data sources ------------------------------
if (amILocal("JCSMR027564ML")){
  setwd("~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/HTSeq/count/"
  pathPrefix <- "~/mount/gduserv/"
} else {
  setwd("~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/HTSeq/count/"
  pathPrefix <-  "~"
}

files <- list.files(path = dataPath, full.names = T)
names(files) <- list.files(path = dataPath, full.names = F)

if (!file.exists("htSeqCountMatrix.rda")){
  htSeqCountMatrix <- deepToolsUtils::makeHTSeqCountMatrix(files)
  # re-order columns into corresponding sample/condition groups
  reOrdCol <- c(grep("wt_HIPPO", colnames(htSeqCountMatrix)),
                grep("hemi_HIPPO", colnames(htSeqCountMatrix)),
                grep("wt_PFC", colnames(htSeqCountMatrix)),
                grep("hemi_PFC", colnames(htSeqCountMatrix)),
                grep("wt_OB", colnames(htSeqCountMatrix)),
                grep("hemi_OB", colnames(htSeqCountMatrix)))
  htSeqCountMatrix <- htSeqCountMatrix[, reOrdCol]
  # remove the HTSeq count summary stats
  startn <- nrow(htSeqCountMatrix)-4
  endn <- nrow(htSeqCountMatrix)
  htSeqCountReport <- htSeqCountMatrix[startn:endn, ]
  htSeqCountMatrix <- htSeqCountMatrix[-(startn:endn), ]
  save(htSeqCountMatrix, file = "htSeqCountMatrix.rda")
} else {
  load("htSeqCountMatrix.rda")
}

# quick summary of the count data
totalCounts <- apply(htSeqCountMatrix, 2, sum)
ERCCCounts <- apply(htSeqCountMatrix[grep("ERCC", rownames(htSeqCountMatrix)), ], 2, sum)
ERCCPerc <- ERCCCounts / totalCounts * 100

# preparing annotation data from Ensembl ----------------------------------
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

ensTranscripts <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                "ensembl_transcript_id",
                                                "transcript_length"),
                                 mart = mouse,
                                 filter = "ensembl_gene_id",
                                 values = ensGenes$ensembl_gene_id)
save(ensTranscripts, file = "ensTranscripts.rda")

mylength <- sapply(ensGenes$ensembl_gene_id, function(x){
  y <- ensTranscripts[which(ensTranscripts$ensembl_gene_id == x), ]
  y <- y[which.max(y$transcript_length), ]$transcript_length
})
save(mylength, file = "mylength.rda")

mygc <- ensGenes$percentage_gc_content
names(mygc) <- ensGenes$ensembl_gene_id
save(mygc, file = "mygc.rda")

mybiotypes <- ensGenes$gene_biotype
names(mybiotypes) <- ensGenes$ensembl_gene_id
save(mybiotypes, file = "mybiotypes.rda")

mychroms <- data.frame(Chr = ensGenes$chromosome_name, GeneStart = ensGenes$start_position, GeneEnd = ensGenes$end_position)
rownames(mychroms) <- ensGenes$ensembl_gene_id
save(mychroms, file = "mychroms.rda")
} else {
  load("ensGenes.rda")
  load("ensTranscripts.rda")
  load("mylength.rda")
  load("mygc.rda")
  load("mybiotypes.rda")
  load("mychroms.rda")
}


# create separate table for each brain region  ----------------------------
htSeqCountMatrix.ob <- htSeqCountMatrix[, grep("OB", colnames(htSeqCountMatrix))]
htSeqCountMatrix.pfc <- htSeqCountMatrix[, grep("PFC", colnames(htSeqCountMatrix))]
htSeqCountMatrix.hippo <- htSeqCountMatrix[, grep("HIPPO", colnames(htSeqCountMatrix))]
htSeqCountMatrix.brain <- list(list(htSeqCountMatrix.ob), list(htSeqCountMatrix.pfc), list(htSeqCountMatrix.hippo))
names(htSeqCountMatrix.brain) <- c("OB", "PFC", "HIPPO")

analysis_version <- 2
analysis_output_file <- paste("DifferentialGeneExpressionAnalysis_", analysis_version, ".rda", sep = "")


# actual processing -------------------------------------------------------
if (!file.exists(analysis_output_file)){
  processedData <- lapply(names(htSeqCountMatrix.brain), function(x){
    original <- htSeqCountMatrix.brain[[x]][[1]]
    filter <- apply(original, 1, function(y) length(y[y>5])>=2)
    filtered <- original[filter, ]
    genes <- rownames(filtered)[grep("ENS", rownames(filtered))]
    spikes <- rownames(filtered)[grep("ERCC", rownames(filtered))]
    #---------------------------------------------
    # create expression set
    condition <- as.factor(unlist(lapply(strsplit(colnames(filtered), "_"), function(z) z[3])))
    condition <- factor(condition, levels(condition)[c(2,1)])
    set <- newSeqExpressionSet(as.matrix(filtered),
                               phenoData = data.frame(condition, row.names=colnames(filtered)))
    #---------------------------------------------
    # data exploration
    rle <- plotRLE(set, outline = FALSE, col = colors[condition])
    pca <- plotPCA(set, col = colors[condition])
    #---------------------------------------------
    # NOISeq analysis
    # create eSet
    tissue <- as.factor(unlist(lapply(strsplit(colnames(original), "_"), function(z) z[4])))
    myfactors <- data.frame(Tissue = tissue, Condition = condition)
    mydata <- NOISeq::readData(data = original, 
                               length = mylength, 
                               gc = mygc, 
                               biotype = mybiotypes, 
                               chromosome = mychroms, 
                               factors = myfactors)
    # actual data exploration
    myexplodata <- dat(mydata, type = "biodetection")
    mybiodetection <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
    mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
    mysaturation <- dat(mydata, k = 0, ndepth = 10, type = "saturation")
    # checking biases
    # length
    mylengthbias <- dat(mydata, factor = "Tissue", type = "lengthbias")
    myGCbias <- dat(mydata, factor = "Tissue", type = "GCbias")
    # compositional bias
    mycd <- dat(mydata, type = "cd", norm = FALSE, refColumn = 1)
    #---------------------------------------------
    # upper quartile normalization
    set <- betweenLaneNormalization(set, which="upper")
    rleUQ <- plotRLE(set, outline = FALSE, col = colors[condition])
    pcaUQ <- plotPCA(set, col = colors[condition])
    #---------------------------------------------
    # RUVg using spike ins
    set1 <- RUVg(set, spikes, k = 1)
    rleRUVg <- plotRLE(set1, outline=FALSE, col=colors[condition])
    pcaRUVg <- plotPCA(set1, col=colors[condition], cex=1.2)
    #---------------------------------------------
    # NOISeq post-normalization & RUV
    mydataNorm <- NOISeq::readData(data = normCounts(set1), 
                               length = mylength, 
                               gc = mygc, 
                               biotype = mybiotypes, 
                               chromosome = mychroms, 
                               factors = myfactors)
    # actual data exploration
    myexplodataNorm <- dat(mydataNorm, type = "biodetection", norm = TRUE)
    mybiodetectionNorm <- dat(mydataNorm, k = 0, type = "biodetection", factor = NULL, norm = TRUE)
    mycountsbioNorm <- dat(mydataNorm, factor = NULL, type = "countsbio",  norm = TRUE)
    mysaturationNorm <- dat(mydataNorm, k = 0, ndepth = 10, type = "saturation", norm = TRUE)
    # checking biases
    # length
    mylengthbiasNorm <- dat(mydataNorm, factor = "Tissue", type = "lengthbias", norm = TRUE)
    myGCbiasNorm <- dat(mydataNorm, factor = "Tissue", type = "GCbias", norm = TRUE)
    # compositional bias
    mycdNorm <- dat(mydataNorm, type = "cd", refColumn = 1, norm = TRUE)
    #---------------------------------------------
    # differential expression analysis using edgeR
    # here including the RUVg factor to account for "unwanted variation"
    design <- model.matrix(~condition + W_1, data = pData(set1))
    y <- DGEList(counts = counts(set1), group = condition)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2)
    tt <- topTags(lrt, n = 5000)
    annotatedTT <- merge(tt[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id")
    #---------------------------------------------
    # differential expression analysis using plain vanilla edgeR
    design1 <- model.matrix(~condition, data = pData(set))
    y1 <- DGEList(counts = counts(set), group = condition)
    y1 <- calcNormFactors(y1, method="upperquartile")
    y1 <- estimateGLMCommonDisp(y1, design1)
    y1 <- estimateGLMTagwiseDisp(y1, design1)
    fit1 <- glmFit(y1, design1)
    lrt1 <- glmLRT(fit1, coef=2)
    tt1 <- topTags(lrt1, n = 5000)
    annotatedTT1 <- merge(tt1[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, sort = F)
    #---------------------------------------------
    # return all objects
    return(list(original = original, 
                filtered = filtered,
                genes = genes,
                spikes = spikes,
                eSet = set,
                RLEplot = rle,
                PCAplot = pca,
                RLEplotUQ = rleUQ,
                RLEplotPCA = pcaUQ,
                eSetRUVg = set1,
                RLEplotRUVg = rleRUVg,
                PCAplotRUVg = pcaRUVg,
                DGEList = y,
                glmFit = fit,
                glmLRT = lrt,
                topTags = tt,
                AnnotatedTopTags = annotatedTT,
                DGEList_noRUV = y1,
                glmFit_noRUV = fit1,
                glmLRT_noRUV = lrt1,
                topTags_noRUV = tt1,
                AnnotatedTopTags_noRUV = annotatedTT1,
                NOISeqRaw = list(mydata = mydata,
                                 myexplodata = myexplodata,
                                 mybiodetection = mybiodetection,
                                 mycountsbio = mycountsbio,
                                 mysaturation = mysaturation,
                                 mylengthbias = mylengthbias,
                                 myGCbias = myGCbias,
                                 mycd = mycd),
                NOISeqNorm = list(mydataNorm = mydataNorm,
                                  myexplodataNorm = myexplodataNorm,
                                  mybiodetectionNorm = mybiodetectionNorm,
                                  mycountsbioNorm = mycountsbioNorm,
                                  mysaturationNorm = mysaturationNorm,
                                  mylengthbiasNorm = mylengthbiasNorm,
                                  myGCbiasNorm = myGCbiasNorm,
                                  mycdNorm = mycdNorm)))
  })
  names(processedData) <- names(htSeqCountMatrix.brain)
  save(processedData, file = analysis_output_file)
} else {
  load(analysis_output_file)
}

# write top table to CSV
sapply(names(processedData), function(x){
  write.csv(processedData[[x]][["AnnotatedTopTags"]][order(processedData[[x]][["AnnotatedTopTags"]]$FDR),], file = paste("TALEN_mouse_", x, "_Differential_Gene_Expression.csv", sep = ""))
})

sapply(names(processedData), function(x){
  table(processedData[[x]][["AnnotatedTopTags"]]$FDR <= 0.1)
})


sapply(names(processedData), function(x){
  table(processedData[[x]][["AnnotatedTopTags_noRUV"]]$FDR <= 0.1)
})

# Plot PCAs pre- and post-RUVg for each group
sapply(names(processedData), function(x){
  outfile <- paste(x, "/RUVg_pre_post_PCA.pdf", sep = "")
  pdf(outfile, paper = "a4r")
  par(mfrow = c(2,1))
  plotPCA(processedData[[x]][["eSet"]],
          col=colors[pData(processedData[[x]][["eSet"]])$condition],
          cex=1.2, 
          main = paste(x, " pre RUVg processing", sep = ""))
  plotPCA(processedData[[x]][["eSetRUVg"]], 
          col=colors[pData(processedData[[x]][["eSetRUVg"]])$condition], 
          cex=1.2, 
          main = paste(x, " post RUVg processing", sep = ""))
  dev.off()
})


# EDA exploration ---------------------------------------------------------
# only for one tissue at the moment
common <- intersect(names(mygc), rownames(processedData[["HIPPO"]][["filtered"]]))
data <- newSeqExpressionSet(counts = as.matrix(processedData[["HIPPO"]][["filtered"]][common,]),
                            featureData = data.frame(gc = mygc[common], length = mylength[common]))
pData(data) <- pData(processedData[["HIPPO"]][["eSet"]])

boxplot(data,col=colors[1:3])
MDPlot(data,c(1,2))
par(mfrow = c(1,2))
meanVarPlot(data[,1:3], log = T)
meanVarPlot(data[,4:6], log = T)
meanVarPlot(data, log = T)
biasPlot(data, "gc", log = T)
biasPlot(data, "length", log = T)

lfc <- log(counts(data)[,2]+0.1) - log(counts(data)[,5]+0.1)
biasBoxplot(lfc, fData(data)$gc)
biasBoxplot(lfc, fData(data)$length)

dataWithin <- withinLaneNormalization(data,"gc", which="upper")
dataNorm <- betweenLaneNormalization(dataWithin, which="full")

biasPlot(dataWithin, "gc", log = T)

dataOffset <- withinLaneNormalization(data,"gc", which="full",offset=TRUE)
dataOffset <- betweenLaneNormalization(dataOffset, which="full",offset=TRUE)
dataOffset <- 

