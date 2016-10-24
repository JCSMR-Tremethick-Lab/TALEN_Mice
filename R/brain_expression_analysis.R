require(rhdf5)
require(sleuth)
require(biomaRt)
require(tidyr)
require(rtracklayer)
require(BiocParallel)
require(tximport)
require(readr)
require(RUVSeq)

source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

# defining global variables
ensemblHost <- "mar2016.archive.ensembl.org"
colors <- RColorBrewer::brewer.pal(3, "Set2")

# setting working directory and data sources ------------------------------
if (amILocal("JCSMR027564ML")){
  mount <- system("mount", intern = T)
  if (length(grep("gduserv", mount)) == 0) {system("sshfs skurscheid@gduserv.anu.edu.au: ~/mount/gduserv/")}
  BPPARAM <- BiocParallel::MulticoreParam(workers = 7)
  mc.cores <- 8L
  pathPrefix <- "~/mount/gduserv/"
  # biomaRt connection
  mouse <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
  attribs <- biomaRt::listAttributes(mouse)
} else {
  BPPARAM <- BiocParallel::MulticoreParam(workers = 16)
  mc.cores <- 16L
  pathPrefix <- "~"
  options(width = 137)
}

lDir <- function(x, y){
  paste(x, y, sep = "/")
}

setwd(lDir(pathPrefix, "Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis"))
dataPathDEXSeq <- lDir(pathPrefix, "Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/DEXSeq/count/")
dataPathHTSeq <- lDir(pathPrefix, "Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/HTSeq/count")
flattenedfile <- lDir(pathPrefix, "Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.DEXSeq.gtf")
kallisto_base_dir <- lDir(pathPrefix, "Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto")


# use Ensembl 84 for annotation
if (!file.exists("t2g.rda")){
  mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
  attribs <- listAttributes(mouse)
  # annotate transcripts
  # Ensembl 84 includes version number in FASTA IDs, therefore have to conactenate them in, other wise mapping does not work
  t2g <- getBM(attributes = c("ensembl_transcript_id", 
                              "ensembl_gene_id", 
                              "external_gene_name", 
                              "version", 
                              "transcript_version"), 
               mart = mouse)
  t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  # add ERCC spike ins
  ercc <- import(paste(pathPrefix, "/Data/References/Transcriptomes/ERCC/ERCC92.gtf", sep = ""))
  ercc.df <- mcols(ercc)
  ercc.df <- data.frame(ercc.df[, c("transcript_id", "gene_id", "gene_id")])
  colnames(ercc.df) <- c("target_id", "ens_gene", "ext_gene")
  ercc.df$version <- 1
  ercc.df$transcript_version <- 1
  ercc.df$target_id <- ercc.df$ens_gene
  t2g <- rbind(t2g, ercc.df)
  rownames(t2g) <- t2g$target_id
  save(t2g, file = "t2g.rda")
} else {
  load("t2g.rda")
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

# transcript level differential expression analysis -----------------------------
# run analysis by brain region so that only pair-wise comparison

base_dirs <- c(OB = paste(pathPrefix, "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_ob", sep = ""),
               PFC = paste(pathPrefix, "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_pfc", sep = ""),
               HIPPO = paste(pathPrefix, "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_hippo", sep = ""))

sleuth_analysis_version <- 4
sleuth_analysis_output <- paste("sleuth_analysis_V", sleuth_analysis_version, "_output.rda", sep = "")
sleuth_analysis_output_combined <- paste("sleuth_analysis_V", sleuth_analysis_version, "_output_combined.rda", sep = "")
sleuth_analysis_outputCompressed <- paste("sleuth_analysis_V", sleuth_analysis_version, "_Compressed_output.rda", sep = "")

previous_sleuth_analysis_version <- sleuth_analysis_version - 1
previous_sleuth_analysis_output <- paste("sleuth_analysis_V", previous_sleuth_analysis_version, "_output.rda", sep = "")
previous_sleuth_analysis_output_combined <- paste("sleuth_analysis_V", previous_sleuth_analysis_version, "_output_combined.rda", sep = "")

if (file.exists(previous_sleuth_analysis_output)) {file.remove(previous_sleuth_analysis_output); file.remove(previous_sleuth_analysis_output_combined)}
  
if (!file.exists(sleuth_analysis_output)){
  sleuthProcessedData <- lapply(names(base_dirs), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    options(mc.cores = mc.cores)
    sample_id <- dir(base_dirs[[x]])
    kal_dirs <- sapply(sample_id, function(id) file.path(base_dirs[[x]], id))
    condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[3], collapse = "_")))
    condition <- as.factor(condition)
    condition <- factor(condition, levels(condition)[c(2,1)])
    s2c <- data.frame(sample = sample_id, condition = condition)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    # transcript level
    so <- sleuth::sleuth_prep(s2c, ~ condition, target_mapping = t2g)
    so <- sleuth::sleuth_fit(so)
    so <- sleuth::sleuth_wt(so, "conditionhemi")
    so <- sleuth::sleuth_fit(so, ~1, "reduced")
    so <- sleuth::sleuth_lrt(so, "reduced", "full")
    kt <- sleuth::kallisto_table(so)
    rt <- sleuth::sleuth_results(so, "conditionhemi")
    rt <- rt[order(rt$qval),]
    kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt_wide) <- kt_wide$target_id
    kt_wide <- kt_wide[,-1]
    kt_wide <- kt_wide[, c(grep("wt", colnames(kt_wide)),
                           grep("hemi", colnames(kt_wide)))]
    kt.pca <- ade4::dudi.pca(t(as.matrix(kt_wide)), center = T, scale = T, scannf = F, nf = 3)
    return(list(sleuth_object = so,
                sleuth_results = rt,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                kallisto_pca = kt.pca))
    })
  names(sleuthProcessedData) <- names(base_dirs)
  save(sleuthProcessedData, file = sleuth_analysis_output)
} else {
  load(sleuth_analysis_output)
}

if (!file.exists(sleuth_analysis_outputCompressed)){
  sleuthProcessedDataCompressed <- lapply(names(sleuthProcessedData), function(x){
    sleuthProcessedData[[x]][grep("sleuth_object", names(sleuthProcessedData[[x]]), invert = T)]
  })
  names(sleuthProcessedDataCompressed) <- names(sleuthProcessedData)
  save(sleuthProcessedDataCompressed, file = sleuth_analysis_outputCompressed)
} else {
  load(sleuth_analysis_outputCompressed)
}

# Use tximport package to create gene level counts from kallisto --------
# then proceed to run RUVseq followed by edgeR

edgeR_analysis_version <- 5
edgeR_analysis_output <- paste("edgeR_analysis_V", edgeR_analysis_version, "_output.rda", sep = "")

# previous_edgeR_analysis_version <- edgeR_analysis_version - 1
# previous_edgeR_analysis_output <- paste("edgeR_analysis_V", previous_edgeR_analysis_version, "_output.rda", sep = "")
# if (file.exists(previous_edgeR_analysis_output)) {file.remove(previous_edgeR_analysis_output)}

# actual edgeR analysis
if (!file.exists(edgeR_analysis_output)){
  edgeRProcessedData <- lapply(names(base_dirs), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    options(mc.cores = mc.cores)
    sample_id <- dir(base_dirs[[x]])
    kal_dirs <- sapply(sample_id, function(id) file.path(base_dirs[[x]], id))
    condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[3], collapse = "_")))
    condition <- as.factor(condition)
    condition <- factor(condition, levels(condition)[c(2,1)])
    s2c <- data.frame(sample = sample_id, condition = condition)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    s2c$files <- paste(s2c$path, "abundance.tsv", sep = "/")
    s2c <- s2c[order(s2c$condition),]
    files <- s2c$files
    names(files) <- s2c$sample
    # read in kallisto data with tximport
    txi <- tximport::tximport(files,
                              type = "kallisto",
                              tx2gene = t2g,
                              geneIdCol = "ens_gene",
                              txIdCol = "target_id",
                              reader = read_tsv)
    # RUVseq analysis starts here
    # have to round the estimated counts to integers
    original <- round(txi$counts, 0)
    filter <- apply(original, 1, function(y) length(y[y>5])>=2)
    filtered <- original[filter, ]
    genes <- rownames(filtered)[grep("ENS", rownames(filtered))]
    spikes <- rownames(filtered)[grep("ERCC", rownames(filtered))]
    #---------------------------------------------
    # create expression set
    if (x == "combined"){
      set <- EDASeq::newSeqExpressionSet(as.matrix(filtered),
                                 phenoData = data.frame(condition, tissue, individual, row.names=colnames(filtered)))
    } else {
      set <- EDASeq::newSeqExpressionSet(as.matrix(filtered),
                                         phenoData = data.frame(condition, row.names=colnames(filtered)))
    }
    #---------------------------------------------
    # data exploration
    rle <- EDASeq::plotRLE(set, outline = FALSE, col = colors[s2c$condition])
    pca <- EDASeq::plotPCA(set, col = colors[s2c$condition])
    #---------------------------------------------
    # RUVseq using spike ins
    set1 <- RUVSeq::RUVg(set, spikes, k = 1)
    rleRUV <- EDASeq::plotRLE(set1, outline=FALSE, col=colors[s2c$condition])
    pcaRUV <- EDASeq::plotPCA(set1, col=colors[s2c$condition], cex=1.2)
    #---------------------------------------------
    # differential expression analysis using edgeR
    # here including the RUVg factor to account for "unwanted variation"
    cts <- txi$counts[filter,]
    normMat <- txi$length[filter,]
    normMat <- normMat/exp(rowMeans(log(normMat)))
    design <- model.matrix(~ condition + W_1, data = pData(set1))
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
    y <- edgeR::DGEList(counts = counts(set1), group = s2c$condition)
    y$offset <- t(t(log(normMat)) + o)
    y <- edgeR::estimateGLMCommonDisp(y, design)
    y <- edgeR::estimateGLMTagwiseDisp(y, design)
    fit <- edgeR::glmFit(y, design)
    lrt <- edgeR::glmLRT(fit, coef=2)
    tt <- edgeR::topTags(lrt, n = 20000)
    annotatedTT <- merge(tt[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id")
    annotatedTT <- annotatedTT[order(annotatedTT$FDR),]
    #---------------------------------------------
    # differential expression analysis using plain vanilla edgeR
    design1 <- model.matrix(~condition, data = pData(set))
    y1 <- edgeR::DGEList(counts = counts(set), group = condition)
    y1 <- edgeR::calcNormFactors(y1, method="upperquartile")
    y1 <- edgeR::estimateGLMCommonDisp(y1, design1)
    y1 <- edgeR::estimateGLMTagwiseDisp(y1, design1)
    fit1 <- edgeR::glmFit(y1, design1)
    lrt1 <- edgeR::glmLRT(fit1, coef=1)
    tt1 <- edgeR::topTags(lrt1, n = 5000)
    annotatedTT1 <- merge(tt1[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, sort = F)
    annotatedTT1 <- annotatedTT1[order(annotatedTT1$FDR), ]
    # return all objects
    return(list(txiGeneCounts = txi,
                condition = condition,
                preRUVSet = set,
                preRUVRle = rle,
                preRUVPca = pca,
                postRUVSet = set1,
                postRUVRle = rleRUV,
                postRUVPca = pcaRUV,
                DGEList = y,
                glmFit = fit,
                glmLRT = lrt,
                topTags = tt,
                AnnotatedTopTags = annotatedTT,
                DGEList_noRUV = y1,
                glmFit_noRUV = fit1,
                glmLRT_noRUV = lrt1,
                topTags_noRUV = tt1,
                AnnotatedTopTags_noRUV = annotatedTT1))
  })
  names(edgeRProcessedData) <- names(base_dirs)
  save(edgeRProcessedData, file = edgeR_analysis_output)
  } else {
  load(edgeR_analysis_output)
}

# see how many genes are DE @ 10%FDR
sapply(edgeRProcessedData, function(x){
  table(x[["AnnotatedTopTags"]]$FDR < 0.1)
})

sapply(edgeRProcessedData, function(x){
  table(x[["AnnotatedTopTags_noRUV"]]$FDR < 0.1)
})

# get txi scaled TPMs for each group
txiGeneCountsScaled <- lapply(names(base_dirs), function(x){
  sample_id <- dir(base_dirs[[x]])
  kal_dirs <- sapply(sample_id, function(id) file.path(base_dirs[[x]], id))
  condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[3], collapse = "_")))
  condition <- as.factor(condition)
  condition <- factor(condition, levels(condition)[c(2,1)])
  s2c <- data.frame(sample = sample_id, condition = condition)
  s2c <- dplyr::mutate(s2c, path = kal_dirs)
  s2c$files <- paste(s2c$path, "abundance.tsv", sep = "/")
  s2c <- s2c[order(s2c$condition),]
  files <- s2c$files
  names(files) <- s2c$sample
  # read in kallisto data with tximport
  txi <- tximport::tximport(files,
                            type = "kallisto",
                            tx2gene = t2g,
                            geneIdCol = "ens_gene",
                            txIdCol = "target_id",
                            reader = read_tsv,
                  
                                      countsFromAbundance = "lengthScaledTPM")
})
names(txiGeneCountsScaled) <- names(edgeRProcessedData)

# run PCA on gene-level abundances
geneLevelPCA <- lapply(txiGeneCountsScaled, function(x) {
  ade4::dudi.pca(t(x[["abundance"]]), scannf = F, nf = 5)
})

# plot s.class diagrams
pdf("Individual_PCAs_gene_level.pdf")
par(mfrow = c(3,1))
lapply(names(geneLevelPCA), function(x){
  ade4::s.class(geneLevelPCA[[x]]$li, 
                fac = as.factor(unlist(lapply(strsplit(rownames(geneLevelPCA[[x]]$tab), "_"), function(y) y[3]))),
                sub = x)
})
dev.off()

masterDFGeneLevel <- cbind(txiGeneCountsScaled[["OB"]][["abundance"]], 
                  txiGeneCountsScaled[["PFC"]][["abundance"]],
                  txiGeneCountsScaled[["HIPPO"]][["abundance"]])
masterPCAGeneLevel <- ade4::dudi.pca(t(masterDF), scannf = F, nf = 5)

pdf("Combined_PCA_gene_level.pdf")
ade4::s.class(masterPCAGeneLevel$li, 
              fac = as.factor(unlist(lapply(strsplit(rownames(masterPCAGeneLevel$tab), "_"), function(y) paste(y[3], y[4], collapse = "_")))))
dev.off()

# plot s.class diagrams on tx-level PCA data
pdf("Individual_PCAs_transcript_level.pdf")
par(mfrow = c(3,1))
lapply(names(sleuthProcessedData), function(x){
  ade4::s.class(sleuthProcessedData[[x]][["kallisto_pca"]]$li, 
                fac = as.factor(unlist(lapply(strsplit(rownames(sleuthProcessedData[[x]][["kallisto_pca"]]$tab), "_"), function(y) y[3]))),
                sub = x)
})
dev.off()

masterDFTxLevel <- cbind(sleuthProcessedData[["OB"]][["kallisto_table_wide"]], 
                         sleuthProcessedData[["PFC"]][["kallisto_table_wide"]],
                         sleuthProcessedData[["HIPPO"]][["kallisto_table_wide"]])
masterPCATxLevel <- ade4::dudi.pca(t(masterDFTxLevel), scannf = F, nf = 5)

pdf("Combined_PCA_transcript_level.pdf")
ade4::s.class(masterPCATxLevel$li, 
              fac = as.factor(unlist(lapply(strsplit(rownames(masterPCATxLevel$tab), "_"), function(y) paste(y[3], y[4], collapse = "_")))))
dev.off()

# get transcript level differential expression results for each group

sleuthResultsList <- lapply(sleuthProcessedData, function(x){
  sleuthResults <- sleuth_results(x[["sleuth_object"]], test = "conditionhemi", test_type = "wald")
})

sapply(sleuthResultsList, function(x){
  table(x$qval < 0.1)
})

#-------------------------------------------------------------------------------
# look at what proportions of reads where used up by ERCC spike-ins
totalCounts <- apply(kt.raw_counts, 2, sum)
ERCCCounts <- apply(kt.raw_counts[grep("ERCC", rownames(kt.raw_counts)), ], 2, sum)
ERCCPerc <- ERCCCounts / totalCounts * 100

kt <- kallisto_table(so)

load("~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/kt.rda")
kt <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
rownames(kt) <- kt$target_id

# check ERCC expression ---------------------------------------------------
ercc_tab <- read.table("~/mount/gduserv/Data/References/Transcriptomes/ERCC/cms_095046.txt", header = T, as.is = T, sep = "\t")
rownames(ercc_tab) <- ercc_tab$ERCC.ID
ercc_exp <- kt[grep("ERCC", kt$target_id),]
i1 <- intersect(rownames(ercc_exp[which(ercc_exp$mean > 0),]), rownames(ercc_tab))

pdf("ERCC_spike_ins.pdf", paper = "a4")
par(mfrow = c(3,2))
lapply(seq_along(1:18), function(x){
  x <- x+1
  s <- colnames(ercc_exp)[x]
  y <- ercc_exp[i1, x]
  names(y) <- i1
  y1 <- which(y > 0.01)
  y <- y[names(y1)]
  e <- unlist(ercc_tab[i1, "concentration.in.Mix.1..attomoles.ul."])
  names(e) <- i1
  e <- e[names(y1)]
  cor <- cor(log2(e), log2(y)) 
  sm <- summary(lm(log2(y) ~ log2(e)))
  plot(log2(e), log2(y),
       main = paste(s, " [#", length(y1), " ERCCs]", sep = ""),
       xlab = "Log2(aM) Mix 1",
       ylab = "Log2(tpm)")
  lines(lowess(log2(y) ~ log2(e)), col="green3")
  abline(lm(log2(y) ~ log2(e)),col="red3")
  legend("bottom", legend = c(paste("cor = ", round(cor, 2), sep = ""), paste("R2 = ", round(sm$adj.r.squared, 2), sep = "")))
})
dev.off()

# perform PCA
kt.geneLog2 <- log2(kt[,c(2:19)] + 1)
pca1 <- ade4::dudi.pca(t(kt.geneLog2))
samples <- s2c$condition
pdf("PCA.pdf", height = 7, width = 7)
s.class(pca1$li, fac = samples)
dev.off()

# look at H2A.Lap1 genes
toGrep <- c("H2afb3", "H2afb2", "Gm14920")
toGrep <- t2g[grep(paste(toGrep, collapse = "|"), t2g$ext_gene),]$target_id
kt.goi <- kt[grep(paste(toGrep, collapse = "|"), kt$target_id),]
