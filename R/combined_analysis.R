require(tximport)
require(readr)
require(ade4)
require(deepToolsUtils)

#-------------------------------------------------
# needs to be run in conjunction with "brain_expression_analysis.R"

combined <- list(combined_hemi = paste(pathPrefix, 
                                       "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_hemi", 
                                       sep = ""),
                 combined_wt = paste(pathPrefix,
                                     "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_wt",
                                     sep = ""))

edgeR_analysis_version <- 4
edgeR_analysis_output_combined <- paste("edgeR_analysis_V", edgeR_analysis_version, "_output_combined.rda", sep = "")


# edgeR Analysis using tximport to derive counts from kallisto transcript quantifcation --------
edgeR_analysis_output_combined <- paste("edgeR_analysis_V", edgeR_analysis_version, "_output_combined.rda", sep = "")

if(!file.exists(edgeR_analysis_output_combined)){
  edgeRProcessedDataCombined <- lapply(names(combined), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    options(mc.cores = mc.cores)
    sample_id <- dir(combined[[x]])
    kal_dirs <- sapply(sample_id, function(id) file.path(combined[[x]], id))
    
    # preparing factors for design matrix
    individual <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[1:2], collapse = "_"))))
    tissue <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[4], collapse = "_"))))
    strain <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[1], collapse = "_"))))
    s2c <- data.frame(sample = sample_id, 
                      individual = individual, 
                      tissue = tissue, 
                      strain = strain)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    # read in kallisto data with tximport
    s2c$files <- paste(s2c$path, "abundance.tsv", sep = "/")
    s2c <- s2c[order(s2c$tissue, s2c$strain),]
    files <- s2c$files
    names(files) <- s2c$sample
    txi <- tximport::tximport(files,
                              type = "kallisto",
                              tx2gene = t2g,
                              geneIdCol = "ens_gene",
                              txIdCol = "target_id",
                              reader = read_tsv)
    txi_lengthScaledTPM <- tximport::tximport(files,
                                              type = "kallisto",
                                              tx2gene = t2g,
                                              geneIdCol = "ens_gene",
                                              txIdCol = "target_id",
                                              reader = read_tsv,
                                              countsFromAbundance = "lengthScaledTPM")
    #---------------------------------------------
    # differential expression analysis using edgeR
    cts <- txi$counts
    normMat <- txi$length
    normMat <- normMat/exp(rowMeans(log(normMat)))
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
    # differential expression analysis using plain vanilla edgeR
    design <- model.matrix(~ tissue, data = s2c)
    y <- edgeR::DGEList(counts = cts, group = s2c$tissue)
    y$offset <- t(t(log(normMat)) + o)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmFit(y, design)
    lrtAll <- edgeR::glmLRT(fit, coef = 2:ncol(design))
    ttAll <- edgeR::topTags(lrtAll, n = 20000)
    annotatedTTAll <- merge(ttAll[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, sort = F)
    annotatedTTAll <- annotatedTTAll[order(annotatedTTAll$FDR), ]
    # return all objects
    return(list(txiGeneCounts = txi,
                txiScaledGeneCounts = txi_lengthScaledTPM,
                DGEList = y,
                glmFit = fit,
                glmLRT = lrtAll,
                topTagsAll = ttAll,
                AnnotatedTopTagsAll = annotatedTTAll))
  })
  names(edgeRProcessedDataCombined) <- names(combined)
  save(edgeRProcessedDataCombined, file = edgeR_analysis_output_combined)
} else {
  load(edgeR_analysis_output_combined)
}

# edgeR Analysis using HTSeq count tables --------
files <- list.files(path = dataPathHTSeq, full.names = T)
names(files) <- list.files(path = dataPathHTSeq, full.names = F)

if (!(!file.exists("htSeqCountMatrix.rda") && !file.exists("htSeqCountReport.rda"))){
  htSeqCountMatrix <- makeHTSeqCountMatrix(files)
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
  save(htSeqCountReport, file = "htSeqCountReport.rda")
} else {
  load("htSeqCountMatrix.rda")
  load("htSeqCountReport.rda")
}

plainEdgeROutput = "plain_edgeR_analysis_output.rda"

htSeqCountMatrixList <- list(wt_tissue = colnames(htSeqCountMatrix)[grep("wt", colnames(htSeqCountMatrix))],
                             hemi_tissue = colnames(htSeqCountMatrix)[grep("hemi", colnames(htSeqCountMatrix))],
                             wt_pfc_vs_ob = colnames(htSeqCountMatrix)[grep("PFC|OB", colnames(htSeqCountMatrix))][grep("wt", colnames(htSeqCountMatrix)[grep("PFC|OB", colnames(htSeqCountMatrix))])],
                             hemi_pfc_vs_ob = colnames(htSeqCountMatrix)[grep("PFC|OB", colnames(htSeqCountMatrix))][grep("hemi", colnames(htSeqCountMatrix)[grep("PFC|OB", colnames(htSeqCountMatrix))])],
                             wt_hippo_vs_pfc = colnames(htSeqCountMatrix)[grep("PFC|HIPPO", colnames(htSeqCountMatrix))][grep("wt", colnames(htSeqCountMatrix)[grep("PFC|HIPPO", colnames(htSeqCountMatrix))])],
                             hemi_hippo_vs_pfc = colnames(htSeqCountMatrix)[grep("PFC|HIPPO", colnames(htSeqCountMatrix))][grep("hemi", colnames(htSeqCountMatrix)[grep("PFC|HIPPO", colnames(htSeqCountMatrix))])],
                             tissue_hippo = colnames(htSeqCountMatrix)[grep("HIPPO", colnames(htSeqCountMatrix))],
                             tissue_pfc = colnames(htSeqCountMatrix)[grep("PFC", colnames(htSeqCountMatrix))],
                             tissue_ob = colnames(htSeqCountMatrix)[grep("OB", colnames(htSeqCountMatrix))])

if(!file.exists(plainEdgeROutput)){
  edgeRAnalysis <- lapply(names(htSeqCountMatrixList), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    countMatrix <- htSeqCountMatrix[, htSeqCountMatrixList[[x]]]
    sample_id <- unlist(lapply(strsplit(colnames(countMatrix), "_"), function(x) paste(x[1:4], collapse = "_")))
    if (unlist(strsplit(x, "_"))[1] == "tissue"){
      p <-3
    } else {
      p <- 4
    }
    condition <- as.factor(unlist(lapply(strsplit(colnames(countMatrix), "_"), function(z) z[p])))
    design <- model.matrix(~condition)
    y <- DGEList(counts = countMatrix, group = condition)
    y <- calcNormFactors(y, method="TMM")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = colnames(fit$design)[grep("condition", colnames(fit$design))])
    lrt <- glmQLFTest(fit, coef = colnames(fit$design)[grep("condition", colnames(fit$design))])
    tt <- topTags(lrt, n = 20000)
    annotatedTT <- merge(tt[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id")
    annotatedTT <- annotatedTT[order(annotatedTT$FDR),]
    # return all objects
    return(list(DGEList = y,
                qlmFit = fit,
                glmQLFt = qlf,
                glmQLRT = lrt,
                topTags = tt,
                AnnotatedTopTagsAll = annotatedTT))
  })
  names(edgeRAnalysis) <- names(htSeqCountMatrixList)
  save(edgeRAnalysis, file = plainEdgeROutput)
} else {
  load(plainEdgeROutput)
}

sapply(edgeRAnalysis, function(x){
  table(x[["AnnotatedTopTagsAll"]]$FDR < 0.1)
})


# edgeR on RUV-Seq normalised data ----------------------------------------
RUVSeqEdgeROutput = "RUVSeq_edgeR_analysis_output.rda"
if(!file.exists(RUVSeqEdgeROutput)){
  RUVSeqEdgeRAnalysisERCC <- lapply(names(htSeqCountMatrixList), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    original <- htSeqCountMatrix[, htSeqCountMatrixList[[x]]]
    filter <- apply(original, 1, function(y) length(y[y>5])>=2)
    filtered <- original[filter, ]
    genes <- rownames(filtered)[grep("ENS", rownames(filtered))]
    spikes <- rownames(filtered)[grep("ERCC", rownames(filtered))]
    #---------------------------------------------
    if (unlist(strsplit(x, "_"))[1] == "tissue"){
      p <-3
    } else {
      p <- 4
      }
    condition <- as.factor(unlist(lapply(strsplit(colnames(filtered), "_"), function(z) z[p])))
    # create expression set
    set <- newSeqExpressionSet(as.matrix(filtered),
                               phenoData = data.frame(condition, row.names=colnames(filtered)))
    #---------------------------------------------
    # data exploration
    rle <- plotRLE(set, outline = FALSE, col = colors[condition])
    pca <- plotPCA(set, col = colors[condition])
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
    # differential expression analysis using edgeR
    # here including the RUVg factor to account for "unwanted variation"
    design <- model.matrix(~condition + W_1, data = pData(set1))
    y <- DGEList(counts = counts(set1), group = condition)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = colnames(fit$design)[grep("condition", colnames(fit$design))])
    lrt <- glmQLFTest(fit, coef= colnames(fit$design)[grep("condition", colnames(fit$design))])
    tt <- topTags(lrt, n = 20000)
    annotatedTT <- merge(tt[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id")
    annotatedTT <- annotatedTT[order(annotatedTT$FDR),]
    #---------------------------------------------
    # differential expression analysis using all genes for RUVs
    differences <- makeGroups(condition)
    set2 <- RUVs(set, genes, k = 1, differences)
    rleRUVs <- plotRLE(set2, outline=FALSE, col=colors[condition])
    pcaRUVs <- plotPCA(set2, col=colors[condition], cex=1.2)
    design1 <- model.matrix(~condition + W_1, data = pData(set2))
    y1 <- DGEList(counts = counts(set2), group = condition)
    y1 <- calcNormFactors(y1, method="upperquartile")
    y1 <- estimateGLMCommonDisp(y1, design)
    y1 <- estimateGLMTagwiseDisp(y1, design)
    fit1 <- glmQLFit(y1, design)
    qlf1 <- glmQLFTest(fit1, coef = colnames(fit$design)[grep("condition", colnames(fit$design))])
    lrt1 <- glmQLFTest(fit1, coef = colnames(fit$design)[grep("condition", colnames(fit$design))])
    tt1 <- topTags(lrt1, n = 20000)
    annotatedTT1 <- merge(tt1[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id")
    annotatedTT1 <- annotatedTT1[order(annotatedTT1$FDR),]
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
                eSetRUVs = set2,
                DGEList = y,
                glmFit = fit,
                glmQLF = qlf,
                glmLRT = lrt,
                topTags = tt,
                AnnotatedTopTags = annotatedTT,
                DGEListRUVs = y1,
                RLEplotRUVs = rleRUVs,
                PCAplotRUVs = pcaRUVs,
                glmFitRUVs = fit1,
                glmLRTRUVs = lrt1,
                topTagsRUVs = tt1,
                AnnotatedTopTagsRUVs = annotatedTT1))
  })
  names(RUVSeqEdgeRAnalysisERCC) <- names(htSeqCountMatrixList)
  save(RUVSeqEdgeRAnalysisERCC, file = RUVSeqEdgeROutput)
} else {
  load(RUVSeqEdgeROutput)
}

sapply(RUVSeqEdgeRAnalysisERCC, function(x){
  table(x[["AnnotatedTopTags"]]$FDR < 0.1)
})

sapply(RUVSeqEdgeRAnalysisERCC, function(x){
  table(x[["AnnotatedTopTagsRUVs"]]$FDR < 0.1)
})

# SLEUTH analysis
sleuth_combined_analysis_version <- 2
sleuth_combined_analysis_output <- paste("sleuth_analysis_V", sleuth_combined_analysis_version, "_output_combined.rda", sep = "")

if(!file.exists(sleuth_combined_analysis_output)){
  sleuthProcessedDataCombined <- lapply(names(combined), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    options(mc.cores = mc.cores)
    sample_id <- dir(combined[[x]])
    kal_dirs <- sapply(sample_id, function(id) file.path(combined[[x]], id))
    condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[3], collapse = "_")))
    condition <- as.factor(condition)
    condition <- factor(condition, levels(condition)[c(2,1)])
    individual <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[1:2], collapse = "_"))))
    tissue <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[4], collapse = "_"))))
    s2c <- data.frame(sample = sample_id, condition = condition, individual = individual, tissue = tissue)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    s2c <- s2c[order(s2c$tissue, s2c$condition),]
    design <- model.matrix(~ tissue, data = s2c)
    
    #-----------------
    # transcript level
    so <- sleuth::sleuth_prep(s2c, ~ tissue, target_mapping = t2g)
    so <- sleuth::sleuth_fit(so, formula = design)
    so <- sleuth::sleuth_wt(so, "tissueOB")
    so <- sleuth::sleuth_wt(so, "tissuePFC")
    so <- sleuth::sleuth_fit(so, ~1, "reduced")
    so <- sleuth::sleuth_lrt(so, "reduced", "full")
    kt <- sleuth::kallisto_table(so)
    rt.ob <- sleuth::sleuth_results(so, "tissueOB")
    rt.pfc <- sleuth::sleuth_results(so, "tissuePFC")
    rt.ob <- rt.ob[order(rt.ob$qval),]
    rt.pfc <- rt.pfc[order(rt.pfc$qval),]
    kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt_wide) <- kt_wide$target_id
    kt_wide <- kt_wide[,-1]
    kt_wide <- kt_wide[, c(grep("HIPPO", colnames(kt_wide)),
                           grep("OB", colnames(kt_wide)),
                           grep("PFC", colnames(kt_wide)))]
    kt.pca <- ade4::dudi.pca(t(as.matrix(kt_wide)), center = T, scale = T, scannf = F, nf = 3)
    
    #-----------
    # gene level
    so.gene <- sleuth::sleuth_prep(s2c, ~ tissue, target_mapping = t2g, aggregation_column = "ens_gene")
    so.gene <- sleuth::sleuth_fit(so.gene, formula = design)
    so.gene <- sleuth::sleuth_wt(so.gene, "tissueOB")
    so.gene <- sleuth::sleuth_wt(so.gene, "tissuePFC")
    so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
    so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
    kt.gene <- sleuth::kallisto_table(so.gene)
    rt.gene.ob <- sleuth::sleuth_results(so.gene, "tissueOB")
    rt.gene.pfc <- sleuth::sleuth_results(so.gene, "tissuePFC")
    rt.gene.ob <- rt.gene.ob[order(rt.gene.ob$qval),]
    rt.gene.pfc <- rt.gene.pfc[order(rt.gene.pfc$qval),]
    kt.gene_wide <- tidyr::spread(kt.gene[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt.gene_wide) <- kt.gene_wide$target_id
    kt.gene_wide <- kt.gene_wide[,-1]
    kt.gene_wide <- kt.gene_wide[, c(grep("HIPPO", colnames(kt_wide)),
                                     grep("OB", colnames(kt_wide)),
                                     grep("PFC", colnames(kt_wide)))]
    kt.gene.pca <- ade4::dudi.pca(t(as.matrix(kt.gene_wide)), center = T, scale = T, scannf = F, nf = 3)
    
    return(list(sleuth_object = so,
                sleuth_results_ob = rt.ob,
                sleuth_results_pfc = rt.pfc,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                kallisto_pca = kt.pca,
                sleuth_object_gene = so.gene,
                sleuth_results_gene_ob = rt.gene.ob,
                sleuth_results_gene_pfc = rt.gene.pfc,
                kallisto_table_gene = kt.gene,
                kallisto_table_wide_gene = kt.gene_wide,
                kallisto_pca = kt.gene.pca))
  })
  names(sleuthProcessedDataCombined) <- names(combined)
  save(sleuthProcessedDataCombined, file = sleuth_combined_analysis_output)
} else {
  load(sleuth_combined_analysis_output)
}

sleuth_combined_analysis_compressed_output <- paste("sleuth_analysis_V", sleuth_combined_analysis_version, "_compressed_output_combined.rda", sep = "")
if (!file.exists(sleuth_combined_analysis_compressed_output)){
  sleuthProcessedDataCombinedCompressed <- lapply(names(sleuthProcessedDataCombined), function(x){
    sleuthProcessedDataCombined[[x]][grep("sleuth_object", names(sleuthProcessedDataCombined[[x]]), invert = T)]
  })
  names(sleuthProcessedDataCombinedCompressed) <- names(sleuthProcessedDataCombined)
  save(sleuthProcessedDataCombinedCompressed, file = sleuth_combined_analysis_compressed_output)
} else {
  load(sleuth_combined_analysis_compressed_output)
}

# extract scaled gene level data and combine for WT and hemi --------------
combinedMatrix <- cbind(edgeRProcessedDataCombined[["combined_wt"]][["txiScaledGeneCounts"]][["abundance"]],
                        edgeRProcessedDataCombined[["combined_hemi"]][["txiScaledGeneCounts"]][["abundance"]])

topGenesWT <- edgeRProcessedDataCombined[["combined_wt"]][["AnnotatedTopTagsAll"]][1:1000,]
topGenesHemi <- edgeRProcessedDataCombined[["combined_hemi"]][["AnnotatedTopTagsAll"]][1:1000,]
intersect(topGenesWT$Row.names, topGenesHemi$Row.names)

sapply(edgeRProcessedDataCombined, function(x){
  table(x[["AnnotatedTopTagsAll"]]$FDR < 0.1)
})

filter <- apply(combinedMatrix, 1, function(y) length(y[y>5])>=1)
filteredMatrix <- combinedMatrix[filter, ]
genes <- rownames(filteredMatrix)[grep("ENS", rownames(filteredMatrix))]
colnames(combinedMatrix) <- unlist(lapply(strsplit(colnames(filteredMatrix), "_"), function(x) paste(x[1:4], collapse = "_")))

sd1 <- apply(filteredMatrix, 1, sd)
sd1 <- sort(sd1, decreasing = T)
heatmap.3(log2(filteredMatrix[sd1 > 15, ] + 1), 
          trace = "none",
          labRow = F,
          cexCol = 0.65)


# perform PCA and plot diagram --------------------------------------------
combinedPCA <- ade4::dudi.pca(t(filteredMatrix), scannf = F, nf = 6)
fact <- unlist(lapply(strsplit(colnames(filteredMatrix), "_"), function(x) paste(x[3:4], collapse = "_")))
pdf("Combined_PCA_gene_level.pdf")
s.class(combinedPCA$li, fac = as.factor(fact), addaxes = T)
dev.off()
screeplot(combinedPCA)
s.arrow(combinedPCA$li)
