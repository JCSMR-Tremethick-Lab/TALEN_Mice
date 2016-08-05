combined <- list(combined_hemi = paste(pathPrefix, "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_hemi", sep = ""),
                 combined_wt = paste(pathPrefix, "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_wt", sep = ""))

edgeR_analysis_output_combined <- paste("edgeR_analysis_V", edgeR_analysis_version, "_output_combined.rda", sep = "")

if(!file.exists(edgeR_analysis_output_combined)){
  edgeRProcessedDataCombined <- lapply(names(combined), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    options(mc.cores = mc.cores)
    sample_id <- dir(combined[[x]])
    kal_dirs <- sapply(sample_id, function(id) file.path(combined[[x]], id))
    individual <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[1:2], collapse = "_"))))
    tissue <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[4], collapse = "_"))))
    s2c <- data.frame(sample = sample_id, individual = individual, tissue = tissue)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    s2c$files <- paste(s2c$path, "abundance.tsv", sep = "/")
    s2c <- s2c[order(s2c$tissue),]
    files <- s2c$files
    names(files) <- s2c$sample
    # read in kallisto data with tximport
    txi <- tximport::tximport(files,
                              type = "kallisto",
                              tx2gene = t2g,
                              geneIdCol = "ens_gene",
                              txIdCol = "target_id",
                              reader = read_tsv)
    original <- round(txi$counts, 0)
    filter <- apply(original, 1, function(y) length(y[y>5])>=2)
    filtered <- original[filter, ]
    genes <- rownames(filtered)[grep("ENS", rownames(filtered))]
    spikes <- rownames(filtered)[grep("ERCC", rownames(filtered))]
    print(paste("After filtering", table(filter)[2], "genes including", length(spikes), "ERCC spike-ins were retained", sep = " "))
    #---------------------------------------------
    # create expression set
    if (grepl("combined", x)){
      set <- EDASeq::newSeqExpressionSet(as.matrix(filtered),
                                         phenoData = data.frame(s2c$tissue, s2c$individual, row.names=colnames(filtered)))
    } else {
      set <- EDASeq::newSeqExpressionSet(as.matrix(filtered),
                                         phenoData = data.frame(s2c$condition, row.names=colnames(filtered)))
    }
    #---------------------------------------------
    # data exploration
    rle <- EDASeq::plotRLE(set, outline = FALSE, col = colors[s2c$tissue])
    pca <- EDASeq::plotPCA(set, col = colors[s2c$tissue])
    #---------------------------------------------
    # RUVseq using spike ins
    set1 <- RUVSeq::RUVg(set, spikes, k = 1)
    rleRUV <- EDASeq::plotRLE(set1, outline=FALSE, col=colors[s2c$tissue])
    pcaRUV <- EDASeq::plotPCA(set1, col=colors[s2c$tissue], cex=1.2)
    #---------------------------------------------
    # differential expression analysis using edgeR
    # here including the RUVg factor to account for "unwanted variation"
    cts <- txi$counts[filter,]
    normMat <- txi$length[filter,]
    normMat <- normMat/exp(rowMeans(log(normMat)))
    design <- model.matrix(~ s2c.tissue + s2c.individual + W_1, data = pData(set1))
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
    y <- edgeR::DGEList(counts = cts, group = s2c$tissue)
    y$offset <- t(t(log(normMat)) + o)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmFit(y, design)
    lrt <- edgeR::glmLRT(fit, coef=4)
    tt <- edgeR::topTags(lrt, n = 20000)
    annotatedTT <- merge(tt[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id")
    annotatedTT <- annotatedTT[order(annotatedTT$FDR),]
    #---------------------------------------------
    # differential expression analysis using plain vanilla edgeR
    design1 <- model.matrix(~ s2c.tissue , data = pData(set))
    y1 <- edgeR::DGEList(counts = counts(set), group = s2c$tissue)
    #y1 <- edgeR::calcNormFactors(y1, method="upperquartile")
    y1$offset <- t(t(log(normMat)) + o)
    y1 <- edgeR::estimateDisp(y1, design1)
    fit1 <- edgeR::glmFit(y1, design1)
    lrtOB_vs_HIPPO <- edgeR::glmLRT(fit1, coef = 2)
    lrtPFC_vs_HIPPO <- edgeR::glmLRT(fit1, coef = 3)
    lrtAll <- edgeR::glmLRT(fit1, coef = 2:3)
    ttAll <- edgeR::topTags(lrtAll, n = 20000)
    annotatedTTAll <- merge(ttAll[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, sort = F)
    annotatedTTAll <- annotatedTTAll[order(annotatedTTAll$FDR), ]
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
  names(edgeRProcessedDataCombined) <- names(combined)
  save(edgeRProcessedDataCombined, file = edgeR_analysis_output_combined)
} else {
  load(edgeR_analysis_output_combined)
}


# TODO edit
if(!file.exists(sleuth_analysis_output_combined)){
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
    # transcript level
    so <- sleuth::sleuth_prep(s2c, ~ tissue + individual + tissue, target_mapping = t2g)
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
  names(sleuthProcessedDataCombined) <- names(combined)
  save(sleuthProcessedDataCombined, file = sleuth_analysis_output_combined)
} else {
  load(sleuth_analysis_output_combined)
}

edgeR_analysis_output_combined <- paste("edgeR_analysis_V", edgeR_analysis_version, "_output_combined.rda", sep = "")

if(!file.exists(edgeR_analysis_output_combined)){
  edgeRProcessedDataCombined <- lapply(names(combined), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    options(mc.cores = mc.cores)
    sample_id <- dir(combined[[x]])
    kal_dirs <- sapply(sample_id, function(id) file.path(combined[[x]], id))
    condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[3], collapse = "_")))
    condition <- as.factor(condition)
    individual <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[1:2], collapse = "_"))))
    tissue <- as.factor(unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[4], collapse = "_"))))
    s2c <- data.frame(sample = sample_id, condition = condition, individual = individual, tissue = tissue)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    s2c$files <- paste(s2c$path, "abundance.tsv", sep = "/")
    s2c$condition <- factor(s2c$condition, levels(s2c$condition)[c(2,1)])
    s2c <- s2c[order(s2c$tissue, s2c$condition),]
    files <- s2c$files
    names(files) <- s2c$sample
    # read in kallisto data with tximport
    txi <- tximport::tximport(files,
                              type = "kallisto",
                              tx2gene = t2g,
                              geneIdCol = "ens_gene",
                              txIdCol = "target_id",
                              reader = read_tsv)
    original <- round(txi$counts, 0)
    filter <- apply(original, 1, function(y) length(y[y>5])>=2)
    filtered <- original[filter, ]
    genes <- rownames(filtered)[grep("ENS", rownames(filtered))]
    spikes <- rownames(filtered)[grep("ERCC", rownames(filtered))]
    #---------------------------------------------
    # create expression set
    if (x == "combined"){
      set <- EDASeq::newSeqExpressionSet(as.matrix(filtered),
                                         phenoData = data.frame(s2c$condition, s2c$tissue, s2c$individual, row.names=colnames(filtered)))
    } else {
      set <- EDASeq::newSeqExpressionSet(as.matrix(filtered),
                                         phenoData = data.frame(s2c$condition, row.names=colnames(filtered)))
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
    design <- model.matrix(~ s2c.condition + s2c.tissue + W_1, data = pData(set1))
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
    y <- edgeR::DGEList(counts = counts(set1), group = s2c$condition)
    y$offset <- t(t(log(normMat)) + o)
    y <- edgeR::estimateGLMCommonDisp(y, design)
    y <- edgeR::estimateGLMTagwiseDisp(y, design)
    fit <- edgeR::glmFit(y, design)
    lrt <- edgeR::glmLRT(fit, coef=5)
    tt <- edgeR::topTags(lrt, n = 20000)
    annotatedTT <- merge(tt[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id")
    annotatedTT <- annotatedTT[order(annotatedTT$FDR),]
    #---------------------------------------------
    # differential expression analysis using plain vanilla edgeR
    design1 <- model.matrix(~condition + tissue , data = pData(set))
    y1 <- edgeR::DGEList(counts = counts(set), group = tissue)
    y1 <- edgeR::calcNormFactors(y1, method="upperquartile")
    y1 <- edgeR::estimateGLMCommonDisp(y1, design1)
    y1 <- edgeR::estimateGLMTagwiseDisp(y1, design1)
    fit1 <- edgeR::glmFit(y1, design1)
    lrt1 <- edgeR::glmLRT(fit1, coef=1)
    tt1 <- edgeR::topTags(lrt1, n = 5000)
    annotatedTT1 <- merge(tt1[[1]], ensGenes, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, sort = F)
    annotatedTT1 <- annotatedTT1[sort(annotatedTT1$FDR), ]
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
  names(edgeRProcessedDataCombined) <- names(combined)
  save(edgeRProcessedDataCombined, file = edgeR_analysis_output_combined)
  } else {
  load(edgeR_analysis_output_combined)
}

