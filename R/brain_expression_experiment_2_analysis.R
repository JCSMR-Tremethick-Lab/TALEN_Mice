library(rhdf5)
library(sleuth)
library(biomaRt)
library(tidyr)
library(rtracklayer)
library(BiocParallel)
library(tximport)
library(readr)
library(RUVSeq)
library(data.table)
library(reshape)
library(reshape2)

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
  option(mc.cores = 8L)
  pathPrefix <- "~/mount/gduserv/"
  # biomaRt connection
  mouse <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
  attribs <- biomaRt::listAttributes(mouse)
} else {
  BPPARAM <- BiocParallel::MulticoreParam(workers = 32)
  options(mc.cores = 32L)
  pathPrefix <- "~"
  options(width = 137)
}

lDir <- function(x, y){
  paste(x, y, sep = "/")
}

if (dir.exists(lDir(pathPrefix, "Data/Tremethick/TALENs/Mus_musculus_brain_experiment_2/R_analysis"))){
  setwd(lDir(pathPrefix, "Data/Tremethick/TALENs/Mus_musculus_brain_experiment_2/R_analysis"))
} else {
  dir.create(lDir(pathPrefix, "Data/Tremethick/TALENs/Mus_musculus_brain_experiment_2/R_analysis"))
  setwd(lDir(pathPrefix, "Data/Tremethick/TALENs/Mus_musculus_brain_experiment_2/R_analysis"))
}

sleuth_analysis_version <- 1
sleuth_analysis_output <- paste("sleuth_analysis_V", sleuth_analysis_version, "_output.rda", sep = "")
sleuth_analysis_output_combined <- paste("sleuth_analysis_V", sleuth_analysis_version, "_output_combined.rda", sep = "")
sleuth_analysis_outputCompressed <- paste("sleuth_analysis_V", sleuth_analysis_version, "_Compressed_output.rda", sep = "")

previous_sleuth_analysis_version <- sleuth_analysis_version - 1
previous_sleuth_analysis_output <- paste("sleuth_analysis_V", previous_sleuth_analysis_version, "_output.rda", sep = "")
previous_sleuth_analysis_output_combined <- paste("sleuth_analysis_V", previous_sleuth_analysis_version, "_output_combined.rda", sep = "")

kallisto_base_dir <- lDir(pathPrefix, "Data/Tremethick/TALENs/Mus_musculus_brain_experiment_2/processed_data/GRCm38_ensembl84/kallisto")


# use Ensembl 84 for annotation
if (!file.exists("t2g.rda")){
  mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
  attribs <- listAttributes(mouse)
  # annotate transcripts
  # Ensembl 84 includes version number in FASTA IDs, therefore have to conactenate them in, other wise mapping does not work
  t2g <- getBM(attributes = c("ensembl_transcript_id", 
                              "ensembl_gene_id", 
                              "external_gene_name",
                              "description",
                              "version", 
                              "transcript_version"), 
               mart = mouse)
  t2g <- as.data.table(t2g)
  t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, description = description)
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
if (file.exists(previous_sleuth_analysis_output)) {file.remove(previous_sleuth_analysis_output); file.remove(previous_sleuth_analysis_output_combined)}


base_dirs <- list(pfc = kallisto_base_dir)
if (!file.exists(sleuth_analysis_output)){
  sleuthProcessedData <- lapply(names(base_dirs), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    sample_id <- dir(base_dirs[[x]])
    kal_dirs <- sapply(sample_id, function(id) file.path(base_dirs[[x]], id))
    condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[5], collapse = "_")))
    condition <- as.factor(condition)
    condition <- relevel(condition, ref = "wt")
    treatment <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[6], collapse = "_")))
    treatment <- as.factor(treatment)
    treatment <- relevel(treatment, ref = "naive")
    s2c <- data.frame(sample = sample_id, condition = condition, treatment = treatment)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    designMatrix <- model.matrix(~ condition * treatment, data = s2c)
    
    # here we test for the interaction between mutant & wt and try to determine which variability is explained by fear-conditionint alone
    so <- sleuth::sleuth_prep(s2c, ~ condition * treatment, target_mapping = t2g)
    so <- sleuth::sleuth_fit(so, designMatrix)
    so <- sleuth::sleuth_wt(so, "conditionmut")
    so <- sleuth::sleuth_wt(so, "treatmentFC")
    so <- sleuth::sleuth_wt(so, "conditionmut:treatmentFC")
    so <- sleuth::sleuth_fit(so, ~1, "reduced")
    so <- sleuth::sleuth_lrt(so, "reduced", "full")
    kt <- sleuth::kallisto_table(so, use_filtered = T)
    rt <- sleuth::sleuth_results(so, "conditionmut:treatmentFC", show_all = F)
    
    # do same analysis on gene level ----------------------------------------
    so.gene <- sleuth::sleuth_prep(s2c, ~ condition * treatment, target_mapping = t2g, aggregation_column = "ens_gene")
    so.gene <- sleuth::sleuth_fit(so.gene, designMatrix)
    so.gene <- sleuth::sleuth_wt(so.gene, "conditionmut")
    so.gene <- sleuth::sleuth_wt(so.gene, "treatmentFC")
    so.gene <- sleuth::sleuth_wt(so.gene, "conditionmut:treatmentFC")
    so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
    so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
    kt.gene <- sleuth::kallisto_table(so.gene, use_filtered = T)
    kt_wide.gene <- tidyr::spread(kt.gene[, c("target_id", "sample", "tpm")], sample, tpm)
    kt_wide.gene <- data.table::as.data.table(merge(kt_wide.gene, subset(t2g, select = c("target_id", "ens_gene")), all.x = TRUE, all.y = FALSE))
    kt_wide.gene <- kt_wide.gene[,lapply(.SD,sum),by=ens_gene, .SDcols = 2:13]
    kt.gene <- melt(kt_wide.gene)
    kt.gene$groups <- unlist(lapply(strsplit(as.character(kt.gene$variable), "_"), function(x) paste(x[5:6], collapse = "_")))
    kt_wide.gene <- data.table::as.data.table(merge(kt_wide.gene, unique(subset(t2g, select =  c("ens_gene", "ext_gene", "description")), by = "ens_gene"), all.x = TRUE, all.y = FALSE))
    rt.gene <- sleuth::sleuth_results(so.gene, "conditionmut:treatmentFC", show_all = F)
    rt.gene <- as.data.table(merge(rt.gene, unique(subset(t2g, select = c("ens_gene", "ext_gene", "description")), by = "ens_gene"), all.x = T, all.y = F, by.x = "target_id", by.y = "ens_gene"))
    
    # alternatively, and statistically probably not so clean, test  --------
    # look at difference between naive WT vs MUT
    s2c.condition.naive <- s2c[grep("naive", s2c$treatment),]
    s2c.condition.naive <- subset(s2c.condition.naive, select = -(treatment))
    so.condition.naive <- sleuth::sleuth_prep(s2c.condition.naive, ~ condition, target_mapping = t2g)
    designMatrix.condition.naive <- model.matrix(~ condition, data = s2c.condition.naive)
    so.condition.naive <- sleuth::sleuth_fit(so.condition.naive, designMatrix.condition.naive)
    so.condition.naive <- sleuth::sleuth_wt(so.condition.naive, "conditionmut")
    so.condition.naive <- sleuth::sleuth_fit(so.condition.naive, ~1, "reduced")
    so.condition.naive <- sleuth::sleuth_lrt(so.condition.naive, "reduced", "full")
    rt.condition.naive <- sleuth::sleuth_results(so.condition.naive, "conditionmut", show_all = F)
    rt.condition.naive <- as.data.table(rt.condition.naive)
    rt.condtion.naive <- rt.condition.naive[!(is.na(qval))]
    
    rt.condition.naive <- merge(rt.condition.naive, unique(subset(t2g, select = c("ens_gene", "ext_gene", "description")), by = "ens_gene"), all.x = T, all.y = F, by.x = "target_id", by.y = "ens_gene")
    
    so.condition.naive.gene <- sleuth::sleuth_prep(s2c.condition.naive, ~ condition, target_mapping = t2g, aggregation_column = "ens_gene")
    so.condition.naive.gene <- sleuth::sleuth_fit(so.condition.naive.gene, designMatrix.condition.naive)
    so.condition.naive.gene <- sleuth::sleuth_wt(so.condition.naive.gene, "conditionmut")
    so.condition.naive.gene <- sleuth::sleuth_fit(so.condition.naive.gene, ~1, "reduced")
    so.condition.naive.gene <- sleuth::sleuth_lrt(so.condition.naive.gene, "reduced", "full")

    # look at difference FC WT vs MUT
    s2c.condition.fc <- s2c[grep("FC", s2c$treatment),]
    s2c.condition.fc <- subset(s2c.condition.fc, select = -(treatment))
    so.condition.fc <- sleuth::sleuth_prep(s2c.condition.fc, ~ condition, target_mapping = t2g)
    designMatrix.condition.fc <- model.matrix(~ condition, data = s2c.condition.fc)
    so.condition.fc <- sleuth::sleuth_fit(so.condition.fc, designMatrix.condition.fc)
    so.condition.fc <- sleuth::sleuth_wt(so.condition.fc, "conditionmut")
    so.condition.fc <- sleuth::sleuth_fit(so.condition.fc, ~1, "reduced")
    so.condition.fc <- sleuth::sleuth_lrt(so.condition.fc, "reduced", "full")
    rt.condition.fc <- sleuth::sleuth_results(so.condition.fc, "conditionmut")
    rt.condition.fc <- as.data.table(rt.condition.fc)
    rt.condition.fc <- rt.condition.fc[!(is.na(qval))]
    
    s2c.treatment.wt <- s2c[grep("wt", s2c$condition),]
    s2c.treatment.wt <- subset(s2c.treatment.wt, select = -(condition))
    so.treatment.wt <- sleuth::sleuth_prep(s2c.treatment.wt, ~ treatment, target_mapping = t2g)
    designMatrix.treatment.wt <- model.matrix(~ treatment, data = s2c.treatment.wt)
    so.treatment.wt <- sleuth::sleuth_fit(so.treatment.wt, designMatrix.treatment.wt)
    so.treatment.wt <- sleuth::sleuth_wt(so.treatment.wt, "treatmentFC")
    so.treatment.wt <- sleuth::sleuth_fit(so.treatment.wt, ~1, "reduced")
    so.treatment.wt <- sleuth::sleuth_lrt(so.treatment.wt, "reduced", "full")
    rt.treatment.wt <- sleuth::sleuth_results(so.treatment.wt, "treatmentFC")
    rt.treatment.wt <- as.data.table(rt.treatment.wt)
    rt.treatment.wt <- rt.treatment.wt[!(is.na(qval))]
    
    s2c.treatment.mut <- s2c[grep("mut", s2c$condition),]
    s2c.treatment.mut <- subset(s2c.treatment.mut, select = -(condition))
    so.treatment.mut <- sleuth::sleuth_prep(s2c.treatment.mut, ~ treatment, target_mapping = t2g)
    designMatrix.treatment.mut <- model.matrix(~ treatment, data = s2c.treatment.mut)
    so.treatment.mut <- sleuth::sleuth_fit(so.treatment.mut, designMatrix.treatment.mut)
    so.treatment.mut <- sleuth::sleuth_wt(so.treatment.mut, "treatmentFC")
    so.treatment.mut <- sleuth::sleuth_fit(so.treatment.mut, ~1, "reduced")
    so.treatment.mut <- sleuth::sleuth_lrt(so.treatment.mut, "reduced", "full")
    rt.treatment.mut <- sleuth::sleuth_results(so.treatment.mut, "treatmentFC")
    rt.treatment.mut <- as.data.table(rt.treatment.mut)
    rt.treatment.mut <- rt.treatment.mut[!(is.na(qval))]
    
    rt <- rt[order(rt$qval),]
    kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt_wide) <- kt_wide$target_id
    kt_wide <- kt_wide[,-1]
    kt_wide <- kt_wide[, c(grep("wt", colnames(kt_wide)),
                           grep("hemi", colnames(kt_wide)))]
    kt.pca <- ade4::dudi.pca(t(as.matrix(kt_wide)), center = T, scale = T, scannf = F, nf = 5)
    return(list(sleuth_object = so,
                sleuth_results = rt,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                kallisto_pca = kt.pca,
                so.treatment.mut,
                so.treatment.wt,
                so.condition.fc,
                so.condition.naive))
    })
  names(sleuthProcessedData) <- names(base_dirs)
  save(sleuthProcessedData, file = sleuth_analysis_output)
} else {
  load(sleuth_analysis_output)
}

write.csv(rt.condition.naive, "CSV_export/sleuth_results_naive_mutant_vs_wt.csv")
write.csv(rt.condition.fc, "CSV_export/sleuth_results_fear_mutant_vs_wt.csv")
write.csv(rt.treatment.mut, "CSV_export/sleuth_results_mutant_naive_vs_fear.csv")
write.csv(rt.treatment.wt, "CSV_export/sleuth_results_wildtype_naive_vs_fear.csv")


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
  condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[5:6], collapse = "_")))
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
                            tx2gene = t2g[,c(1,3)],
                            reader = read_tsv)
  pca1 <- ade4::dudi.pca(t(txi$abundance), scannf = F, nf = 6)
  ade4::s.class(pca1$li, fac = condition )
  ade4::s.arrow(pca1$li)
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
