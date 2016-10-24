require(ade4)
require(deepToolsUtils)
require(rhdf5)
require(sleuth)
require(biomaRt)
require(tidyr)
require(rtracklayer)
require(BiocParallel)

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
  # biomaRt connection
  mouse <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
  attribs <- biomaRt::listAttributes(mouse)
  pathPrefix = "~/mount/gduserv"
} else {
  BPPARAM <- BiocParallel::MulticoreParam(workers = 16)
  mc.cores <- 16L
  pathPrefix <- "~"
}

setwd(paste(pathPrefix, 
            "Data/Tremethick/Brain/SRA_mmus_brain/SRP047108/R_Analysis",
            sep = "/"))

kallisto_base_dir <- paste(pathPrefix,
                           "Data/Tremethick/Brain/SRA_mmus_brain/SRP047108/processed_data/GRCm38_ensembl84_ERCC/kallisto_se",
                           sep = "/")

t2gFile <- paste(pathPrefix,
                 "Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/t2g.rda",
                 sep = "/")

# use Ensembl 84 for annotation
if (!file.exists(t2gFile)){
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
  load(t2gFile)
}

experiments <- list(SRP047108 = paste(pathPrefix, "Data/Tremethick/Brain/SRA_mmus_brain/SRP047108/processed_data/GRCm38_ensembl84_ERCC/kallisto_se/", sep = "/"))

# actual analysis ---------------------------------------------------------

sleuth_output = "sleuthOutput.rda"
if(!file.exists(sleuth_output)){
  sleuthProcessedData <- lapply(names(experiments), function(x){
    print(paste("Processing", x, "samples", sep = " "))
    options(mc.cores = mc.cores)
    sample_id <- dir(experiments[[x]])
    kal_dirs <- sapply(sample_id, function(id) file.path(experiments[[x]], id))
    condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[1], collapse = "_")))
    condition <- as.factor(condition)
    s2c <- data.frame(sample = sample_id, condition = condition)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    s2c <- s2c[order(s2c$condition),]
    design <- model.matrix(~ condition, data = s2c)
    # transcript level
    so <- sleuth::sleuth_prep(s2c, ~ condition, target_mapping = t2g)
    so <- sleuth::sleuth_fit(so)
    so <- sleuth::sleuth_wt(so, "conditionPFC")
    so <- sleuth::sleuth_fit(so, ~1, "reduced")
    so.lrt <- sleuth::sleuth_lrt(so, "reduced", "full")
    kt <- sleuth::kallisto_table(so)
    rt <- sleuth::sleuth_results(so, "conditionPFC")
    rt <- rt[order(rt$qval),]
    kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt_wide) <- kt_wide$target_id
    kt_wide <- kt_wide[,-1]
    kt.pca <- ade4::dudi.pca(t(as.matrix(kt_wide)), center = T, scale = T, scannf = F, nf = 6)
    
    # gene level analysis
    so.gene <- sleuth::sleuth_prep(s2c, ~ condition, target_mapping = t2g, aggregation_column = "ens_gene")
    so.gene <- sleuth::sleuth_fit(so.gene)
    so.gene <- sleuth::sleuth_wt(so.gene, "conditionPFC")
    so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
    so.gene.lrt <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
    kt.gene <- sleuth::kallisto_table(so.gene)
    rt.gene <- sleuth::sleuth_results(so.gene, "conditionPFC")
    rt.gene <- rt[order(rt.gene$qval),]
    kt.gene_wide <- tidyr::spread(kt.gene[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt.gene_wide) <- kt.gene_wide$target_id
    kt.gene_wide <- kt.gene_wide[,-1]
    kt.gene.pca <- ade4::dudi.pca(t(as.matrix(kt.gene_wide)), center = T, scale = T, scannf = F, nf = 6)
    
    # return to the list
    return(list(sleuth_object = so,
                sleuth_results = rt,
                sleuth_lrt = so.lrt,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                kallisto_pca = kt.pca,
                sleuth_object_gene = so.gene,
                sleuth_results_gene = rt.gene,
                sleuth_lrt_gene = so.gene.lrt,
                kallisto_table_gene = kt.gene,
                kallisto_table_wide_gene = kt.gene_wide,
                kallisto_pca = kt.gene.pca))
  })
  names(sleuthProcessedData) <- names(experiments)
  save(sleuthProcessedData, file = sleuth_output)
} else {
  load(sleuth_output)
}

ade4::s.class(kt.pca$li, fac = s2c$condition)
