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

# setting working directory and data sources ------------------------------
if (amILocal("JCSMR027564ML")){
  mount <- system("mount", intern = T)
  if (length(grep("gduserv", mount)) == 0) {system("sshfs skurscheid@gduserv.anu.edu.au: ~/mount/gduserv/")}
  BPPARAM <- BiocParallel::MulticoreParam(workers = 7)
  mc.cores <- 8
  setwd("~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/DEXSeq/count/"
  kallisto_base_dir <- "~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto"
  # biomaRt connection
  mouse <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
  attribs <- biomaRt::listAttributes(mouse)
  pathPrefix = "~/mount/gduserv"
  if (file.exists("~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.DEXSeq.gtf")) {
    flattenedfile <- "~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.DEXSeq.gtf"
  } else {
    stop("GTF file missing!")
  }
} else {
  BPPARAM <- BiocParallel::MulticoreParam(workers = 16)
  mc.cores <- 16
  pathPrefix = "~"
  setwd("~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/DEXSeq/count/"
  flattenedfile <- "~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.DEXSeq.gtf"
  kallisto_base_dir <- "~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto"
}

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

# transcript level differential expression analysis -----------------------------
# run analysis by brain region so that only pair-wise comparison

base_dirs <- c(OB = paste(pathPrefix, "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_ob", sep = ""),
               PFC = paste(pathPrefix, "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_pfc", sep = ""),
               HIPPO = paste(pathPrefix, "/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto_hippo", sep = ""))

sleuth_analysis_version <- 2
sleuth_analysis_output <- paste("sleuth_analysis_V", sleuth_analysis_version, "_output.rda", sep = "")

if (!file.exists(sleuth_analysis_output)){
  sleuthProcessedData <- lapply(names(base_dirs), function(x){
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
    kt <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt) <- kt$target_id
    kt <- kt[,-1]
    kt.pca <- ade4::dudi.pca(t(as.matrix(kt)), center = T, scale = T, scannf = F, nf = 3)
    return(list(sleuth_object = so,
                kallisto_table = kt,
                kallisto_pca = kt.pca))
    })
  names(sleuthProcessedData) <- names(base_dirs)
  save(sleuthProcessedData, file = sleuth_analysis_output)
} else {
  load(sleuth_analysis_output)
}


s.class(kt_gene_pca$li, fac = as.factor(unlist(lapply(strsplit(rownames(kt_gene_pca$tab), "_"), function(x) x[3]))))

# run sleuth on all samples -----------------------------------------------

# creating data.frame for experimental design and file names --------------
sample_id <- dir(kallisto_base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(kallisto_base_dir, id))
condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) x[3]))
tissue <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) x[4]))
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
x <- s2c$condition
x <- factor(x, levels(x)[c(4,5,6,1,2,3)])
s2c$condition <- x
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)

# get raw data - count table needed for RUV analysis etc.
kt.raw <- kallisto_table(so, normalize = F)
kt.raw_counts <- tidyr::spread(kt.raw[, c("target_id", "sample", "est_counts")], sample, est_counts)
rownames(kt.raw_counts) <- unlist(lapply(strsplit(kt.raw_counts$target_id, "\\."), function(x) x[1]))
kt.raw_counts <- kt.raw_counts[, -1]
save(kt.raw_counts, file = "kt.raw_counts.rda")
so <- sleuth_fit(so)
so <- sleuth_wt(so, "conditionhemi_HIPPO")
so <- sleuth_fit(so, ~1, "reduced")


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
