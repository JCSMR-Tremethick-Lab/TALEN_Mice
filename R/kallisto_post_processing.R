require(rhdf5)
require(sleuth)
require(biomaRt)

source("~/Development/GeneralPurpose/R/heatmap.3.R")


# setting working directory and data sources ------------------------------
setwd("~/Data/Tremethick/TALENs/NB501086_0047_TSoboleva_JCSMR_stranded_RNASeq/R_analysis")
base_dir <- "~/Data/Tremethick/TALENs/NB501086_0047_TSoboleva_JCSMR_stranded_RNASeq/processed_data"

# creating data.frame for experimental design and file names --------------
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
s2c <- data.frame(sample = sample_id, condition = c("KO", "WT", "WT", "KO", "WT", "KO"))
s2c <- dplyr::mutate(s2c, path = kal_dirs)
x <- s2c$condition
x <- factor(x, levels(x)[c(2,1)])
s2c$condition <- x

# transcript level differential expression analysis -----------------------
so <- sleuth_prep(s2c, ~ condition)
so <- sleuth_fit(so)
so <- sleuth_wt(so, "conditionKO")


# gene level differential expression analysis -----------------------------
#use Ensembl 74 for annotation
mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2013.archive.ensembl.org")
attribs <- listAttributes(mouse)
# annotate transcripts
t2g <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_id"), mart = mouse)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_id)
# re-run sleuth
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_wt(so, "conditionKO")
results_table <- sleuth_results(so, "conditionKO")
results_table$FC_estimated <- log2(exp(results_table$b))

kt <- kallisto_table(so)

toGrep <- c("H2afb3", "H2afb2", "Gm14920")
toGrep <- t2g[grep(paste(toGrep, collapse = "|"), t2g$ext_gene),]$target_id
kt.goi <- kt[grep(paste(toGrep, collapse = "|"), kt$target_id),]
kt.goi$sample <- factor(kt.goi$sample)
kt.goi <- spread(kt.goi[, c("target_id", "sample", "tpm")], sample, tpm)
rownames(kt.goi) <- kt.goi$target_id
kt.goi$gene <- t2g[rownames(kt.goi),]$ext_gene
kt.goi <- kt.goi[, c("NMG3-62wt_S2", "NMG3-74wt_S3", "NMG3-76wt_S5", "NMG3-60hemi_S1", "NMG3-75hemi_S4", "NMG3-77hemi_S6", "gene")]

write.csv(kt.goi, "KO_Gene_expression_levels.csv")
write.csv(results_table, "Differential_Expression_WT_vs_KO.csv")
