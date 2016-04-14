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
levels(s2c$condition) <- levels(s2c$condition)[c(2,1)]

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


