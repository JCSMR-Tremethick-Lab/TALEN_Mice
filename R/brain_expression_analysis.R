require(rhdf5)
require(sleuth)
require(biomaRt)
require(tidyr)
require(rtracklayer)

source("~/Development/GeneralPurpose/R/heatmap.3.R")


# setting working directory and data sources ------------------------------
setwd("~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
base_dir <- "~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto"

# creating data.frame for experimental design and file names --------------
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[3:4], collapse = "_")))
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
x <- s2c$condition
x <- factor(x, levels(x)[c(4,5,6,1,2,3)])
s2c$condition <- x

# transcript level differential expression analysis -----------------------------
# use Ensembl 84 for annotation
mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "asia.ensembl.org")
attribs <- listAttributes(mouse)
# annotate transcripts
# Ensembl 84 includes version number in FASTA IDs, therefore have to conactenate them in, other wise mapping does not work
t2g <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "version", "transcript_version"), mart = mouse)
t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# add ERCC spike ins
ercc <- import("~/mount/gduserv/Data/References/Transcriptomes/ERCC/ERCC92.gtf")
ercc.df <- mcols(ercc)
ercc.df <- data.frame(ercc.df[, c("transcript_id", "gene_id", "gene_id")])
colnames(ercc.df) <- c("target_id", "ens_gene", "ext_gene")
t2g <- rbind(t2g, ercc.df)

rownames(t2g) <- t2g$target_id

save(t2g, file = "~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/kallisto/t2g.rda")

load("~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
# re-run sleuth
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_wt(so, "conditionhemi_HIPPO")
so <- sleuth_fit(so, ~1, "reduced")

results_table.hippo <- sleuth_results(so, "conditionhemi_HIPPO")
results_table$FC_estimated <- log2(exp(results_table$b))

kt <- kallisto_table(so)

load("~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis/kt.rda")
kt <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
rownames(kt) <- kt$target_id
# check ERCC expression
ercc_tab <- read.table("~/mount/gduserv/Data/References/Transcriptomes/ERCC/cms_095046.txt", header = T, as.is = T, sep = "\t")
rownames(ercc_tab) <- ercc_tab$ERCC.ID
ercc_exp <- kt[grep("ERCC", kt$target_id),]
i1 <- intersect(rownames(ercc_exp), rownames(ercc_tab))

pdf("ERCC_spike_ins.pdf", paper = "a4")
par(mfrow = c(3,2))
lapply(seq_along(1:18), function(x){
  x <- x+1
  s <- colnames(ercc_exp)[x]
  cor <- cor(log2(ercc_tab[i1, "concentration.in.Mix.1..attomoles.ul."] + 1), log2(ercc_exp[i1, x] + 1)) 
  sm <- summary(lm(log2(ercc_exp[i1, x] + 1) ~ log2(ercc_tab[i1, "concentration.in.Mix.1..attomoles.ul."] + 1)))
  plot(log2(ercc_tab[i1, "concentration.in.Mix.1..attomoles.ul."] + 1), log2(ercc_exp[i1, x] + 1),
       main = paste(s),
       xlab = "Log2(aM) Mix 1",
       ylab = "Log2(tpm + 1)")
  lines(lowess(log2(ercc_exp[i1, x] + 1) ~ log2(ercc_tab[i1, "concentration.in.Mix.1..attomoles.ul."] + 1)),col="green3")
  abline(lm(log2(ercc_exp[i1, x] + 1) ~ log2(ercc_tab[i1, "concentration.in.Mix.1..attomoles.ul."] + 1)),col="red3")
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


  
# run analysis by brain region so that only pair-wise comparison will be done
s2c.pfc <- s2c[grep("PFC", s2c$condition),]
s2c.ob <- s2c[grep("OB", s2c$condition),]
s2c.hippo <- s2c[grep("HIPPO", s2c$condition),]

so.pfc <- sleuth_prep(s2c.pfc, ~ condition, target_mapping = t2g)
so.pfc <- sleuth_fit(so.pfc)
so.pfc <- sleuth_wt(so.pfc, "conditionhemi_PFC")
so.pfc <- sleuth_fit(so.pfc, ~1, "reduced")

so.pfc <- sleuth_prep(s2c.pfc, ~ condition, target_mapping = t2g)
so.pfc <- sleuth_fit(so.pfc)
so.pfc <- sleuth_wt(so.pfc, "conditionhemi_PFC")
so.pfc <- sleuth_fit(so.pfc, ~1, "reduced")

so.pfc <- sleuth_prep(s2c.pfc, ~ condition, target_mapping = t2g)
so.pfc <- sleuth_fit(so.pfc)
so.pfc <- sleuth_wt(so.pfc, "conditionhemi_PFC")
so.pfc <- sleuth_fit(so.pfc, ~1, "reduced")





