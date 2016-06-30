require(rhdf5)
require(sleuth)
require(biomaRt)
require(tidyr)

source("~/Development/GeneralPurpose/R/heatmap.3.R")


# setting working directory and data sources ------------------------------
setwd("~/Data/Tremethick/TALENs/NB501086_0047_TSoboleva_JCSMR_stranded_RNASeq/R_analysis")
base_dir <- "~/Data/Tremethick/TALENs/NB501086_0047_TSoboleva_JCSMR_stranded_RNASeq/processed_data/kallisto"

# creating data.frame for experimental design and file names --------------
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
s2c <- data.frame(sample = sample_id, condition = c("KO", "WT", "WT", "KO", "WT", "KO"))
s2c <- dplyr::mutate(s2c, path = kal_dirs)
x <- s2c$condition
x <- factor(x, levels(x)[c(2,1)])
s2c$condition <- x

# transcript level differential expression analysis -----------------------------
#use Ensembl 74 for annotation
mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2013.archive.ensembl.org")
attribs <- listAttributes(mouse)
# annotate transcripts
t2g <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_id"), mart = mouse)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_id)
rownames(t2g) <- t2g$target_id
# re-run sleuth
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_wt(so, "conditionKO")
results_table <- sleuth_results(so, "conditionKO")
results_table$FC_estimated <- log2(exp(results_table$b))


# gene level analysis -----------------------------------------------------
options(mc.cores = 8L)
so.gene <- sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column = 'ens_gene')
so.gene <- sleuth_fit(so.gene)
so.gene <- sleuth_wt(so.gene, "conditionKO")
results_table.gene <- sleuth_results(so.gene, "conditionKO")
results_table.gene$FC_estimated <- log2(exp(results_table.gene$b))
rownames(results_table.gene) <- results_table.gene$target_id

# results table contains gene level differential expression, but kallisto_table is still transcript level
kt.gene <- kallisto_table(so.gene, normalized = T)
kt.gene <- tidyr::spread(kt.gene[, c("target_id", "sample", "tpm")], sample, tpm)

# do gene level summarizing
t <- as.matrix(kt.gene[,c(2:7)])
rownames(t) <- kt.gene$target_id
ensembl_gene_IDs <- as.data.frame(so.gene$target_mapping)
colnames(ensembl_gene_IDs) <- c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol")

require(snowfall)
sfInit( parallel = T, cpus = 8)
sfExport("ensembl_gene_IDs")
l1 <- sfLapply(unique(ensembl_gene_IDs$ensembl_gene_id), function(x) unique(ensembl_gene_IDs[which(ensembl_gene_IDs$ensembl_gene_id == x), "ensembl_transcript_id"])) 
names(l1) <- unique(ensembl_gene_IDs$ensembl_gene_id)
sfExport("l1")
sfExport("t")

l2 <- sfLapply(l1, function(x) {
  r1 <- as.character(unlist(x))
  i1 <- as.character(intersect(rownames(t), r1))
  print(i1)
  if(length(i1) == 1) {
    r <- t[i1, ]
  } else if(length(i1) > 1) {
    r <- apply(t[i1, ], 2, sum)
  }
})
df1 <- t(data.frame(l2[lapply(l2, length)>0]))
sfStop()

# perform PCA
kt.geneLog2 <- log2(df1 + 1)
pca1 <- dudi.pca(t(kt.geneLog2))
samples <- as.factor(c("NMG3-60hemi_S1" = "KO", 
                       "NMG3-62wt_S2" = "WT",
                       "NMG3-74wt_S3" = "WT", 
                       "NMG3-75hemi_S4" = "KO", 
                       "NMG3-76wt_S5" = "WT", 
                       "NMG3-77hemi_S6" = "KO"))

pdf("PCA.pdf", height = 7, width = 7)
s.class(pca1$li, fac = samples)
dev.off()

# integrating H2A.Lap1 dependent gene expression data ---------------------
lap1Dep <- read.table("../../Lap1_dependent_genes/RS_lap_data.txt", 
                      header = F,
                      as.is = T, 
                      sep ="\t")
i1 <- intersect(results_table$ens_gene, lap1Dep$V1)
write.csv(results_table[which(results_table$ens_gene %in% i1),], file = "Lap1_dependent_genes.csv")

# extracting expression levels of sets of target genes --------------------
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



# analysis of Lap1-driven genes (based on mouse testis data) --------------
lap1Genes <- read.table("~/Data/Tremethick/TALENs/Lap1_dependent_genes/RS_lap_data.txt", header = F, as.is = T)
rownames(lap1Genes) <- lap1Genes[,1]
lap1Genes <- lap1Genes[,-1]

lap1Transcripts <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_id"), mart = mouse, filters = "ensembl_gene_id", values = rownames(lap1Genes))
kt.lap1 <- kt[grep(paste(lap1Transcripts$ensembl_transcript_id, collapse = "|"), kt$target_id),]
kt.lap1$sample <- factor(kt.lap1$sample)
kt.lap1 <- tidyr::spread(kt.lap1[, c("target_id", "sample", "tpm")], sample, tpm)
rownames(kt.lap1) <- kt.lap1$target_id
kt.lap1 <- kt.lap1[,-1]

m <- as.matrix(kt.lap1)
ens_genes <- lap1Transcripts
require(snowfall)
sfInit( parallel = T, cpus = 8)
sfExport("ens_genes")
l1 <- sfLapply(unique(ens_genes$ensembl_gene_id), function(x) unique(ens_genes[which(ens_genes$ensembl_gene_id == x), "ensembl_transcript_id"])) 
names(l1) <- unique(ens_genes$ensembl_gene_id)
sfExport("l1")
sfExport("m")

l2 <- sfLapply(l1, function(x) {
  r1 <- as.character(unlist(x))
  i1 <- as.character(intersect(rownames(m), r1))
  #print(i1)
  if(length(i1) == 1) {
    r <- m[i1, ]
  } else if(length(i1) > 1) {
    r <- apply(m[i1, ], 2, sum)
  }
})
df1 <- t(data.frame(l2[lapply(l2, length)>0]))
sfStop()

kt.lap1 <- as.data.frame(df1)
kt.lap1$gene <- unique(lap1Transcripts$external_gene_id)
kt.lap1 <- kt.lap1[, c("NMG3-62wt_S2", "NMG3-74wt_S3", "NMG3-76wt_S5", "NMG3-60hemi_S1", "NMG3-75hemi_S4", "NMG3-77hemi_S6", "gene")]
kt.lap1$target_id <- rownames(kt.lap1)

# add logFC, p-val, etc
kt.lap1 <- merge(kt.lap1, results_table.gene[, c("pval", "qval", "b", "FC_estimated", "target_id")], by.x = "target_id", by.y = "target_id")


# analysis of H2AZ-driven (heterotypic) genes (based on mouse testis data) --------------
hetGenes <- read.table("~/Data/Tremethick/TALENs/Lap1_dependent_genes/RS_het_data.txt", header = F, as.is = T, sep = "\t")
rownames(hetGenes) <- hetGenes[,1]
hetGenes <- hetGenes[,-1]

hetTranscripts <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_id"), mart = mouse, filters = "ensembl_gene_id", values = rownames(hetGenes))
kt.het <- kt[grep(paste(hetTranscripts$ensembl_transcript_id, collapse = "|"), kt$target_id),]
kt.het$sample <- factor(kt.het$sample)
kt.het <- tidyr::spread(kt.het[, c("target_id", "sample", "tpm")], sample, tpm)
rownames(kt.het) <- kt.het$target_id
kt.het <- kt.het[,-1]

m <- as.matrix(kt.het)
ens_genes <- hetTranscripts
require(snowfall)
sfInit( parallel = T, cpus = 8)
sfExport("ens_genes")
l1 <- sfLapply(unique(ens_genes$ensembl_gene_id), function(x) unique(ens_genes[which(ens_genes$ensembl_gene_id == x), "ensembl_transcript_id"])) 
names(l1) <- unique(ens_genes$ensembl_gene_id)
sfExport("l1")
sfExport("m")

l2 <- sfLapply(l1, function(x) {
  r1 <- as.character(unlist(x))
  i1 <- as.character(intersect(rownames(m), r1))
  #print(i1)
  if(length(i1) == 1) {
    r <- m[i1, ]
  } else if(length(i1) > 1) {
    r <- apply(m[i1, ], 2, sum)
  }
})
df1 <- t(data.frame(l2[lapply(l2, length)>0]))
sfStop()

kt.het <- as.data.frame(df1)
kt.het$gene <- unique(hetTranscripts$external_gene_id)
kt.het <- kt.het[, c("NMG3-62wt_S2", "NMG3-74wt_S3", "NMG3-76wt_S5", "NMG3-60hemi_S1", "NMG3-75hemi_S4", "NMG3-77hemi_S6", "gene")]
kt.het$target_id <- rownames(kt.het)

# add logFC, p-val, etc
kt.het <- merge(kt.het, results_table.gene[, c("pval", "qval", "b", "FC_estimated", "target_id")], by.x = "target_id", by.y = "target_id")

# boxplots of lap1/h2az-driven genes in WT vs KO
pdf("Boxplot_Lap1_vs_H2AZ_driven_genes.pdf", height = 5, width = 7.5)
par(mfrow = c(1,2))
boxplot(log2(unlist(kt.lap1[,c(2:4)])), 
        log2(unlist(kt.het[,c(2:4)])),
        names = c("Lap1 genes", "H2AZ genes"),
        main = "WT mice")
boxplot(log2(unlist(kt.lap1[,c(5:7)])),
        log2(unlist(kt.het[,c(5:7)])),
        names = c("Lap1 genes", "H2AZ genes"),
        main = "TALEN KO mice")
dev.off()
        


