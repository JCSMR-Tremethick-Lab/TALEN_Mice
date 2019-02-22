library(RUVSeq)
library(tximport)
library(data.table)
library(EDASeq)
library(RColorBrewer)
library(ggplot2)
library(sleuth)
library(GenomicRanges)

# comments: remove KO_19_26 & WT_46_47
# libPaths "/home/sebastian/miniconda3/envs/r_351/lib/R/library"

# load config.json
runInfo <- jsonlite::fromJSON("~/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/configs/config_RNA-Seq_round_spermatids.json")
runInfo$samples$`RNA-Seq`$NB501086_0219_TSoboleva_JCSMR_RNAseq

# import kallisto data using tximport -------------------------------------
setwd("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/R_analysis/")
base_dir <- "/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/kallisto/quant/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq"
load("t2g.rda")

ERCCs <- rtracklayer::import("~/Data/References/Transcriptomes/ERCC/ERCC92.gtf")
t2g <- rbind(t2g, data.table(target_id = as.character(seqnames(ERCCs)),
                             ens_gene = ERCCs$gene_id, 
                             ext_gene = ERCCs$gene_id, 
                             description = ERCCs$gene_id,
                             version = 1,
                             transcript_version = 1))

sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) x[1]))
names(condition) <- sample_id
files <- paste(kal_dirs, "abundance.h5", sep = "/")
names(files) <- sample_id
#files <- files[-which(names(files) %in% c("KO_19_26", "WT_46_47"))]
#sample_id <- sample_id[-c(1,6)]

s2c <- data.table::data.table(sample = sample_id, condition = condition)
s2c$path <- kal_dirs
s2c$sample <- as.character(s2c$sample)
s2c$condition <- factor(s2c$condition)
s2c$condition <- relevel(s2c$condition, ref = "WT")
design <- model.matrix(~ condition, data = s2c)

#files[7] <- "/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_brain_experiment_2/processed_data/GRCm38_ensembl84/kallisto/mmus_pfc_10_51_mut_FC/abundance.h5"

txi <- tximport::tximport(files, 
                          type = "kallisto",
                          geneIdCol = gene_id,
                          txIdCol = target_id,
                          tx2gene = t2g[,c(1,2)],
                          ignoreTxVersion = F,
                          countsFromAbundance = "no")
geneCounts <- txi$counts
geneCounts <- round(geneCounts, digits = 0)
geneLengths <- txi$length
dim(geneCounts)

filter <- apply(geneCounts, 1, function(x) length(x[x > 1])>=1)
table(filter)

filtered <- geneCounts[filter,]
head(filtered)
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
head(spikes)

x <- s2c$condition
set <- EDASeq::newSeqExpressionSet(filtered, phenoData = data.frame(x, row.names = colnames(filtered)))

# exploratory diagnostics -------------------------------------------------
colors <- RColorBrewer::brewer.pal(3, "Set2")
EDASeq::plotRLE(set, outline=FALSE, ylim=c(-1, 1), col=colors[x])
EDASeq::plotPCA(set, col=colors[x], cex=1.2)

set <- EDASeq::betweenLaneNormalization(set, which="upper")
EDASeq::plotRLE(set, outline=FALSE, ylim=c(-1, 1), col=colors[x])
EDASeq::plotPCA(set, col=colors[x], cex=1.2)
# using spike-ins for normalisation ---------------------------------------
set1 <- RUVSeq::RUVg(set, spikes, k=1)
set2 <- RUVSeq::RUVg(set, spikes, k=2)
set3 <- RUVSeq::RUVg(set, spikes, k=3)
set4 <- RUVSeq::RUVg(set, spikes, k=4)
pData(set1)
pData(set2)
pData(set3)
pData(set4)
EDASeq::plotRLE(set1, outline=FALSE, ylim=c(-1, 1), col=colors[x])
EDASeq::plotPCA(set1, col=colors[x], cex=1.2)
EDASeq::plotRLE(set2, outline=FALSE, ylim=c(-1, 1), col=colors[x])
EDASeq::plotPCA(set2, col=colors[x], cex=1.2)
EDASeq::plotRLE(set3, outline=FALSE, ylim=c(-1, 1), col=colors[x])
EDASeq::plotPCA(set3, col=colors[x], cex=1.2)
EDASeq::plotRLE(set4, outline=FALSE, ylim=c(-1, 1), col=colors[x])
EDASeq::plotPCA(set4, col=colors[x], cex=1.2)

w1 <- pData(set1)$W_1
w2 <- pData(set2)$W_2
w3 <- pData(set3)$W_3
names(w1) <- rownames(pData(set1))
names(w2) <- rownames(pData(set2))
names(w3) <- rownames(pData(set3))
# differential gene expression analysis using edgeR -----------------------
library(edgeR)
#design <- model.matrix(~x + W_1:W_2, data=pData(set2))
design <- model.matrix(~x + W_1 + W_2 + W_3, data=pData(set3))
y <- DGEList(counts=counts(set1), group=x)
y <- calcNormFactors(y, method="upper")

# diagnostic plot ---------------------------------------------------------
library(affy)
mat1 <- txi$abundance + 1
spikeRows <- which(rownames(mat1) %in% spikes)
normMat1 <- normalize.loess(mat1, subset = spikeRows, log.it = F)

par(mfrow = c(2,ncol(normMat1)/2))
for (i in (1:ncol(normMat1))){
  plotMD(cpm(y, log = T), column = i)
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

par(mfrow = c(2,ncol(normMat1)/2))
for (i in (1:ncol(normMat1))){
  plotMD(log2(normMat1), column = i)
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

points <- c(0,1,2,15,16,17)
colors <- c("red", "blue")
plotMDS(y, col=colors[x], pch=points[x])
legend("topleft", legend=levels(x), pch=points, col=colors, ncol=2)

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

mva.pairs(cpm(y, log = T), log.it = F)

topTags(lrt)

eResults <- as.data.table(topTags(lrt, n = nrow(lrt$table)))
eResults$ensGene <- rownames(topTags(lrt, n = nrow(lrt$table)))
table(eResults$table.FDR < 0.1)

# edgeR - merge in meta-data ------------------------------------------------------
load("ensGenes.rda")
ensGenes <- as.data.table(ensGenes)
setkey(ensGenes, ensembl_gene_id)
setkey(eResults, ensGene)

eResults <- ensGenes[eResults]
eResults[table.FDR < 0.1]
write.csv(eResults, file = "edgeR_DGE_results_3PCs.csv")
write.csv(eResults[table.FDR <= 0.1], file = "edgeR_DGE_results_FDR0.1_3PCs.csv")
eResults <- data.table::fread("edgeR_DGE_results_3PCs.csv")

eResults[table.FDR < 0.1][table.logFC < 0]


# volcano plots from edgeR results ----------------------------------------
vp1 <- ggplot(data = eResults, aes(x = table.logFC, y = -log10(table.FDR))) +
  geom_point() +
  xlim(-8,8) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_point(data = eResults[table.FDR < 0.1], aes(x = table.logFC, y =-log10(table.FDR), color = "red"))
vp1
eResults[abs(table.logFC) > 1][table.FDR < 0.1]



# Run basic sleuth analysis of the data in order to get normalised --------
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
so <- sleuth::sleuth_fit(so, formula = design)
so <- sleuth::sleuth_wt(so, "conditionKO")

so.gene <- sleuth::sleuth_prep(s2c, ~ condition, target_mapping = t2g, aggregation_column = "ext_gene")
so.gene <- sleuth::sleuth_fit(so.gene, formula = design)
so.gene <- sleuth::sleuth_wt(so.gene, "conditionKO")
sleuth_live(so.gene)

