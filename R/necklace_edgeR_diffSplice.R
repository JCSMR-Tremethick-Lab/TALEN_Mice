library(data.table)
library(edgeR)

pathPrefix <- "~"
if (!dir.exists(file.path(pathPrefix, "Data/Tremethick/Hodgkins-Lymphoma/RNA-Seq/necklace_results/R_analysis"))) {
  dir.create("R_analysis")
  setwd(file.path(pathPrefix, "Data/Tremethick/Hodgkins-Lymphoma/RNA-Seq/necklace_results/R_analysis"))
} else {
  setwd(file.path(pathPrefix, "Data/Tremethick/H
                  
                  odgkins-Lymphoma/RNA-Seq/necklace_results/R_analysis"))
}
runConfig <- jsonlite::fromJSON(file.path(pathPrefix, "Development/JCSMR-Tremethick-Lab/Hodgkins-Lymphoma/snakemake/configs/config_RNA-Seq.json"))

dataDir <- "/home/sebastian/Data/Tremethick/Hodgkins-Lymphoma/RNA-Seq/necklace_results/counts"
blockCounts <- data.table::fread(file.path(dataDir, "HodgLymph_block.counts"))
blockCounts[, grep("L1236", colnames(blockCounts)):=NULL]
colnames(blockCounts)[7:18] <- unlist(lapply(strsplit(unlist(lapply(strsplit(colnames(blockCounts)[7:18], "\\."), function(x) x[1])), "/"), function(x) x[2]))
colnames(blockCounts)[7:18] <- names(runConfig$samples$`RNA-Seq`$NB501086_0143_TSoboleva_JCSMR_RNAseq1_12)[sapply(colnames(blockCounts)[7:18], function(cond) grep(cond, runConfig$samples$`RNA-Seq`$NB501086_0143_TSoboleva_JCSMR_RNAseq1_12)) ]

# prepare DGEList object for edgeR input
y <- edgeR::DGEList(counts = blockCounts[, c(7:18)], genes = blockCounts[, c(1:6)])
y$samples
condition <- as.factor(unlist(lapply(strsplit(colnames(y), "-"), function(x) x[2])))
condition <- relevel(condition, ref = "NTC")
timepoint <- as.factor(unlist(lapply(strsplit(colnames(y), "_"), function(x) x[2])))
timepoint <- relevel(timepoint, ref = "48h")
design <- model.matrix( ~ condition + timepoint)

keep <- rowSums(cpm(y) > 1) >=3
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
plotMDS(y)

y <- edgeR::estimateDisp(y, design, robust=TRUE)
y$common.dispersion
edgeR::plotBCV(y)

fit <- edgeR::glmQLFit(y, design, robust=TRUE)
edgeR::plotQLDisp(fit)

qlf <- edgeR::glmQLFTest(fit, coef=2)
topTags(qlf)
is.de <- decideTests(qlf, p.value=0.05)
summary(is.de)
sp <- diffSpliceDGE(fit, coef=2, geneid="Geneid", exonid="Start")
topSpliceDGE(sp, test="gene", n=20)
topSpliceDGE(sp, test="Simes", n=20)

topSpliceDGE(sp, test="gene", n=20)[4,]$Geneid

par(mfrow=c(1,2))
plotSpliceDGE(sp, geneid="ENSG00000113360:MSTRG.40096", genecol="Geneid") 
plotSpliceDGE(sp, geneid="ENSG00000141068:ENSG00000141068:ENSG00000266728:MSTRG.21482", genecol="Geneid") 
plotSpliceDGE(sp, geneid=topSpliceDGE(sp, test="gene", n=20)[4,]$Geneid, genecol="Geneid") 


