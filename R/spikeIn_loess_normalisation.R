# ERCC based normalisation as described in Loven et al. 2012
# use the abundance data as imported from kallisto using tximport
mat1 <- txi$abundance

# filter out lowly expressed genes [sum(tpm) < 1]
filter <- apply(mat1, 1, function(x) length(x[x > 1])>=3)
table(filter)
mat1 <- mat1[filter,]
filter <- apply(mat1, 1, function(x) length(x[x > 6000])>=3)
table(filter)
mat1 <- mat1[-filter,]
spikes <- grep("ERCC", rownames(mat1))
normMat1 <- normalize.loess(mat1, subset = spikes, log.it = F)

ensGenes[grep("H2afb3", external_gene_name)]

normMat1["ENSMUSG00000083616",]
mva.pairs(normMat1, log.it = F)
