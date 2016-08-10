require(NOISeq)
require(biomaRt)

source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

ensemblHost <- "mar2016.archive.asia.ensembl.org"

# setting working directory and data sources ------------------------------
if (amILocal("JCSMR027564ML")){
  setwd("~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/mount/gduserv/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/HTSeq/count/"
} else {
  setwd("~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  load("~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/t2g.rda")
  dataPath <- "~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/HTSeq/count/"
}

# load prepared count data (HTS count output)
load("htSeqCountMatrix.rda")


condition <- as.factor(unlist(lapply(strsplit(colnames(htSeqCountMatrix), "_"), function(z) z[3])))
tissue <- as.factor(unlist(lapply(strsplit(colnames(htSeqCountMatrix), "_"), function(z) z[4])))
myfactors <- data.frame(Tissue = tissue, Condition = condition)

# preparing annotation data from Ensembl
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
                           mart = mouse,
                           filter = "ensembl_gene_id",
                           values = rownames(htSeqCountMatrix))

ensTranscripts <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                "ensembl_transcript_id",
                                                "transcript_length"),
                                 mart = mouse,
                                 filter = "ensembl_gene_id",
                                 values = ensGenes$ensembl_gene_id)

mylength <- sapply(ensGenes$ensembl_gene_id, function(x){
  y <- ensTranscripts[which(ensTranscripts$ensembl_gene_id == x), ]
  y <- y[which.max(y$transcript_length), ]$transcript_length
})

mygc <- ensGenes$percentage_gc_content
names(mygc) <- ensGenes$ensembl_gene_id

mybiotypes <- ensGenes$gene_biotype
names(mybiotypes) <- ensGenes$ensembl_gene_id

mychroms <- data.frame(Chr = ensGenes$chromosome_name, GeneStart = ensGenes$start_position, GeneEnd = ensGenes$end_position)
rownames(mychroms) <- ensGenes$ensembl_gene_id

# create eSet
mydata <- NOISeq::readData(data = htSeqCountMatrix, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms, factors = myfactors, norm)

# actual data exploration
myexplodata <- dat(mydata, type = "biodetection")
explo.plot(myexplodata, plottype = "persample")

mybiodetection <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1,2))
explo.plot(mybiodetection, samples = c(1,2), plottype = "persample")

mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
pdf("NOISeq_analysis_output/Biotype_counts.pdf", paper = "a4")
par(mfrow = c(2,2))
explo.plot(mycountsbio, toplot = 1, samples = 3, plottype = "boxplot")
explo.plot(mycountsbio, toplot = 1, samples = 11, plottype = "boxplot", add = T)
dev.off()

mysaturation <- dat(mydata, k = 0, ndepth = 10, type = "saturation")

# checking biases
# length
mylengthbias <- dat(mydata, factor = "Tissue", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")

# gc
myGCbias <- dat(mydata, factor = "Tissue", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")

# compositional bias
mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 1)




