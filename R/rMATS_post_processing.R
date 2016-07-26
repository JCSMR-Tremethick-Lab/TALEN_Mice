# rMATS_post_processing.R

# load external functions
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

getExonStructure <- function(ensembl_gene_id = NULL, mart = NULL){
  if(is.null(ensembl_gene_id)) {stop("Ensembl Gene ID missing")}
  if(is.null(mart)) {stop("BiomaRt object missing")}
  exonStructure <- getBM(attributes = c("ensembl_gene_id",
                                        "ensembl_exon_id",
                                        "chromosome_name",
                                        "exon_chrom_start",
                                        "exon_chrom_end",
                                        "is_constitutive",
                                        "rank"),
                         values = ensembl_gene_id,
                         filter = "ensembl_gene_id",
                         mart = mart)
  return(exonStructure)
}

if (amILocal("JCSMR027564ML")){
  setwd("~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  dataPath <- "~/mount/gduserv/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/rMATS/"
} else {
  setwd("~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/R_analysis")
  dataPath <- "~/Data/Tremethick/TALENs/NB501086_0063_TSoboleva_JCSMR_standed_RNAseq/processed_data/GRCm38_ensembl84_ERCC/rMATS/"
}

if (!file.exists("ensGenes.rda")){
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
  save(ensGenes, file = "ensGenes.rda")
  
  ensTranscripts <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                  "ensembl_transcript_id",
                                                  "transcript_length"),
                                   mart = mouse,
                                   filter = "ensembl_gene_id",
                                   values = ensGenes$ensembl_gene_id)
  save(ensTranscripts, file = "ensTranscripts.rda")
  
  mylength <- sapply(ensGenes$ensembl_gene_id, function(x){
    y <- ensTranscripts[which(ensTranscripts$ensembl_gene_id == x), ]
    y <- y[which.max(y$transcript_length), ]$transcript_length
  })
  save(mylength, file = "mylength.rda")
  
  mygc <- ensGenes$percentage_gc_content
  names(mygc) <- ensGenes$ensembl_gene_id
  save(mygc, file = "mygc.rda")
  
  mybiotypes <- ensGenes$gene_biotype
  names(mybiotypes) <- ensGenes$ensembl_gene_id
  save(mybiotypes, file = "mybiotypes.rda")
  
  mychroms <- data.frame(Chr = ensGenes$chromosome_name, GeneStart = ensGenes$start_position, GeneEnd = ensGenes$end_position)
  rownames(mychroms) <- ensGenes$ensembl_gene_id
  save(mychroms, file = "mychroms.rda")
} else {
  load("ensGenes.rda")
  load("ensTranscripts.rda")
  load("mylength.rda")
  load("mygc.rda")
  load("mybiotypes.rda")
  load("mychroms.rda")
}

tissues <- list.dirs(dataPath, recursive = F, full.names = F)

rMATSResults <- sapply(tissues, function(x){
  files <- list.files(paste(dataPath, x, "WT_vs_HEMI", "MATS_output", sep = "/"), recursive = T, full = T)
  names <- list.files(paste(dataPath, x, "WT_vs_HEMI", "MATS_output", sep = "/"), recursive = T, full = F)
  names <- unlist(lapply(strsplit(names, "\\."), function(x) paste(x[1:3], collapse = ".")))
  results <- sapply(files, function(y){
    tab <- read.table(y, header = T, as.is = T, sep = "\t")
  })
  names(results) <- names
  return(list(results = results))
})

names(rMATSResults) <- tissues

lapply(rMATSResults, function(x){
  lapply(x, function(y) {
    table(y$FDR < 0.2)
  })
})



