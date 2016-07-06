require(stringr)
require(GenomicRanges)
require(rtracklayer)
require(biomaRt)
require(TFBSTools)
require(BSgenome.Mmusculus.UCSC.mm10)
require(JASPAR2014)

genome <- BSgenome.Mmusculus.UCSC.mm10
seqlevels(genome) <- gsub("chr", "", seqlevels(genome))


fantom5 <- read.table("~/Downloads/mm9.cage_peak_phase1and2combined_tpm_ann_decoded.osc.txt.gz.extract.tsv", header = T, as.is = T, sep = "\t")
fantom5$chrom <- unlist(lapply(str_split(fantom5$X00Annotation, ":"), function(x) x[1]))
fantom5$strand <- unlist(lapply(str_split(fantom5$X00Annotation, ","), function(x) x[2]))
fantom5$start <- as.integer(unlist(lapply(str_split(unlist(lapply(str_split(fantom5$X00Annotation, ":"), function(x) x[2])), "\\."), function(y) y[1])))
fantom5$end <- as.integer(unlist(lapply(str_split(unlist(lapply(str_split(unlist(lapply(str_split(fantom5$X00Annotation, ":"), function(x) x[2])), "\\."), function(y) y[3])), ","), function(z) z[1])))
fantom5 <- fantom5[,-1]
gr.fantom5 <- GRanges(fantom5$chrom,
                      IRanges(fantom5$start, fantom5$end),
                      strand = fantom5$strand,
                      mcols = fantom5[,c(1:6)])

chain <- import.chain("/Users/u1001407/Data/References/Annotations/Mus_musculus/mm9/UCSC/mm9ToMm10.over.chain")
gr.fantom5.mm10 <- unlist(liftOver(gr.fantom5, chain))
seqlevels(gr.fantom5.mm10) <- gsub("chr", "", seqlevels(gr.fantom5.mm10))
save(gr.fantom5.mm10, file = "~/Data/References/Annotations/Mus_musculus/GRCm38_mm10/FANTOM5/gr.fantom5.mm10.rda")

mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2013.archive.ensembl.org")
ensGenes <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_id"), mart = mouse)

toGrep <- c("H2afb3", "H2afb2", "Gm14920")
toGrep <- ensGenes[grep(paste(toGrep, collapse = "|"), ensGenes$external_gene_id),]$ensembl_gene_id

ensGenes <- getBM(attributes = c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
                  filter = "ensembl_gene_id",
                  values = toGrep,
                  mart = mouse)

ensGenes <- GRanges(ensGenes$chromosome_name,
                    IRanges(start = ensGenes$start_position,
                            end = ensGenes$end_position,
                            names = ensGenes$ensembl_gene_id),
                    strand = c("+", "-")[match(ensGenes$strand, c(1, -1))],
                    external_gene_id = ensGenes$external_gene_id)

db <- file.path(system.file("extdata", package = "JASPAR2014"), "JASPAR2014.sqlite")
opts <- list()
opts[["collection"]] <- "CORE"
opts[["matrixtype"]] <- "PWM"
ms1 <- getMatrixSet(db, opts)

set1 <- getSeq(genome, promoters(ensGenes, upstream = 1500, downstream = 500))
set1[1]

searchSetList <- searchSeq(ms1, set1[[1]], seqname = "ENSMUSG00000067441", min.score="80%", strand="*")
gff3 <- writeGFF3(searchSetList)
pvals <- pvalues(searchSetList, type = "sampling")





