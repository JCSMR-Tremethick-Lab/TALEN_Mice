require(rtracklayer)
require(GenomicFeatures)
require(gdata)

mmusEnsemblTxDB <- GenomicFeatures::makeTxDbFromBiomart(dataset = "mmusculus_gene_ensembl")
chromInfo <- GenomicFeatures::getChromInfoFromBiomart(dataset="mmusculus_gene_ensembl")

# create BED file for all genes -------------------------------------------
grGenes <- genes(mmusEnsemblTxDB)
deepToolsUtils::WriteGRangesToBED(grGenes, out_file = "/Data/References/Annotations/Mus_musculus/GRCh38_ensembl93/allGenes.bed")
exons <- GenomicFeatures::exonsBy(mmusEnsemblTxDB, by = "gene")
grExons <- unlist(exons)
grExons <- reduce(grExons)
deepToolsUtils::WriteGRangesToBED(grExons, out_file = "/Data/References/Annotations/Mus_musculus/GRCh38_ensembl93/allExons.bed")
grIntrons <- intronsByTranscript(mmusEnsemblTxDB, use.names = T)
grIntrons <- unlist(grIntrons)
grIntrons <- reduce(grIntrons)
deepToolsUtils::WriteGRangesToBED(grIntrons, out_file = "/Data/References/Annotations/Mus_musculus/GRCh38_ensembl93/allIntrons.bed")
