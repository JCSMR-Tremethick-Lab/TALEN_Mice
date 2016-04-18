# map_transcriptome_to_genome_coordinates.R

# required libraries
require(GenomicFeatures)
require(GenomicRanges)

# First get transcriptome annotation
mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2013.archive.ensembl.org")
mmusEnsemblTxDB <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                                       dataset="mmusculus_gene_ensembl",
                                       transcript_ids=NULL,
                                       circ_seqs=DEFAULT_CIRC_SEQS,
                                       filters="",
                                       id_prefix="ensembl_",
                                       host="dec2013.archive.ensembl.org",
                                       port=80,
                                       taxonomyId=NA,
                                       miRBaseBuild=NA)
mmusTranscripts <- transcripts(mmusEnsemblTxDB)
