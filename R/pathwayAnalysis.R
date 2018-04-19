library(clusterProfiler)
library(org.Mm.eg.db)

# pathway/GO analysis on FC DE genes
mouse <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = ensemblHost)
human <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = ensemblHost)

attribs <- biomaRt::listAttributes(mouse)

fcGenes <- unique(sr[qval < 0.1]$ens_gene)

fcGenesHsap <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), 
                              filters = "ensembl_gene_id",
                              mart = mouse,
                              values = fcGenes)

fcGenesHsapEntreIDs <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene"), 
                                      filters = "ensembl_gene_id",
                                      mart = human,
                                      values = fcGenesHsap)

Universe <- biomaRt::getBM(attributes = "entrezgene", 
                           mart = human)

#
minGSSize <- 2
qvalueCutoff <- 0.2
OrgDb <- org.Mm.eg.db

egoMF <- enrichGO(gene          = fcGenesHsapEntreIDs$entrezgene[!is.na(fcGenesHsapEntreIDs$entrezgene)],
                  universe      = Universe,
                  OrgDb         = OrgDb,
                  keytype       = 'ENTREZID',
                  ont           = "MF",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = qvalueCutoff,
                  minGSSize     = minGSSize)

egoBP <- enrichGO(gene          = fcGenesHsapEntreIDs$entrezgene[!is.na(fcGenesHsapEntreIDs$entrezgene)],
                  universe      = Universe,
                  OrgDb         = OrgDb,
                  keytype       = 'ENTREZID',
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = qvalueCutoff,
                  minGSSize     = minGSSize)

egoCC <- enrichGO(gene          = fcGenesHsapEntreIDs$entrezgene[!is.na(fcGenesHsapEntreIDs$entrezgene)],
                  universe      = Universe,
                  OrgDb         = OrgDb,
                  keytype       = 'ENTREZID',
                  ont           = "CC",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = qvalueCutoff,
                  minGSSize     = minGSSize)

gmtfiles <- list.files("~/Data/References/Annotations/MSigDB/msigdb_v6.0_GMTs", pattern = "symbols", full.names = T)
gmtfiles <- gmtfiles[grep("all", gmtfiles)]
gmtfilesNames <- unlist(lapply(strsplit(gmtfiles, "/"), function(x) x[9]))

gmtList <- lapply(gmtfiles, function(y){
  gmt <- read.gmt(y)
  egmt <- GSEA(geneList = geneList, 
               TERM2GENE=gmt, 
               nPerm = 1000, 
               exponent = 1, 
               seed = T, 
               pAdjustMethod = "fdr")
})
names(gmtList) <- gmtfilesNames
