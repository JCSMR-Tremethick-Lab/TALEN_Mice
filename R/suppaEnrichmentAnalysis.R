library(clusterProfiler)
library(org.Mm.eg.db)
library(data.table)

setwd("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/R_Analysis/GRCm38_ensembl93_ERCC/NB501086_0219_TSoboleva_JCSMR_RNAseq")
DTU <- data.table::fread("DifferentialTranscriptUsage_suppa_results.csv")
DGE <- data.table::fread("edgeR_DGE_results.csv")


# Gene ontology analysis of SUPPA data --------------------------------------
keyType <- "SYMBOL"
ontologies <- c("MF", "BP", "CC")
pAdjust <- "fdr"
pCutoff <- 0.10
analysis <- "clusterProfiler"

geneUniverse <- DGE$external_gene_name

suppaFiles <- list.files(pattern = "suppa")
suppaFiles <- suppaFiles[grep("intersect", suppaFiles, invert = T)]
suppaGO <- lapply(suppaFiles, function(x){
  tab <- data.table::fread(x)
  tab <- merge(tab, 
               unique(t2g[,c("ens_gene", "ext_gene")]), 
               by.x = "ensembl_gene_id", 
               by.y = "ens_gene", 
               all.x = T, 
               all.y = F)
  geneList <- unique(tab$ext_gene)
  MF <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Mm.eg.db,
                 keyType = keyType,
                 ont = "MF",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  BP <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Mm.eg.db,
                 keyType = keyType,
                 ont = "BP",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  CC <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Mm.eg.db,
                 keyType = keyType,
                 ont = "CC",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  return(list(MF, BP, CC))
})
names(suppaGO) <- suppaFiles
dotplot(suppaGO$MX_suppa_results.csv[[3]])
tab <- fread(suppaFiles[6])
tab <- tab[, -c("V1", "qval", "lfdr")]
tab <- merge(tab, 
             unique(t2g[,c("ens_gene", "ext_gene")]), 
             by.x = "ensembl_gene_id", 
             by.y = "ens_gene", 
             all.x = T, 
             all.y = F)
save(tab, file = "MX_suppa_results.csv")
tab1 <- as.data.table(suppaGO$MX_suppa_results.csv[[3]])
save(tab1, file = "MX_suppa_results_GO_CC.csv")


# mapping of mouse to human gene IDs to facilitate GSEA analysis ----------
library(biomaRt)
dataset <- "mmusculus_gene_ensembl"
biomart <- "ensembl"
ensemblHost <- "www.ensembl.org"

mart <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, host = ensemblHost)
attribs <- biomaRt::listAttributes(mart)
ensGenes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                          "hsapiens_homolog_ensembl_gene",
                                          "hsapiens_homolog_associated_gene_name"), 
                           values = DGE$ensembl_gene_id,
                           mart = mart)
ensGenes <- data.table::as.data.table(ensGenes)
table(duplicated(ensGenes[!hsapiens_homolog_ensembl_gene == ""]$ensembl_gene_id))

merge(DGE, ensGenes[!hsapiens_homolog_ensembl_gene == ""], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
mappedDGE <- merge(DGE, ensGenes[!hsapiens_homolog_ensembl_gene == ""], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
rankedList <- mappedDGE[order(table.logCPM, decreasing = T)]$table.logCPM
rankedList <- mappedDGE[order(table.logCPM, decreasing = T)]$hsapiens_homolog_associated_gene_name

hsapDataset <- "hsapiens_gene_ensembl"
mart <- biomaRt::useEnsembl(biomart = biomart, dataset = hsapDataset, host = ensemblHost)

hsapGeneIDs <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                             "entrezgene"), 
                              values = mappedDGE$hsapiens_homolog_ensembl_gene,
                              mart = mart)
hsapGeneIDs <- data.table::data.table(hsapGeneIDs)
hsapGeneIDs <- hsapGeneIDs[!is.na(hsapGeneIDs$entrezgene)]
mappedDGE <- merge(mappedDGE, hsapGeneIDs, by.x = "hsapiens_homolog_ensembl_gene", by.y = "ensembl_gene_id")


rankedList <- mappedDGE[order(table.logFC, decreasing = T)]$table.logFC
names(rankedList) <- mappedDGE[order(table.logFC, decreasing = T)]$entrezgene.y
rankedList <- rankedList[!duplicated(rankedList)]

# GSEA of all expression data ---------------------------------------------
c2.all <- read.gmt(gmtfile)
gmtfile <- "/Data/References/Annotations/MSigDB/msigdb_v6.0_GMTs/c2.cp.v6.0.entrez.gmt"
c2.cp <- read.gmt(gmtfile)

egmt.c2.all <- GSEA(rankedList, TERM2GENE=c2.all)
egmt.c2.cp <- GSEA(rankedList, TERM2GENE = c2.cp, pvalueCutoff = 0.2, pAdjustMethod = "fdr")
egmt.c2.cp.results <- as.data.table(egmt.c2.cp)
egmt.c2.cp.results$ID

gseaplot(egmt.c2.cp, "REACTOME_NONSENSE_MEDIATED_DECAY_ENHANCED_BY_THE_EXON_JUNCTION_COMPLEX")
gse

gmtfile <- "/Data/References/Annotations/MSigDB/msigdb_v6.0_GMTs/msigdb.v6.0.entrez.gmt"
MisgDBV6 <- read.gmt(gmtfile)
gsea.MisgDBV6 <- GSEA(rankedList, TERM2GENE = MisgDBV6, pvalueCutoff = 1, pAdjustMethod = "fdr")
gsea.MisgDBV6.results <- as.data.table(gsea.MisgDBV6)
gsea.MisgDBV6.results[order(qvalues)][order(NES, decreasing = T)]
gseaplot(gsea.MisgDBV6, "KEGG_AUTOIMMUNE_THYROID_DISEASE")

# iterate over all individual sets to reduce multiple testing penalty
gmtFilesList <- list.files("/Data/References/Annotations/MSigDB/msigdb_v6.0_GMTs", pattern = "entrez", full.names = T)
# remove the full set
gmtFilesList <- gmtFilesList[-21]

pvalueCutoff <- 0.1
gseaList <- lapply(gmtFilesList, function(gmtfile){
  set <- read.gmt(gmtfile = gmtfile)
  gsea <- GSEA(rankedList, TERM2GENE=set, pAdjustMethod = "fdr", pvalueCutoff = pvalueCutoff)
  return(gsea)
})

gmtFilesNames <- list.files("/Data/References/Annotations/MSigDB/msigdb_v6.0_GMTs", pattern = "entrez", full.names = F)
gmtFilesNames <- gmtFilesNames[-21]
names(gseaList) <- gmtFilesNames
lapply(gseaList, function(x) {table(as.data.table(x)$qvalues < 0.1)})

# running enrichment analysis on MSigDB -----------------------------------
library(snowfall)
sfStop()
sfInit(cpus = 8, parallel = T)
sfExport(list = c("suppaFiles", "t2g", "hsapGeneIDs", "ensGenes", "gmtFilesList", "pvalueCutoff", "gmtFilesNames"))
sfLibrary(clusterProfiler)
sfLibrary(data.table)
sfLibrary(org.Mm.eg.db)

suppaDataList <- lapply(suppaFiles, function(x){
  print(x)
  tab <- data.table::fread(x)
  tab <- merge(tab, 
               unique(t2g[,c("ens_gene", "ext_gene")]), 
               by.x = "ensembl_gene_id", 
               by.y = "ens_gene", 
               all.x = T, 
               all.y = F)
  geneList <- unique(tab$ensembl_gene_id)
  geneListHsapEntrez <- hsapGeneIDs[ensGenes[geneList][!hsapiens_homolog_ensembl_gene == ""]$hsapiens_homolog_ensembl_gene]$entrezgene
  return(list(tab = tab, geneListHsapEntrez = geneListHsapEntrez))
})

sfStop()
sfInit(cpus = 8, parallel = T)
sfExport(list = c("suppaDataList", "gmtFilesList", "pvalueCutoff", "gmtFilesNames"))
sfLibrary(clusterProfiler)
sfLibrary(data.table)

suppaEnricher <- sfLapply(suppaDataList, function(x){
  geneListHsapEntrez <- x$geneListHsapEntrez
  enrichList <- lapply(gmtFilesList, function(gmtfile){
    set <- read.gmt(gmtfile = gmtfile)
    enrich <- enricher(geneListHsapEntrez, TERM2GENE=set, pAdjustMethod = "fdr", pvalueCutoff = pvalueCutoff)
    return(enrich)
  })
  names(enrichList) <- gmtFilesNames
  return(enrichList)
})
sfStop()
names(suppaEnricher) <- unlist(lapply(strsplit(suppaFiles, "\\."), function(x) x[1]))

lapply(suppaEnricher, function(x) {
  lapply(x, function(y){
    table(as.data.frame(y)$qvalue < 0.1)
  })
})

lapply(names(suppaEnricher), function(x) {
  lapply(names(suppaEnricher[[x]]), function(y){
    df <- as.data.frame(suppaEnricher[[x]][[y]])
    if (nrow(df) > 0){
      fn <- paste(x,y,sep="_")
      fn <- paste(fn,"csv",sep=".")
      entrezIds <- unlist(strsplit(df$geneID, "/"))
      df$gene_name <- paste(unique(t2g[ens_gene %in% ensGenes[hsapiens_homolog_ensembl_gene %in% hsapGeneIDs[entrezgene %in% entrezIds]$ensembl_gene_id]$ensembl_gene_id, c("ens_gene", "ext_gene")]$ext_gene), collapse = "/")
      write.csv(df, file = fn)
    }
  })
})

