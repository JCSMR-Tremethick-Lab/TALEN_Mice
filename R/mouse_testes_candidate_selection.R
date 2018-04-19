if (dir.exists(lDir(pathPrefix, "Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/R_analysis"))){
  setwd(lDir(pathPrefix, "Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/R_analysis"))
} else {
  dir.create(lDir(pathPrefix, "Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/R_analysis"))
  setwd(lDir(pathPrefix, "Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/R_analysis"))
}

sleuth_analysis_version <- 1
sleuth_analysis_output <- paste("sleuth_analysis_V", sleuth_analysis_version, "_output.rda", sep = "")
sleuth_analysis_output_combined <- paste("sleuth_analysis_V", sleuth_analysis_version, "_output_combined.rda", sep = "")
sleuth_analysis_outputCompressed <- paste("sleuth_analysis_V", sleuth_analysis_version, "_Compressed_output.rda", sep = "")

ERCCs <- rtracklayer::import("~/Data/References/Transcriptomes/ERCC/ERCC92.gtf")
t2g <- rbind(t2g, data.table(target_id = as.character(seqnames(ERCCs)), ens_gene = ERCCs$gene_id, ext_gene = ERCCs$gene_id, description = ERCCs$gene_id))

rownames(ercc_tab) <- ercc_tab$ERCC.ID
ercc_exp <- kt[grep("ERCC", kt$target_id),]
i1 <- intersect(rownames(ercc_exp[which(ercc_exp$mean > 0),]), rownames(ercc_tab))


kallisto_base_dir <- lDir(pathPrefix, "Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/processed_data/GRCm38_ensembl84_ERCC/kallisto")

sample_id <- dir(kallisto_base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(kallisto_base_dir, id))
condition <- c("hemi", "wt", "wt", "hemi", "wt", "hemi")
condition <- as.factor(condition)
condition <- relevel(condition, ref = "wt")
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
design <- model.matrix(~ condition, data = s2c)

# load kallisto data with tximport and inspect via PCA -------------------------
if (length(kallisto_base_dir) == 1){
  if (dir.exists(kallisto_base_dir)){
    sample_id <- dir(kallisto_base_dir)
    kal_dirs <- sapply(sample_id, function(id) file.path(kallisto_base_dir, id))
    condition <- c("hemi", "wt", "wt", "hemi", "wt", "hemi")
    names(condition) <- sample_id
    files <- paste(kal_dirs, "abundance.h5", sep = "/")
    names(files) <- sample_id
  } else { #end of directory IF
    stop("Directory with kallisto input data is missing!")
  }
} else {
  kal_dirs <- unlist(lapply(runID, function(x) {
    paste(base_dir[grep(x, base_dir)], names(runConfig$samples$`RNA-Seq`[[x]]), sep = "/")
  }))
  sample_id <- unlist(lapply(runID, function(x) {
    names(runConfig$samples$`RNA-Seq`[[x]])
  }))
  condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) paste(x[1:2], collapse = "_")))
  names(kal_dirs) <- sample_id
  names(condition) <- sample_id
  files <- paste(kal_dirs, "abundance.h5", sep = "/")
  names(files) <- sample_id
}

geneIdCol <- "ens_gene"
targetIdCol <- "target_id"

txi <- tximport::tximport(files, 
                          type = "kallisto",
                          tx2gene = subset(t2g, select = c(targetIdCol, geneIdCol)), geneIdCol = geneIdCol, txIdCol = targetIdCol,
                          ignoreTxVersion = F)

# perform PCA for first inspection of data --------------------------------
sd1 <- apply(log2(txi$abundance + 0.5), 1, sd)
summary(sd1)
pca1 <- ade4::dudi.pca(t(log2(txi$abundance[sd1 > 1, ] + 0.5)), scannf = F, nf = 6)
pdf(paste("PCA_TALEN_Mice_testes", annotationVersion, ".pdf", sep = ""))
ade4::s.arrow(pca1$li, boxes = F)
ade4::s.class(pca1$li, fac = as.factor(condition))
dev.off()

selected <- grep("ERCC", rownames(txi$abundance))
pca2 <- ade4::dudi.pca(t(log2(txi$abundance[selected, ] + 0.5)), scannf = F, nf = 6)

## for now, I will remove NMG3-77hemi_S6 and NMG3-76wt_S5
## unfortunately, this results in no DEG...
#s2c <- s2c[!s2c$sample %in% c("NMG3-77hemi_S6", "NMG3-76wt_S5"),]
#design <- model.matrix(~ condition, data = s2c)

so <- sleuth::sleuth_prep(s2c, ~ condition, target_mapping = t2g, extra_bootstrap_summary = T, read_bootstrap_tpm = T)
so <- sleuth::sleuth_fit(so, design)
so <- sleuth::sleuth_wt(so, "conditionhemi")
so <- sleuth::sleuth_fit(so, ~1, "reduced")
so <- sleuth::sleuth_lrt(so, "reduced", "full")
kt <- sleuth::kallisto_table(so, use_filtered = T)
kt <- data.table::data.table(kt)
kt <- merge(kt, subset(t2g, select = c("target_id", "ext_gene")), by.x = "target_id", by.y = "target_id", all.x = T, all.y = F)
kt.gene <- kt[,.(sum_est_counts = sum(est_counts), sum_tpm = sum(tpm)), by = .(ext_gene, sample, condition)]
kt_wide.gene <- tidyr::spread(kt.gene[, c("ext_gene", "sample", "sum_tpm")], sample, sum_tpm)
write.csv(kt_wide.gene, "testes_gene_expression.csv")

setkey(kt.gene, "ext_gene")
subset(kt.gene, ext_gene %in% shortListed$FeatureID)

# loading candidate genes from da Cruz et al. 2016 ------------------------
fn <- "~/Data/Publications/daCruz_etal_2016/12864_2016_2618_MOESM2_ESM.csv"
dT1 <- data.table::fread(fn, verbose = T)
colnames(dT1) <- gsub(" ", "", colnames(dT1))
dT2 <- subset(dT1, `2CRPKM` < 40 & LZRPKM < 100 & PSRPKM > 10 & RSRPKM > 100)
dT2[order(-PSRPKM, RSRPKM)]
setkey(dT2, "FeatureID")
dT2["Cypt12"]
subset(kt.gene, ext_gene %in% dT2$FeatureID)

dT3 <- subset(dT1, `2CRPKM` < 2 & LZRPKM < 2 & PSRPKM > 2 & RSRPKM > 10)
so.results <- sleuth_results(so, "conditionhemi", show_all = F)
so.results <- data.table(so.results)
setkey(so.results, "ext_gene")
so.results.sub <- subset(so.results, ext_gene %in% dT3$FeatureID)
so.results.sub$fdr <- p.adjust(so.results.sub$pval, "fdr")
summary(so.results.sub$fdr)
GOIs <- subset(so.results.sub, pval < 0.1)$ext_gene
hist(so.results.sub$fdr)

# short-listed genes ------------------------------------------------------
# from Tanya's list
fn <- "~/Data/Publications/daCruz_etal_2016/shortListed.csv"
shortListed <- data.table::fread(fn, verbose = T)
colnames(shortListed) <- gsub(" ", "", colnames(shortListed))
kt.gene$log2_sum_tpm <- log2(kt.gene$sum_tpm + 1)
p1 <- ggplot(data = subset(kt.gene, ext_gene %in% shortListed$FeatureID), aes(condition, log2_sum_tpm, colour = condition)) + 
  geom_jitter(width = 0.1) +
  facet_wrap(~ ext_gene, scales = "free_x")
p2 <- ggplot(data = subset(kt.gene, ext_gene %in% GOIs), aes(condition, log2_sum_tpm, colour = condition)) + 
  geom_jitter(width = 0.1) +
  facet_wrap(~ ext_gene, scales = "free_x")
ggsave("shortlisted_expression.pdf", plot = p1, height = 20, width = 20, units = "cm")
