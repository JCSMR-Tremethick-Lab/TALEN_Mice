# tissue map of transcribed histones in mouse/human,
# with transcribed promoter activity levels
# 1) compile list of all histone genes
# 2) obtain FANTOM5 sequencing data
# 3) obtain matching RNA-Seq data from mouse
# 4) create tissue specific coverage maps for each histone gene

# setting working directory and data sources ------------------------------
setwd("~/mount/gduserv/Data/Tremethick/TALENs/meta_analysis/R_analysis")
base_dir <- "~/mount/gduserv/Data/Tremethick/TALENs/meta_analysis/WT"

# creating data.frame for experimental design and file names --------------
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
condition <- unlist(lapply(strsplit(names(kal_dirs), "_"), function(x) paste(x[3:4], collapse = "_")))
condition[10:12] <- "wt_testis"
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
x <- s2c$condition
x <- factor(x, levels(x)[c(4,1,2,3)])
s2c$condition <- x
