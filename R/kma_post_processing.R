# intron retention analysis
# RNA-Seq data was aligned using bowtie2 against a custom index of transcripts and introns
# quantification was performed using eXpress

library(kma)

intron_to_trans <- data.table::fread(file.path("/home/sebastian/Data/References/Transcriptomes/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.introns/",
                                               "intron_to_transcripts.txt"), data.table = FALSE)
head(intron_to_trans)

base_dir <- "/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/processed_data/GRCm38_ensembl84/KMA_analysis/"
xprs_fnames <- Sys.glob(file.path(base_dir, "experiment/*/*/express/results.xprs"))

sample_names <- sub(file.path(base_dir, "experiment/[a-z]+/"), "", xprs_fnames) %>%
  sub("express/results.xprs", "", .) %>%
  gsub("/", "", .)
sample_names

condition_names <- sub("[0-9]+", "", sample_names)
condition_names

xprs <- read_express(xprs_fnames, sample_names, condition_names)
names(xprs)

ir <- debugonce(newIntronRetention(xprs$tpm, intron_to_trans, xprs$condition,
                         xprs$uniq_counts))

ir <- ir %>%
  filter_low_tpm(1) %>%
  filter_perfect_psi() %>%
  filter_low_frags(3)