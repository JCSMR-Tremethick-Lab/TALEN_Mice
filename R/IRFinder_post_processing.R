# IRFinder post-processing
library(DESeq2)
IRFinderDir <- "/home/sebastian/miniconda3/envs/irfinder/opt/irfinder-1.2.3/bin"
source(paste(IRFinderDir, "DESeq2Constructor.R", sep = "/")) 

filePaths <- list.files("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/processed_data/GRCm38_ensembl84/IRFinder_analysis/experiment", 
                        pattern = "IRFinder-IR-dir.txt",
                        recursive = T, 
                        full.names = T)

filePaths <- filePaths[grep("pooled", filePaths, invert = T)]

experiment <- data.frame(SampleNames = unlist(lapply(strsplit(filePaths, "/"), function(x) x[14])),
                         condition = factor(unlist(lapply(strsplit(filePaths, "/"), function(x) x[13])), levels = c("wt", "hemi")))

metaList <- DESeqDataSetFromIRFinder(filePaths=filePaths, designMatrix=experiment, designFormula=~1)
dds <- metaList$DESeq2Object
design(dds) <- ~condition + condition:IRFinder
dds <- DESeq(dds)
resultsNames(dds)

res.WT <- results(dds, name = "conditionwt.IRFinderIR")

testIR.WT <- results(dds, name = "conditionwt.IRFinderIR")
# This tests if the number of IR reads are significantly different from the sum of normal spliced reads and intronic reads, in the WT samples.
# We might only be interested in the "log2FoldChange" column, instead of the significance.
# This is because "log2FoldChange" represents log2(number of intronic reads/the sum of normal spliced reads and intronic reads).
# It actually equals log2(IR ratio) in the WT samples.
# So the IR ratio for each intron in WT samples can be easily extracted by the following line
IRratio.WT <- 2^testIR.WT$log2FoldChange    

res.diff <- results(dds, contrast=list("conditionhemi.IRFinderIR","conditionwt.IRFinderIR"))     
# This actually returns the Wald test result of each intron for differential IR analysis.
# This is because it is testing for if the two log2-transformed fold changes are different from each other of not.
# It equals to test if the log2(IRratio.WT) is the same as log2(IRratio.KO), which actually compares the IR ratio in two conditions.

WT.IR_vs_Splice <- 2^res.WT$log2FoldChange
IRratio.WT = WT.IR_vs_Splice/(1+WT.IR_vs_Splice)

res.HEMI = results(dds, name = "conditionhemi.IRFinderIR")
HEMI.IR_vs_Splice=2^res.HEMI$log2FoldChange
IRratio.HEMI = HEMI.IR_vs_Splice/(1+HEMI.IR_vs_Splice)
res.diff = results(dds, contrast=list("conditionhemi.IRFinderIR","conditionwt.IRFinderIR"))  

IR.change = IRratio.HEMI - IRratio.WT
plot(IR.change,col=ifelse(res.diff$padj < 0.05 & abs(IR.change)>=0.1, "red", "black"))

# inspect raw data --------------------------------------------------------
rawData <- lapply(filePaths, function(x) {data.table::fread(x)})
names(rawData) <- unlist(lapply(strsplit(filePaths, "/"), function(x) x[14]))


# results from pooled replicates analysis ---------------------------------
pooledData <- data.table::fread("/home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_testes/processed_data/GRCm38_ensembl84/IRFinder_analysis/experiment/wt_vs_hemi.tab", header = T)

