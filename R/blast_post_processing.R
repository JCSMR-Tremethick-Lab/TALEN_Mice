require(data.table)
require(ggplot2)
blastResults <- data.table::fread("~/Data/Tremethick/TALENs/RNA-Seq/necklace_results/testis_rs/superT/novelTranscriptCandidates_vs_ercc.blast")
erccInfo <- data.table::fread("/Data/References/Transcriptomes/ERCC/cms_095046.txt")
erccInfo[order(erccInfo$`concentration in Mix 1 (attomoles/ul)`)]
erccInfo <- erccInfo[order(erccInfo$`concentration in Mix 1 (attomoles/ul)`, decreasing = T)]
setkey(blastResults, "V2")
setkey(erccInfo, "ERCC ID")

ggplot(erccInfo[blastResults][order(erccInfo[blastResults]$`concentration in Mix 1 (attomoles/ul)`, decreasing = T)], 
       aes(x = 1:27, y = log10(`concentration in Mix 1 (attomoles/ul)`))) + 
  geom_line() + 
  ggtitle("Concentrations observed [log10(attomoles/uL)]")


ggplot(erccInfo[order(erccInfo$`concentration in Mix 1 (attomoles/ul)`, decreasing = T)], 
       aes(x = 1:nrow(erccInfo), y = log10(`concentration in Mix 1 (attomoles/ul)`))) + 
  geom_line() + 
  ggtitle("Concentrations expected [log10(attomoles/uL)]")

spikeIns <- blastResults$V1
length(unique(erccInfo[blastResults]$`ERCC ID`))
