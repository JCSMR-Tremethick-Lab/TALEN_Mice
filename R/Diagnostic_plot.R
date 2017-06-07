ERCCs <- sleuthProcessedDataCompressed[["HIPPO"]][["kallisto_table"]][grep("ERCC", sleuthProcessedDataCompressed[["HIPPO"]][["kallisto_table"]]$target_id),]
ERCCs <- rbind(ERCCs, sleuthProcessedDataCompressed[["PFC"]][["kallisto_table"]][grep("ERCC", sleuthProcessedDataCompressed[["PFC"]][["kallisto_table"]]$target_id),])
ERCCs <- rbind(ERCCs, sleuthProcessedDataCompressed[["OB"]][["kallisto_table"]][grep("ERCC", sleuthProcessedDataCompressed[["OB"]][["kallisto_table"]]$target_id),])
Transcripts <- sleuthProcessedDataCompressed[["HIPPO"]][["kallisto_table"]][grep("ENS", sleuthProcessedDataCompressed[["HIPPO"]][["kallisto_table"]]$target_id),]
Transcripts <- rbind(Transcripts, sleuthProcessedDataCompressed[["PFC"]][["kallisto_table"]][grep("ENS", sleuthProcessedDataCompressed[["PFC"]][["kallisto_table"]]$target_id),])
Transcripts <- rbind(Transcripts, sleuthProcessedDataCompressed[["OB"]][["kallisto_table"]][grep("ENS", sleuthProcessedDataCompressed[["OB"]][["kallisto_table"]]$target_id),])

pdf("Diagnostic_boxplots.pdf", paper = "a4r")
p1 <- ggplot(data = ERCCs, mapping = aes(x = sample, y = log2(tpm + 1), fill = condition))
p1 + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "ERCC spike ins") # looks like spike ins have massive variability
p2 <- ggplot(data = Transcripts, mapping = aes(x = sample, y = log2(tpm + 1), fill = condition))
p2 + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "All transcripts")
dev.off()
