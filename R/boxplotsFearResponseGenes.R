BDNF - ENSMUSG00000048482
Arc - ENSMUSG00000022602
Fos - ENSMUSG00000021250
Egr1/ZIF268 - ENSMUSG00000038418
Egr2 - ENSMUSG00000037868
Npas4 - ENSMUSG00000045903
Ppp3ca - ENSMUSG00000028161
Ppp1cc - ENSMUSG00000004455

H2afb3 - ENSMUSG00000083616
H2afb1/Gm14920 - ENSMUSG00000067441

c("ENSMUSG00000048482", 'ENSMUSG00000022602', "ENSMUSG00000021250", "ENSMUSG00000038418", "ENSMUSG00000037868", "ENSMUSG00000045903", "ENSMUSG00000028161", "ENSMUSG00000004455")

ggplot(kt.gene[c("ENSMUSG00000048482", 'ENSMUSG00000022602', "ENSMUSG00000021250", "ENSMUSG00000038418", "ENSMUSG00000037868", "ENSMUSG00000045903", "ENSMUSG00000028161", "ENSMUSG00000004455")], aes(groups, value)) + geom_boxplot() + facet_wrap("ens_gene")

ggplot(kt.gene["ENSMUSG00000021250"], aes(groups, value)) + geom_boxplot() + labs(title = "Fos - ENSMUSG00000021250")
ggplot(kt.gene["ENSMUSG00000004455"], aes(groups, value)) + geom_boxplot() + labs(title = "Ppp1cc - ENSMUSG00000004455")
ggplot(kt.gene["ENSMUSG00000028161"], aes(groups, value)) + geom_boxplot() + labs(title = "Ppp3ca - ENSMUSG00000028161")
ggplot(kt.gene["ENSMUSG00000038418"], aes(groups, value)) + geom_boxplot() + labs(title = "Egr1/ZIF268 - ENSMUSG00000038418")
ggplot(kt.gene["ENSMUSG00000037868"], aes(groups, value)) + geom_boxplot() + labs(title = "Egr2 - ENSMUSG00000037868")
ggplot(kt.gene["ENSMUSG00000083616"], aes(groups, value)) + geom_boxplot() + labs(title = "H2afb3 - ENSMUSG00000083616")
ggplot(kt.gene["ENSMUSG00000067441"], aes(groups, value)) + geom_boxplot() + labs(title = "H2afb1/Gm14920 - ENSMUSG00000067441")
