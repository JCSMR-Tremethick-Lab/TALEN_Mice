import json
from pprint import pprint
from snakemake.io import expand

REF_VERSION = "blablub"

with open("config_RNA-Seq_round_spermatids.json") as data_file:
    config = json.load(data_file)


expand("{assayType}/kallisto/{reference_version}/{runID}/{library}",
       assayType = "RNA-Seq",
       reference_version = "GRCm38_ensembl84_ERCC",
       runID = "NB501086_0219_TSoboleva_JCSMR_RNAseq",
       library = [y for y in config["samples"]["NB501086_0219_TSoboleva_JCSMR_RNAseq"].keys()])

expand("{assayType}/bowtie2/{reference_version}/{runID}/{library}.{suffix1}",
        assayType = "CutRun",
        reference_version = REF_VERSION,
        runID = "180731_NB501086_0217_CutandRun_Tanya",
        library = ["WT_01_H2AL2_3_7_18", "WT_01_IGG_3_7_18", "KO_01_H2AL2_3_7_18", "KO_02_H2AL2_24_6_18", "WT_02_H2AL2_24_6_18", "WT_01_H3K27me3_23_5_18", "WT_01_H3K36me3_23_5_18"],
        suffix1 = ["bam"],
        suffix2 = ["unmapped_r1.fastq.gz", "unmapped_r2.fastq.gz"])


with open("config_RNA-Seq_Brain.json") as data_file:
    config = json.load(data_file)

wildcards = dict()
wildcards = {"outdir" : "processed_data",
             "reference_version" : "hg38",
             "unit" : "10_9_hemi_OB_CCGTCC",
             "tissue" : "OB",
             "condition" : "WT"}

def getFASTQ(wildcards, assayID):
    fn = []
    for i in wildcards[assayID]:
        fn.append("./" + i)
    return(fn)

def getTPMs(wildcards):
    fn = []
    for i in config["Groups"][wildcards["tissue"]][wildcards["condition"]]:
        fn.append("/".join([wildcards["outdir"], wildcards["reference_version"], "suppa", wildcards["tissue"], wildcards["condition"], i, "abundance.tpm.tsv"]))
    return(",".join(fn))


def writeBAMFilesList(wildcards):
    f = open('workfile', 'w+')
    f.write(wildcards["unit"] +
              "\t" +
              "./" +
              wildcards["outdir"] +
              "/" +
              wildcards["reference_version"] +
              "/STAR/full/" +
              wildcards["unit"] +
              ".aligned.bam" +
              "\tNA\n")
    f.close()



def getGroups(wildcards):
    cond1 = []
    cond2 = []
    c1 = wildcards.condition.split("_vs_")[0]
    c2 = wildcards.condition.split("_vs_")[1]
    for i in config["groups"][wildcards.tissue][c1]:
        cond1.append("./" + wildcards.outdir + "/" + wildcards.reference_version + "/STAR/full/" + i + ".aligned.bam")
    for i in config["groups"][wildcards.tissue][c2]:
        cond2.append("./" + wildcards.outdir + "/" + wildcards.reference_version + "/STAR/full/" + i + ".aligned.bam")
    return(cond1, cond2)


with open("/Users/u1001407/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/configs/config_RNA-Seq.json") as data_file:
    config = json.load(data_file)

def dummFunc():
    cond1 = []
    cond2 = "condition2"
    for i in ("c1", "c2", "c3"):
        cond1.append(i)
    return(cond1, cond2)
