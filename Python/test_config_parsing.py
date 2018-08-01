import json
from pprint import pprint

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
