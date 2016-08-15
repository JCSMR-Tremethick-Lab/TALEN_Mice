import json
from pprint import pprint

def getFASTQ(wildcards, assayID):
    fn = []
    for i in wildcards[assayID]:
        fn.append("./" + i)
    return(fn)

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

wildcards = dict()
wildcards = {"processed_dir" : "processed_data",
             "genome_version" : "hg38",
             "units" : "10_9_hemi_OB_CCGTCC",
             "tissue" : "OB",
             "condition" : "WT_vs_HEMI"}

with open("config.json") as data_file:
    config = json.load(data_file)

def dummFunc():
    cond1 = []
    cond2 = "condition2"
    for i in ("c1", "c2", "c3"):
        cond1.append(i)
    return(cond1, cond2)
