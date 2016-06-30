import json
from pprint import pprint

def getFASTQ(wildcards):
    fn = []
    for i in config["units"][wildcards["units"]]:
        fn.append("./" + wildcards["processed_dir"] + "/" + wildcards["genome_version"] + "/fastq/" + i)
    return(fn)

wildcards = dict()
wildcards = {"processed_dir" : "processed_data", "genome_version" : "hg38", "units" : "10_9_hemi_OB_CCGTCC"}

with open("config_RNA-Seq_Brain.json") as data_file:
    config = json.load(data_file)
