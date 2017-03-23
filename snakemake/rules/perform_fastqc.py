__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-25"

from snakemake.exceptions import MissingInputException
import os

def getAllFASTQ(wildcards):
    fn = []
    for i in config["units"]:
        for j in config["units"][i]:
            fn.append("fastq/" + i + "/" + j)
    return(fn)

rule dummy:
    input:
        expand("{rdir}", rdir = config["reports_dir"])

rule fastqc:
    params:
        rdir = config["reports_dir"],
        data_dir = config["raw_dir"]
    input:
        getAllFASTQ
    output:
        "reports"
    shell:
        "fastqc {input} --noextract --outdir  {output}"
