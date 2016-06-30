__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-25"

from snakemake.exceptions import MissingInputException
import os

def getAllFASTQ():
    fn = []
    for i in config["units"]:
        for j in config["units"][i]:
            fn.append("./fastq/" + j)
    return(fn)

rule dummy:
    input:
        expand("{rdir}", rdir = config["reports_dir"])

rule fastqc:
    message: "Performing FastQC..."
    params:
        rdir = config["reports_dir"],
        data_dir = config["raw_dir"]
    input:
        getAllFASTQ
    output:
        "{rdir}/"
    shell:
        "/usr/local/bin/fastqc {input} --noextract --outdir  {params.rdir}"
