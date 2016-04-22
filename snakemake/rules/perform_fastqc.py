__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-25"

from snakemake.exceptions import MissingInputException
import os

rule fastqc:
    message: "Performing FastQC..."
    params:
        rdir = config["reports_dir"],
        data_dir = config["raw_dir"]
    input:
        "fastqc/{unit}_R2_001.fastq.gz",
        "fastqc/{unit}_R1_001.fastq.gz"
    output:
        "{rdir}/{unit}"
    shell:
        "/usr/local/bin/fastqc {params.data_dir}/{input}  --noextract --outdir {output}"

rule all:
    input:
        expand("{rdir}/{unit}", unit = config["units"], rdir = config["reports_dir"])
