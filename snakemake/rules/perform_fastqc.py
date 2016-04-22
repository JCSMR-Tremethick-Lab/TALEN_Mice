__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-25"

from snakemake.exceptions import MissingInputException

rule fastqc:
    message: "Performing FastQC..."
    params:
        rdir = config["reports_dir"],
        data_dir = config["raw_dir"]
    input:
        expand("{data}/{sample}{suffix}.fastq.gz", sample = config["units"], data = config["raw_dir"], suffix = ("_R1_001", "_R2_001")),
    output:
        expand("{rdir}/{sample}{suffix}_fastqc.zip", sample = config["units"], rdir = config["reports_dir"], suffix = ("_R1_001", "_R2_001"))
    shell:
        "/usr/local/bin/fastqc {input} --noextract --outdir  {params.rdir}"

# rule all:
#     input:
#         expand("{rdir}/{sample}{suffix}_fastqc.zip", sample = config["units"], rdir = config["reports_dir"], suffix = ("_R1_001", "_R2_001"))
