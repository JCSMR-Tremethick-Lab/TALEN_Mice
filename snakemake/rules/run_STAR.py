__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-22"

from snakemake.exceptions import MissingInputException

rule all:
    input:
        expand("{outdir}/STAR/{unit}", outdir = config["processed_dir"], unit = config["units"])

rule star_align:
    params:
        genomeDir = config["STAR"]["genomeDir"],
        runThreadN = config["STAR"]["runThreadN"]
    input:
        "fastq/subsets/{unit}_subset_R1_001.fastq.gz",
        "fastq/subsets/{unit}_subset_R2_001.fastq.gz"
    output:
        "{outdir}/STAR/{unit}"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {params.genomeDir} \
                 --readFilesIn {input[0]} {input[1]} \
                 --outFileNamePrefix {output[0]}
        """
