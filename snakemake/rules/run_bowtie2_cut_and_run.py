__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-01"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for aligning paired-end reads using bowtie2.

For use, include in your workflow.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException

# set local variables
home = os.environ['HOME']
REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"][0]


rule bowtie2_pe:
    version:
        "1"
    params:
        max_in = config["program_parameters"]["bt2_params"]["max_insert"],
        bt2_index = config["references"][REF_GENOME]["bowtie2"][REF_VERSION]
    threads:
        lambda wildcards: int(str(config["program_parameters"]["bt2_params"]["threads"]).strip("['']"))
    input:
        trimmed_read1 = "{assayType}/trimmed/{runID}/{library}.end1.fastq.gz",
        trimmed_read2 = "{assayType}/trimmed/{runID}/{library}.end2.fastq.gz"
    output:
        protected("{assayType}/bowtie2/{reference_version}/{runID}/{library}.bam")
    shell:
        """
            bowtie2 \
            -x {params.bt2_index}\
            --no-mixed \
            --no-discordant \
            --maxins {params.max_in} \
            --threads {threads}\
            --rg-id '{wildcards.library}' \
            --rg 'LB:{wildcards.library}' \
            --rg 'SM:{wildcards.library}' \
            --rg 'PL:Illumina' \
            --rg 'PU:NA' \
            -1 {input.trimmed_read1} \
            -2 {input.trimmed_read2} \
            | samtools view -Sb - > {output}
        """


rule bam_stats:
    version:
        "1"
    input:
        rules.bowtie2_pe.output
    output:
        protected("{assayType}/bowtie2/{reference_version}/{runID}/{library}.bam.stats.txt")
    shell:
        """
            samtools flagstat {input} > {output}
        """
