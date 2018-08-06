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

rule bam_extract_unmapped_reads:
    version:
        "1"
    input:
        rules.bowtie2_pe.output
    output:
        temp("{assayType}/bowtie2/{reference_version}/{runID}/{library}.unmapped.bam")
    shell:
        "samtools view -f 4 -b {input} > {output}"


rule bam_sort_unmapped_reads:
    version:
        "1"
    input:
        rules.bam_extract_unmapped_reads.output
    output:
        temp("{assayType}/bowtie2/{reference_version}/{runID}/{library}.unmapped.sorted.bam")
    shell:
        "samtools sort -n {input} -T {wildcards.library}.sorted -o {output}"


rule unmapped_reads_to_pe_fastq:
    version:
        "1"
    input:
        rules.bam_sort_unmapped_reads.output
    output:
        temp("{assayType}/unmapped_reads/{runID}/{library}.unmapped_r1.fastq"),
        temp("{assayType}/unmapped_reads/{runID}/{library}.unmapped_r2.fastq")
    shell:
        """
            bedtools bamtofastq -i {input} \
                                -fq {output[0]} \
                                -fq2 {output[1]}
        """


rule gzip_unmapped_fastq:
    version:
        "1"
    input:
        rules.unmapped_reads_to_pe_fastq.output
    output:
        protected("{assayType}/unmapped_reads/{runID}/{library}.unmapped_r1.fastq.gz"),
        protected("{assayType}/unmapped_reads/{runID}/{library}.unmapped_r1.fastq.gz")
    shell:
        "gzip {input[0]}; gzip {input[1]}"
