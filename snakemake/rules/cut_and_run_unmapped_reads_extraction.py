__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-06"

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
        temp("{assayType}/unmapped_reads/{reference_version}/{runID}/{library}.unmapped_r1.fastq"),
        temp("{assayType}/unmapped_reads/{reference_version}/{runID}/{library}.unmapped_r2.fastq")
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
        "{assayType}/unmapped_reads/{reference_version}/{runID}/{library}.unmapped_r1.fastq.gz",
        "{assayType}/unmapped_reads/{reference_version}/{runID}/{library}.unmapped_r2.fastq.gz"
    shell:
        "gzip {input[0]}; gzip {input[1]}"
