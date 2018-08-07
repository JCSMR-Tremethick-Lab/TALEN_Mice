__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-07"

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
SPIKEIN_GENOME = config["references"]["spikeIn"]
SPIKEIN_VERSION = config["references"][SPIKEIN_GENOME]["version"][0]


rule bowtie2_pe_spikeIn:
    version:
        "1"
    params:
        max_in = config["program_parameters"]["bt2_params"]["max_insert"],
        bt2_index = config["references"][SPIKEIN_GENOME]["bowtie2"][SPIKEIN_VERSION]
    threads:
        lambda wildcards: int(str(config["program_parameters"]["bt2_params"]["threads"]).strip("['']"))
    input:
        read1 = "{assayType}/unmapped_reads/" + REF_VERSION + "/{runID}/{library}.unmapped_r1.fastq.gz",
        read2 = "{assayType}/unmapped_reads/" + REF_VERSION + "/{runID}/{library}.unmapped_r2.fastq.gz"
    output:
        protected("{assayType}/bowtie2/{spikein_version}/{runID}/{library}.bam")
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
        protected("{assayType}/bowtie2/{spikein_version}/{runID}/{library}.bam.stats.txt")
    shell:
        """
            samtools flagstat {input} > {output}
        """


# rules
rule bam_quality_filter:
    version:
        "1.0"
    params:
        qual = config["alignment_quality"]
    input:
        rules.bowtie2_pe.output
    output:
        temp("{assayType}/samtools/quality_filtered/{spikein_version}/{runID}/{library}.bam")
    shell:
        "samtools view -b -h -q {params.qual} {input} > {output}"


rule bam_sort:
    version:
        "1.0"
    threads:
        4
    input:
        rules.bam_quality_filter.output
    output:
        temp("{assayType}/samtools/sorted/{spikein_version}/{runID}/{library}.bam")
    shell:
        "samtools sort -@ {threads} {input} -T {wildcards.library}.sorted -o {output}"


rule bam_mark_duplicates:
    version:
        "1.0"
    params:
        qual = config["alignment_quality"],
        picard = home + config["program_parameters"]["picard_tools"]["jar"],
        temp = home + config["temp_dir"]
    input:
        rules.bam_sort.output
    output:
        out = temp("{assayType}/picardTools/MarkDuplicates/{spikein_version}/{runID}/{library}.bam"),
        metrics = protected("{assayType}/picardTools/MarkDuplicates/{spikein_version}/{runID}/{library}.metrics.txt")
    shell:
        """
            java -Djava.io.tmpdir={params.temp} \
            -Xmx24G \
            -jar {params.picard} MarkDuplicates \
            INPUT={input}\
            OUTPUT={output.out}\
            ASSUME_SORTED=TRUE\
            METRICS_FILE={output.metrics}
        """


rule bam_rmdup:
    input:
        rules.bam_mark_duplicates.output.out
    output:
        temp("{assayType}/samtools/rmdup/{spikein_version}/{runID}/{library}.bam")
    shell:
        "samtools rmdup {input} {output}"


rule bam_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_rmdup.output
    output:
        protected("{assayType}/samtools/rmdup/{spikein_version}/{runID}/{library}.bam.bai")
    shell:
        "samtools index {input} {output}"
