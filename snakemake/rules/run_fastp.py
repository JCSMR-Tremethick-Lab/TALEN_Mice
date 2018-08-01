__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-05-18"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with fastq
(https://github.com/OpenGene/fastp)

For usage, include this in your workflow.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException

home = os.environ['HOME']

configfile: home + "/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/configs/config_CUT_and_RUN.json"

rule run_fastp:
    version:
        0.1
    threads:
        4
    input:
        read1 = lambda wildcards: wildcards["assayType"] + "/fastq/" + wildcards["runID"] + "/" + config["samples"][wildcards["assayType"]][wildcards["runID"]][wildcards["library"]][0],
        read2 = lambda wildcards: wildcards["assayType"] + "/fastq/" + wildcards["runID"] + "/" + config["samples"][wildcards["assayType"]][wildcards["runID"]][wildcards["library"]][1]
    output:
        trimmed_read1 = "{assayType}/trimmed/{runID}/{library}.end1.fastq.gz",
        trimmed_read2 = "{assayType}/trimmed/{runID}/{library}.end2.fastq.gz",
        report_html = "{assayType}/trimmed/{runID}/{library}.report.html",
        report_json = "{assayType}/trimmed/{runID}/{library}.report.json"
    shell:
        "fastp -i {input.read1} -I {input.read2} -o {output.trimmed_read1} -O {output.trimmed_read2} --html {output.report_html} --json {output.report_json} --thread {threads}"
