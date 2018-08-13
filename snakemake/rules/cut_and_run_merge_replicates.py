__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-13"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for merging BAM files.

For use, include in your workflow.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException
from snakemake.io import expand


# set local variables
home = os.environ['HOME']
REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"][0]

def getReplicates(wildcards):
    fn = []
    for i in config["samples"][wildcards["assayType"]]["replicates"][wildcards["runID"]][wildcards["replicate"]]:
        fn.append("/".join([wildcards["assayType"], "deepTools", "bamCoverage", wildcards["reference_version"], wildcards["runID"], i + wildcards["suffix"] + ".bw"]))
    return(" ".join(fn))


rule all:
    input:
        expand("{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{replicate}_{suffix}.bw",
              assayType = "CutRun",
              reference_version = "GRCm38_ensembl93_ERCC",
              runID = "180731_NB501086_0217_CutandRun_Tanya",
              replicate =  [x for x in config["samples"]["CutRun"]["replicates"]["180731_NB501086_0217_CutandRun_Tanya"].keys()],
              suffix = ["RPKM", "1xgenome"]),
        expand("{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{replicate}_{suffix}.bw",
              assayType = "CutRun",
              reference_version = "GRCm38_ensembl93_ERCC",
              runID = "NB501086_0221_TSoboleva_JCSMR_CutandRun",
              replicate =  [x for x in config["samples"]["CutRun"]["replicates"]["NB501086_0221_TSoboleva_JCSMR_CutandRun"].keys()],
              suffix = ["RPKM", "1xgenome"])


rule bigWig_merge:
    version:
        "1"
    params:
        bin = home + "/miniconda3/envs/deepTools/bin/bigWigMerge"
    threads:
        1
    input:
        bigWigFiles = getReplicates
    output:
        mergedBedGraph = "{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{replicate}_{suffix}.bdg"
    shell:
        """
           {params.bin} {input.bigWigFiles} {output.mergedBedGraph}
        """

rule bedGraph_to_bigWig:
    version:
        "1"
    params:
        bin = home + "/bin/bedGraphToBigWig",
        chromSizes = "/Data/References/Genomes/Mus_musculus/GRCm38_mm10/mm10chromInfo.txt"
    threads:
        1
    input:
        rules.bigWig_merge.output.mergedBedGraph
    output:
        bigWig = "{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{replicate}_{suffix}.bw"
    shell:
        """
           {params.bin} {input} {params.chromSizes} {output.bigWig}
        """
