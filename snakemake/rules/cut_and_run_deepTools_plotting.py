__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-10"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools analysis on ChIP-Seq data
For usage, include this in your workflow.
"""

rules all:
    input:
        expand("{assayType}/deepTools/computeCoverage/{reference_version}/{runID}/{library}_RPKM.bw",
               assayType = "CutRun",
               reference_version = "GRCm38_ensembl93",
               runID = "NB501086_0221_TSoboleva_JCSMR_CutandRun",
               library [x for x in config["samples"][wildcards["assayType"]][wildcards["runID"]].keys()])

def get_computeMatrix_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayType"],
                     "samtools",
                     "rmdup",
                     wildcards["reference_version"],
                     wildcards["runID"]))
    for i in config["samples"][wildcards["assayType"]][wildcards["runID"]]:
        fn.append("/".join((path, "_".join((i, wildcards["mode"], "RPKM.bw")))))
    return(fn)


def cli_parameters_computeMatrix(wildcards):
    a = config["program_parameters"]["deepTools"][wildcards["tool"]][wildcards["command"]]
    if wildcards["command"] == "reference-point":
        a["--referencePoint"] = wildcards.referencePoint
    return(a)

def cli_parameters_normalization(wildcards):
    if wildcards["norm"] == "RPKM":
        a = "--normalizeUsingRPKM"
    elif wildcards["norm"] == "1xcoverage":
        a = " ".join(("--normalizeTo1x", config["references"][REF_GENOME]["effectiveSize"]))
    return(a)

def cli_parameters_bamCoverage(wildcards):
    a = config["program_parameters"]["deepTools"]["bamCoverage"]["normal"]
    b = str()
    for (key, val) in a.items():
        if val == " ":
            f = key + " "
            b = b + f
        else:
            f = key + "=" + val + " "
            b = b + f
    return(b.rstrip())

def get_computeMatrix_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayID"],
                     wildcards["runID"],
                     config["processed_dir"],
                     config["references"][REF_GENOME]["version"][0],
                     wildcards["application"],
                     "bamCoverage",
                     wildcards["mode"],
                     wildcards["duplicates"]))
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        fn.append("/".join((path, "_".join((i, wildcards["mode"], "RPKM.bw")))))
    return(fn)


rule bamCoverage_normal:
    version:
        1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"]
    threads:
        32
    input:
        bam = "{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.bam",
        index = "{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.bam.bai"
    output:
        "{assayType}/deepTools/computeCoverage/{reference_version}/{runID}/{library}_RPKM.bw"
    shell:
        """
        {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                           --outFileName {output} \
                                           --outFileFormat bigwig \
                                           --centerReads,
                                           --binSize 25,
                                           --smoothLength "75",
                                           --numberOfProcessors {threads} \
                                           --normalizeUsingRPKM \
                                           --ignoreForNormalization {params.ignore}
        """

# rule computeMatrix:
#     version:
#         "1"
#     params:
#         deepTools_dir = home + config["deepTools_dir"],
#         program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
#     threads:
#         32
#     input:
#         file = get_computeMatrix_input,
#         region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards.region]
#     output:
#         matrix_gz = "{assayType}/deepTools/computeMatrix/{reference_version}/{runID}/{region}/matrix.gz"
#     wrapper:
#         "file://" + wrapper_dir + "/deepTools/computeMatrix/wrapper.py"
