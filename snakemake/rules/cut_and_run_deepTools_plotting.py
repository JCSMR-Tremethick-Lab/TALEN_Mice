__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-10"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os
home = os.environ['HOME']

"""
Rules for running deepTools analysis on ChIP-Seq data
For usage, include this in your workflow.
"""

rules all:
    input:
        expand("{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{library}_RPKM.bw",
               assayType = "CutRun",
               reference_version = "GRCm38_ensembl93",
               runID = "NB501086_0221_TSoboleva_JCSMR_CutandRun",
               library [x for x in config["samples"][wildcards["assayType"]][wildcards["runID"]].keys()]),
        expand("{assayType}/deepTools/computeMatrix/scale-region/{reference_version}/{runID}/{region}/matrix.gz",
               assayType = "CutRun",
               reference_version = "GRCm38_ensembl93",
               runID = ["NB501086_0221_TSoboleva_JCSMR_CutandRun", "180731_NB501086_0217_CutandRun_Tanya"],
               region = "allGenes")


def get_computeMatrix_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayType"],
                     "deepTools",
                     "bamCoverage",
                     wildcards["reference_version"],
                     wildcards["runID"]))
    for i in config["samples"][wildcards["assayType"]][wildcards["runID"]]:
        fn.append("/".join((path, "_".join((i, "RPKM.bw")))))
    return(fn)

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
        "{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{library}_RPKM.bw"
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

rule computeMatrix_scaled:
    version:
        "1"
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = ,
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"]
    threads:
        32
    input:
        file = get_computeMatrix_input,
        region = lambda wildcards: config["program_parameters"]["deepTools"]["regionFiles"][wildcards["region"]]
    output:
        matrix_gz = "{assayType}/deepTools/computeMatrix/scale-region/{reference_version}/{runID}/{region}/matrix.gz"
    shell:
        """
        {params.deepTools_dir}/computeCoverage scale-regions --numberOfProcessors {threads} \
                                                             --smartLabels \
                                                             --missingDataAsZero \
                                                             --regionBodyLength 5000 \
                                                             --beforeRegionStartLength 2000 \
                                                             --afterRegionStartLength \
                                                             --unscaled5prime 350 \
                                                             --unscaled3prime 350 \
                                                             --regionsFileName {input.region} \
                                                             --scoreFileName {input.file} \
                                                             --outFileName {output.matrix_gz}
        """