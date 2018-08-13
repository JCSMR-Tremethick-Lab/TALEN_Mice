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

# set local variables
home = os.environ['HOME']
REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"][0]


rule bigWig_merge:
    version:
        "1"
    params:
        bin = home + "/miniconda3/envs/deepTools/bin/bigWigMerge"
    threads:
    input:
        bigWigFiles = expand("{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{library}_{suffix}.gz",
                             assayType = wildcards["assayType"],
                             reference_version = wildcards["reference_version"],
                             runID = wildcards["runID"],
                             library = [x for x in config["samples"]["replicates"][wildcards["runID"]][wildcards["replicate"]]],
                             suffix = wildcards["suffix"])
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
        bin = home + "/bin/bedGraphToBigWig"
        chromSizes = "/Data/References/Genomes/Mus_musculus/GRCm38_mm10/mm10chromInfo.txt"
    threads:
    input:
        rules.output.bigWig_merge.mergedBedGraph
    output:
        bigWig = "{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{replicate}_{suffix}.bw"
    shell:
    """
    {params.bin} {input} {params.chromSizes} {output.bigWig}
    """
