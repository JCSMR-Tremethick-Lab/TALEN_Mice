__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-01"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for performing macs2 peak calling.

For use, include in your workflow.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException

# set local variables
home = os.environ['HOME']
REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"][0]


rule run_macs2_narrow:
    threads:
        1
    params:
        macs2_bin = home + "/miniconda3/envs/py27/bin/macs2",
        seed = "1234"
    input:
        bam = "{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.bam",
        bai = "{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.bam.bai"
    output:
        directory("{assayType}/macs2/callpeak/narrow/{reference_version}/{runID}/{library}")
    shell:
        """
            {params.macs2_bin} callpeak -f AUTO\
                               --seed {params.seed}\
                               --cutoff-analysis\
                               -g mm\
                               --bdg \
                               --SPMR \
                               --treatment {input.bam} \
                               --name {wildcards.library} \
                               --call-summits \
                               --outdir {output}
        """

rule run_macs2_broad:
    threads:
        1
    params:
        macs2_bin = home + "/miniconda3/envs/py27/bin/macs2",
        seed = "1234"
    input:
        bam = "{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.bam",
        bai = "{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.bam.bai"
    output:
        directory("{assayType}/macs2/callpeak/broad/{reference_version}/{runID}/{library}")
    shell:
        """
            {params.macs2_bin} callpeak -f AUTO\
                               --seed {params.seed}\
                               --cutoff-analysis\
                               -g mm\
                               --bdg \
                               --SPMR \
                               --treatment {input.bam} \
                               --name {wildcards.library} \
                                --broad \
                               --outdir {output}
        """

