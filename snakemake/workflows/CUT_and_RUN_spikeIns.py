__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-01"

from snakemake.exceptions import MissingInputException
import os

REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"]

SPIKEIN_GENOME = config["references"]["spikeIn"]
SPIKEIN_VERSION = config["references"][SPIKEIN_GENOME]["version"][0]

rule:
    version:
        "1.0"

localrules:
    all

home = os.environ['HOME']

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/rules/"

include:
    include_prefix + "cut_and_run_second_alignment.py"

rule all:
    input:
        expand("{assayType}/samtools/rmdup/{reference_version}/{runID}/spikeIn/{library}.{suffix}",
                assayType = "CutRun",
                reference_version = SPIKEIN_VERSION,
                runID = "180731_NB501086_0217_CutandRun_Tanya",
                library = ["WT_01_H2AL2_3_7_18", "WT_01_IGG_3_7_18", "KO_01_H2AL2_3_7_18", "KO_02_H2AL2_24_6_18", "WT_02_H2AL2_24_6_18", "WT_01_H3K27me3_23_5_18", "WT_01_H3K36me3_23_5_18"],
                suffix = ["bam.bai"]),
        expand("{assayType}/samtools/flagstat/{reference_version}/{runID}/spikeIn/{library}.{suffix}",
                assayType = "CutRun",
                reference_version = SPIKEIN_VERSION,
                runID = "180731_NB501086_0217_CutandRun_Tanya",
                library = ["WT_01_H2AL2_3_7_18", "WT_01_IGG_3_7_18", "KO_01_H2AL2_3_7_18", "KO_02_H2AL2_24_6_18", "WT_02_H2AL2_24_6_18", "WT_01_H3K27me3_23_5_18", "WT_01_H3K36me3_23_5_18"],
                suffix = ["bam.stats.txt"])

