__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-01"

from snakemake.exceptions import MissingInputException
import os

REF_GENOME = config["references"]["active"]

rule:
    version:
        "1.0"

localrules:
    all

home = os.environ['HOME']

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/rules/"

include:
    include_prefix + "run_fastp.py"
include:
    include_prefix + "run_bowtie2_cut_and_run.py"


rule all:
    input:
        expand("{assayType}/bowtie2/{reference_version}/{runID}/{library}.{suffix}",
                assayType = "CutRun",
                reference_version = config["references"][REF_GENOME]["version"],
                runID = "180731_NB501086_0217_CutandRun_Tanya",
                library = ["WT_01_H2AL2_3_7_18", "WT_01_IGG_3_7_18", "KO_01_H2AL2_3_7_18", "KO_02_H2AL2_24_6_18", "WT_02_H2AL2_24_6_18", "WT_01_H3K27me3_23_5_18", "WT_01_H3K36me3_23_5_18"],
                suffix = ["bam"])
