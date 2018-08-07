_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-07"

from snakemake.exceptions import MissingInputException
import os

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
    include_prefix + "run_kallisto.py"
include:
    include_prefix + "run_STAR.py"

rule run_kallisto:
    input:
        expand("{assayType}/kallisto/{reference_version}/{runID}/{library}",
               assayType = "RNA-Seq",
               reference_version = "GRCm38_ensembl84_ERCC",
               runID = "NB501086_0219_TSoboleva_JCSMR_RNAseq",
               library = [y for y in config["samples"]["NB501086_0219_TSoboleva_JCSMR_RNAseq"].keys()])

rule all:
    input:
        expand("./{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               unit = config["units"],
               trim_data = config["trim_dir"],
               suffix = ["R1_001", "R2_001"]),
        expand("{outdir}/{reference_version}/kallisto/{unit}",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["units"]),
        expand("{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["units"]),
        expand("{outdir}/{reference_version}/PICARD/insert_size_metrics/{unit}.insert_size_metrics.{suffix}",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["units"],
               suffix = ("pdf", "txt")),
        expand("{outdir}/{reference_version}/DEXSeq/count/{unit}.txt",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["units"]),
