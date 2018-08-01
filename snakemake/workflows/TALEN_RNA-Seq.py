_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.3

localrules:
    all

home = os.environ['HOME']

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/rules/"

include:
    include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_kallisto.py"
include:
    include_prefix + "run_STAR.py"
include:
    include_prefix + "run_bowtie2.py"
include:
    include_prefix + "run_express.py"

rule run_kallisto:
    input:
        expand("{outdir}/{reference_version}/kallisto/{unit}",
               outdir = config["processed_dir"],
               reference_version = "GRCm38_ensembl84_cDNA",
               unit = config["units"])

rule run_bowtie2:
    input:
        expand("{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{condition}{unit}/hits.bam",
                outdir = config["processed_dir"],
                reference_version="GRCm38_ensembl84",
                condition = "wt",
                unit = ["1", "2", "3"]),
        expand("{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{condition}{unit}/hits.bam",
                outdir = config["processed_dir"],
                reference_version="GRCm38_ensembl84",
                condition = "hemi",
                unit = ["1", "2", "3"]),

rule run_express:
    input:
        expand("{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{condition}{unit}/express/results.xprs",
                outdir = config["processed_dir"],
                reference_version="GRCm38_ensembl84",
                condition = "wt",
                unit = ["1", "2", "3"]),
        expand("{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{condition}{unit}/express/results.xprs",
                outdir = config["processed_dir"],
                reference_version="GRCm38_ensembl84",
                condition = "hemi",
                unit = ["1", "2", "3"]),

rule convert_bam_to_bw:
    input:
        expand("{outdir}/{reference_version}/deepTools/bamCoverage/{unit}.bw",
               outdir=config["processed_dir"],
               reference_version="GRCm38_ensembl84",
               unit=config["units"])

rule convert_rMATS_bam_to_bw:
    input:
        expand("{outdir}/{reference_version}/rMATS/BWs/{unit}.bw",
                outdir = config["processed_dir"],
                reference_version = "GRCm38_ensembl84",
                unit = config["units"])

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
