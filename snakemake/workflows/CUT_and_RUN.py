__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-01"

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
    include_prefix + "run_bowtie2.py"

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

rule convert_bam_to_bw:
    input:
        expand("{outdir}/{reference_version}/deepTools/bamCoverage/{unit}.bw",
               outdir=config["processed_dir"],
               reference_version="GRCm38_ensembl84",
               unit=config["units"])
rule all:
    input:
        expand("{assayType}/trimmed/{runID}/{library}.{suffix}",
            assayType = "CutRun",
            runID = "180731_NB501086_0217_CutandRun_Tanya",
            library = [""],
            suffix = ["end1.fastq.gz", "end2.fastq.gz", "report.html", "report.json"])
