_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException

rule:
    version: 0.3

localrules:
    all

include_prefix="/home/skurscheid/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/rules/"

# include:
#     include_prefix + "perform_fastqc.py"
include:
    include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_kallisto.py"
include:
    include_prefix + "run_STAR.py"

rule all:
    input:
        expand("{outdir}/{reference_version}/kallisto_se/{unit}",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["units"])
