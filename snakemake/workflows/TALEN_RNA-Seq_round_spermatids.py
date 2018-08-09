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
#include:
#    include_prefix + "run_STAR.py"

rule run_kallisto:
    input:
        expand("{assayType}/kallisto/{reference_version}/{runID}/{library}",
               assayType = "RNA-Seq",
               reference_version = ["GRCm38_ensembl84_ERCC", "GRCm38_gencodeM18", "GRCm38_ensembl93_ERCC"],
               runID = "NB501086_0219_TSoboleva_JCSMR_RNAseq",
               library = [y for y in config["samples"]["RNA-Seq"]["NB501086_0219_TSoboleva_JCSMR_RNAseq"].keys()])

