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

rule run_STAR_full:
    input:
        expand("{assayType}/igvtools/count/{reference_version}/{runID}/{library}/{library}.tdf",
               assayType = "RNA-Seq",
               reference_version = ["GRCm38_ensembl93_ERCC"],
               runID = "NB501086_0219_TSoboleva_JCSMR_RNAseq",
               library = [y for y in config["samples"]["RNA-Seq"]["NB501086_0219_TSoboleva_JCSMR_RNAseq"].keys()])

rule run_kallisto:
    input:
        expand("{assayType}/kallisto/quant/{reference_version}/{runID}/{library}",
               assayType = "RNA-Seq",
               reference_version = ["GRCm38_ensembl84_ERCC", "GRCm38_gencodeM18", "GRCm38_ensembl93_ERCC"],
               runID = "NB501086_0219_TSoboleva_JCSMR_RNAseq",
               library = [y for y in config["samples"]["RNA-Seq"]["NB501086_0219_TSoboleva_JCSMR_RNAseq"].keys()])

rule run_kallisto_pseudoalignment:
    input:
        expand("{assayType}/kallisto/genomebam/{reference_version}/{runID}/{library}",
               assayType = "RNA-Seq",
               reference_version = "GRCm38_ensembl93_ERCC",
               runID = "NB501086_0219_TSoboleva_JCSMR_RNAseq",
               library = [y for y in config["samples"]["RNA-Seq"]["NB501086_0219_TSoboleva_JCSMR_RNAseq"].keys()])

rule run_pseudoalignment_to_bigwig:
    input:
        expand("{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{library}_RPKM.bw",
               assayType = "RNA-Seq",
               reference_version = "GRCm38_ensembl93_ERCC",
               runID = "NB501086_0219_TSoboleva_JCSMR_RNAseq",
               library = [y for y in config["samples"]["RNA-Seq"]["NB501086_0219_TSoboleva_JCSMR_RNAseq"].keys()])
