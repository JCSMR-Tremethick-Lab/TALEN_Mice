_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-04-19"

from snakemake.exceptions import MissingInputException

rule:
    version: 0.3

localrules:
    all

home = os.environ['HOME']

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/rules/"

include:
    include_prefix + "prepare_suppa_input.py"

rule make_tpm_tsv_files:
    input:
        expand("{outdir}/{reference_version}/suppa/{unit}/abundance.tpm",
               outdir = config["processed_dir"],
               reference_version = "GRCm38_ensembl84_cDNA",
               unit = config["units"])

rule prepare_suppa:
    input:
        expand("{outdir}/{reference_version}/suppa/pooled/{tissue}/{condition}/abundance.tpm",
                outdir = config["processed_dir"],
                reference_version = "GRCm38_ensembl84_cDNA",
                tissue = ["PFC", "OB", "HIPPO"],
                condition = ["WT", "HEMI"])
