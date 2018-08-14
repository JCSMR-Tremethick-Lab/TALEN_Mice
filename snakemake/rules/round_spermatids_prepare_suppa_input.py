__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-04-19"

from snakemake.exceptions import MissingInputException

def getTPMs(wildcards):
    fn = []
    for i in config["samples"][wildcards["assayType"]]["condition"][wildcards["runID"]][wildcards["condition"]]:
        fn.append("/".join([wildcards["assayType"], "suppa", wildcards["reference_version"], i, "abundance.tpm"]))
    return(" ".join(fn))


rule all:
    input:
        expand("{assayType}/suppa/pooled/{reference_version}/{runID}/{condition}",
                assayType = "RNA-Seq",
                reference_version = "GRCm38_ensembl93_ERCC",
                runID = "NB501086_0219_TSoboleva_JCSMR_RNAseq",
                condition = ["KO", "WT"])

rule make_tpm_tsv:
    threads:
        1
    params:
        sample = lambda wildcards: wildcards["unit"]
    input:
        "{assayType}/kallisto/genomebam/{reference_version}/{runID}/{library}/abundance.tsv"
    output:
        "{assayType}/suppa/{reference_version}/{runID}/{library}/abundance.tpm"
    shell:
        """
            awk '{{if(NR==1) {{print "{params.sample}";}} else {{split($1,a,"."); print a[1]"\t"$5}}}}' < {input} > {output}
        """

rule collate_samples:
    threads:
        1
    params:
        suppa_dir = home + config["suppa_dir"]
    input:
        getTPMs
    output:
        directory("{assayType}/suppa/pooled/{reference_version}/{runID}/{condition}")
    shell:
        """
            suppa.py joinFiles --file-extension tpm\
                               --input-files {input}\
                               --output {output}
        """
