__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-04-19"

from snakemake.exceptions import MissingInputException

def getTPMs(wildcards):
    fn = []
    for i in config["Groups"][wildcards["tissue"]][wildcards["condition"]]:
        fn.append("/".join([wildcards["outdir"], wildcards["reference_version"], "kallisto", i, "abundance.tpm.tsv"]))
    return(" ".join(fn))

rule make_tpm_tsv:
    threads:
        1
    input:
        "{outdir}/{reference_version}/kallisto/{unit}/abundance.tsv"
    output:
        "{outdir}/{reference_version}/kallisto/{unit}/abundance.tpm.tsv"
    shell:
        """
            awk '{split($1,a,"."); print a[1]"\t"$5}' < {input}/abundance.tsv > {output}
        """

rule collate_samples:
    threads:
        1
    params:
        suppa_dir = home + config["suppa_dir"]
    input:
        replicates = getTPMs
    output:
        "{outdir}/{reference_version}/suppa/{tissue}/{condition}/abundances.tpm",
        "{outdir}/{reference_version}/suppa/{tissue}/{condition}/abundances"
    shell:
        """
            suppa.py joinFiles --file-extension tpm\
                               --input-files {input.replicates}\
                               --output {output[1]}
        """
