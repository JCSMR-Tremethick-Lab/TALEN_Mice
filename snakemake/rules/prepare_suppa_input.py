__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-04-19"

from snakemake.exceptions import MissingInputException

def getTPMs(wildcards):
    fn = []
    for i in config["Groups"][wildcards["tissue"]][wildcards["condition"]]:
        fn.append("/".join([wildcards["outdir"], wildcards["reference_version"], "suppa", i, "abundance.tpm.tsv"]))
    return(" ".join(fn))

rule collate_samples:
    threads:
        1
    params:
        suppa_dir = home + config["suppa_dir"]
    input:
        replicates = getTPMs
    output:
        "{outdir}/{reference_version}/suppa/pooled/{tissue}/{condition}/abundances.tpm",
        "{outdir}/{reference_version}/suppa/pooled/{tissue}/{condition}/abundances"
    shell:
        """
            suppa.py joinFiles --file-extension tpm\
                               --input-files {input.replicates}\
                               --output {output[1]}
        """
