__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-04-19"

from snakemake.exceptions import MissingInputException

def getTPMs(wildcards):
    fn = []
    for i in config["Groups"][wildcards["tissue"]][wildcards["condition"]]:
        fn.append("/".join([wildcards["outdir"], wildcards["reference_version"], "suppa", i, "abundance.tpm"]))
    return(" ".join(fn))

rule make_tpm_tsv:
    threads:
        1
    params:
        sample = lambda wildcards: wildcards["library"]
    input:
        "{assayType}/kallisto/{reference_version}/{runID}/{library}/abundance.tsv"
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
        suppa_dir = home + "/miniconda3/bin/"
    input:
        getTPMs
    output:
        "{outdir}/{reference_version}/suppa/pooled/{tissue}/{condition}/abundance.tpm",
        directory("{outdir}/{reference_version}/suppa/pooled/{tissue}/{condition}/abundance")
    shell:
        """
            suppa.py joinFiles --file-extension tpm\
                               --input-files {input}\
                               --output {output[1]}
        """
