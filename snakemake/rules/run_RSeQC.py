__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-08-19"

from snakemake.exceptions import MissingInputException

rule read_distribution:
    input:
        expand("{outdir}/{reference_version}/RSeQC/read_distribution/{unit}.txt",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["units"])

rule run_read_distributions:
    message:
        "Running RSeQC read_distribution.py"
    params:
        unit = lambda wildcards: wildcards["unit"],
        binary = config["RSeQC"]["binaries"]["read_distribution.py"]
    input:
        bam = "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
        gene_models = config["RSeQC"]["bedFiles"]["genes"]
    output:
        "{outdir}/{reference_version}/RSeQC/read_distribution/{unit}.txt"
    shell:
        """
            {params.binary} -i {input.bam} -r {input.gene_models} > {output}
        """
