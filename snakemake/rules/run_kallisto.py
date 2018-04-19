__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

rule kallisto_quant_se:
    threads:
        4
    params:
        raw_data = config["raw_dir"],
        outdir = config["processed_dir"],
        bootstraps = config["kallisto"]["bootstraps"],
    input:
        read1 = lambda wildcards: "fastq/" + config["units"][wildcards.unit],
        ki = lambda wildcards: config["kallisto_index"][wildcards.reference_version]
    output:
        "{outdir}/{reference_version}/kallisto_se/{unit}"
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads={threads} \
                           --single \
                           --fragment-length=200 \
                           --sd=50 \
                           --bootstrap-samples={params.bootstraps} \
                           {input.read1}
        """

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

rule kallisto_quant:
    threads:
        4
    params:
        raw_data = config["raw_dir"],
        outdir = config["processed_dir"],
        bootstraps = config["kallisto"]["bootstraps"],
    input:
        rules.cutadapt_pe.output,
        ki = lambda wildcards: config["kallisto_index"][wildcards.reference_version]
    output:
        "{outdir}/{reference_version}/kallisto/{unit}"
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads={threads} \
                           --bootstrap-samples={params.bootstraps} \
                           {input[0]} {input[1]}
        """
