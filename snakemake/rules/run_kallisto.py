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

rule kallisto_quant:
    threads:
        4
    params:
        bootstraps = config["kallisto"]["bootstraps"],
    input:
        trimmed_read1 = rules.run_fastp.output.trimmed_read1,
        trimmed_read2 = rules.run_fastp.output.trimmed_read2,
        ki = lambda wildcards: config["kallisto_index"][wildcards.reference_version]
    output:
        "{assayType}/kallisto/{reference_version}/{runID}/{unit}"
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads={threads} \
                           --bootstrap-samples={params.bootstraps} \
                           {input[0]} {input[1]}
        """
