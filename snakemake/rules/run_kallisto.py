__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

rule kallisto_quant:
    threads:
        4
    params:
        bootstraps = config["kallisto"]["bootstraps"],
	    kallisto_bin = home + "/miniconda3/envs/kallisto/bin/kallisto"
    input:
        trimmed_read1 = rules.run_fastp.output.trimmed_read1,
        trimmed_read2 = rules.run_fastp.output.trimmed_read2,
        ki = lambda wildcards: config["kallisto_index"][wildcards.reference_version]
    output:
        directory("{assayType}/kallisto/quant/{reference_version}/{runID}/{library}")
    shell:
        """
            {params.kallisto_bin} quant --index={input.ki} \
                           --output-dir={output} \
                           --threads={threads} \
                           --bootstrap-samples={params.bootstraps} \
                           {input[0]} {input[1]}
        """

rule kallisto_quant_pseudo:
    threads:
        4
    params:
        bootstraps = config["kallisto"]["bootstraps"],
	    kallisto_bin = home + "/miniconda3/envs/kallisto/bin/kallisto"
    input:
        trimmed_read1 = rules.run_fastp.output.trimmed_read1,
        trimmed_read2 = rules.run_fastp.output.trimmed_read2,
        ki = lambda wildcards: config["kallisto_index"][wildcards.reference_version],
        gtf = lambda wildcards: config["STAR"][wildcards.reference_version]["GTF"]
    output:
        directory("{assayType}/kallisto/genomebam/{reference_version}/{runID}/{library}")
    shell:
        """
            {params.kallisto_bin} quant --index={input.ki} \
                                        --output-dir={output} \
                                        --threads={threads} \
                                        --bootstrap-samples={params.bootstraps} \
                                        --genomebam \
                                        --gtf {input.gtf} \
                                        {input[0]} {input[1]}
        """


rule pseudobam_to_bigwig:
    version:
        1
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        outFileFormat = "bigwig",
        binSize = 10,
        smoothLength = 30,
        normalizeUsing = "RPKM",
        maxFragmentLength = "50000"
    threads:
        32
    input:
        rules.kallisto_quant_pseudo.output
    output:
        "{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{library}_RPKM.bw"
    shell:
        """
        {params.deepTools_dir}/bamCoverage --bam {input}/pseudoalignments.bam \
                                           --outFileName {output} \
                                           --outFileFormat {params.outFileFormat} \
                                           --binSize {params.binSize} \
                                           --smoothLength {params.smoothLength}\
                                           --numberOfProcessors {threads} \
                                           --normalizeUsing RPKM \
                                           --maxFragmentLength {params.maxFragmentLength} \
                                           --ignoreForNormalization {params.ignore}
        """
