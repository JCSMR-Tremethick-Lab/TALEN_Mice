__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException
import os
import subprocess

def getFragmentSize(wildcards):
    insertSize = os.popen("/home/sebastian/miniconda3/envs/pizzly/bin/pizzly_get_fragment_length.py "\
                 + wildcards["assayType"] + "/kallisto/fusion/"\
                 + wildcards["reference_version"] + "/"\
                 + wildcards["runID"] + "/"\
                 + wildcards["library"] + "/"\
                 "abundance.h5").read()
    return(insertSize)

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

rule kallisto_fusion:
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
        directory("{assayType}/kallisto/fusion/{reference_version}/{runID}/{library}")
    shell:
        """
            {params.kallisto_bin} quant --index={input.ki} \
                                        --output-dir={output} \
                                        --threads={threads} \
                                        --bootstrap-samples={params.bootstraps} \
                                        --fusion \
                                        --gtf {input.gtf} \
                                        {input[0]} {input[1]}
        """


rule run_pizzly:
    threads:
        1
    params:
        pizzly_bin = home + "/miniconda3/envs/pizzly/bin/pizzly",
        k = 31,
        alignScore = 2,
        insertSize = getFragmentSize
    input:
        abundance = "{assayType}/kallisto/fusion/{reference_version}/{runID}/{library}",
        gtf = lambda wildcards: config["STAR"][wildcards.reference_version]["GTF"],
        fasta = lambda wildcards: config["STAR"][wildcards.reference_version]["fasta"]
    output:
        directory("{assayType}/pizzly/{reference_version}/{runID}/{library}/")
    shell:
        """
            {params.pizzly_bin} -k {params.k}\
                                --gtf {input.gtf}\
                                --cache {input.gtf}.cache.txt\
                                --fasta {input.fasta}\
                                --output {output}/{wildcards.library}\
                                {input.abundance}/fusion.txt
        """
