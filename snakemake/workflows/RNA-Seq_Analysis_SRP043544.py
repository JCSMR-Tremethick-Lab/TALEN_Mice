_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException
from os.path import join
import os

rule:
    version: 0.4

localrules:
    all


home = os.environ['HOME']


include_prefix= os.environ['HOME'] + "/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/rules/"

rule star_align_full:
    version:
        0.5
    threads:
        lambda wildcards: int(str(config["STAR"]["runThreadN"]).strip("['']"))
    params:
        tmp_dir = home + "/tmp"
    input:
        expand("{raw_dir}/{unit}",
                raw_dir = config["raw_dir"],
                unit = [ y \
                            for x in config["units"].keys() \
                                for y in config["units"][x]]),
        index = lambda wildcards: config["STAR"][wildcards.reference_version]
    output:
        "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {threads} \
                 --genomeDir {input.index} \
                 --readFilesIn {input[0]} {input[1]} \
                 --readFilesCommand zcat \
                 --outTmpDir {params.tmp_dir}/{wildcards.unit} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 --quantMode GeneCounts \
                 > {output}
        """

rule bam_index_STAR_output:
    version:
        0.2
    input:
        "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    output:
        "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai"
    wrapper:
        "file://" + wrapper_dir + "/samtools/index/wrapper.py"

rule kallisto_quant:
    message:
        "Running kallisto with paired end data..."
    params:
        raw_data = config["raw_dir"],
        outdir = config["processed_dir"],
        bootstraps = config["kallisto"]["bootstraps"],
    input:
        expand("{raw_dir}/{unit}",
                raw_dir = config["raw_dir"],
                unit = [ y \
                            for x in config["units"].keys() \
                                for y in config["units"][x]]),
        ki = lambda wildcards: home + config["kallisto_index"][wildcards.reference_version]
    output:
        "{outdir}/{reference_version}/kallisto/{unit}"
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                           {input[0]} {input[1]}
        """

rule all:
    input:
        expand("{outdir}/{reference_version}/kallisto/{unit}",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["units"])
