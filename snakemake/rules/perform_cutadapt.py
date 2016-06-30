i__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-27"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with cutadapt
(http://cutadapt.readthedocs.org/en/latest/guide.html#illumina-truseq)

For usage, include this in your workflow.
"""

UNIT_TO_SAMPLE = {
    unit: sample for sample, units in config["units"].items()
    for unit in units}

# configfile: "config.json"

rule cutadapt_pe:
    """Trims given paired-end reads with given parameters"""
    params:
        trim_params = config["trim_params"],
        trim_data = config["trim_dir"],
        raw_data = config["raw_dir"]
    input:
        "./NB501086_0034_MNekrasov_JCSMR_ChIPseq/{units}_R1_001.fastq.gz",
        "./NB501086_0034_MNekrasov_JCSMR_ChIPseq/{units}_R2_001.fastq.gz"
    output:
        "./trimmed_data/{units}_R1_001.QT.CA.fastq.gz",
        "./trimmed_data/{units}_R2_001.QT.CA.fastq.gz"
    shell:
        """
            cutadapt {params.trim_params} \
            -o {output[0]} -p {output[1]} \
            {input[0]} {input[1]}
        """

rule all:
     """Trim all reads with all supplied trimming parameters"""
     input:
         expand("./{trim_data}/{units}_R1_001.QT.CA.fastq.gz", units = config["units"], trim_data = config["trim_dir"]),
         expand("./{trim_data}/{units}_R2_001.QT.CA.fastq.gz", units = config["units"], trim_data = config["trim_dir"])
