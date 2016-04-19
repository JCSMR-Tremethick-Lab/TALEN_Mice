__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-18"

from snakemake.exceptions import MissingInputException

rule tophat_align:
    params:
        bt_index = config["references"]["genome"],
        ref_trome = config["tophat"]["GTF"],
        del_length = config["tophat"]["max-deletion-length"],
        lib_type = config["tophat"]["library-type"],
        sample=lambda wildcards: wildcards.unit
    input:
        "fastq/filtered/{unit}_R1_001.fastq.gz",
        "fastq/filtered/{unit}_R2_001.fastq.gz"
    output:
        "processed_data/tophat2/{unit}"
    shell:
        """
            tophat2 --output-dir {ouput} \
                    --GTF {params.ref_trome} \
                    --max-deletion-length {params.del_length} \
                    --library-type {params.lib_type} \
                    {params.bt_index} \
                    {input[0]} {input[1]}
        """

rule all:
    input:
        expand("{outdir}/tophat2/{unit}/accepted_hits.bam", outdir = config["processed_dir"], unit = config["units"])
