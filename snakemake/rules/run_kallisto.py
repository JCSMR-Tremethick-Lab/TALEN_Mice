__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

#configfile: "~/Development/JCSMR-Tremethick-Lab/Hodgkins-Lymphoma/snakemake/configs/config.json"

rule kallisto_quant:
    message:
        "Running kallisto..."
    params:
        ki = config["references"]["transcriptome_kallisto"],
        raw_data = config["raw_dir"],
        outdir = config["processed_dir"]
    input:
        "fastq/{unit}_R1_001.fastq.gz",
        "fastq/{unit}_R2_001.fastq.gz"
    output:
        "processed_data/{unit}"
    shell:
        """
            kallisto quant --index={params.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples=100 \
                           {input[0]} {input[1]}
        """

rule all:
    input:
        #expand("{output}/{sample}" , output = config["processed_dir"], sample = config["units"])
        expand("{outdir}/{unit}", outdir = config["processed_dir"], unit = config["units"])
        # expand("{outdir}/L12362715qia_S9", outdir = config["processed_dir"]),
        # expand("{outdir}/L12363-6-15_S7", outdir = config["processed_dir"])
