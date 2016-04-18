__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

#configfile: "~/Development/JCSMR-Tremethick-Lab/Hodgkins-Lymphoma/snakemake/configs/config.json"

# REF_TO_PATH = {
#     path: paths in config["kallisto_index"].items()
#     for path in paths}

# rule kallisto_quant:
#     message:
#         "Running kallisto..."
#     params:
#         raw_data = config["raw_dir"],
#         outdir = config["processed_dir"],
#         bootstraps = config["kallisto"]["bootstraps"],
#         ki=lambda wildcards: config["kallisto_index"][wildcards.ref]
#     input:
#         "fastq/{unit}_R1_001.fastq.gz",
#         "fastq/{unit}_R2_001.fastq.gz"
#     output:
#         "processed_data/{ref}/{unit}"
#     shell:
#         """
#             kallisto quant --index={params.ki} \
#                            --output-dir={output} \
#                            --threads=4 \
#                            --bootstrap-samples={params.bootstraps} \
#                            {input[0]} {input[1]}
#         """

rule kallisto_quant_pseudobam:
    message:
        "Running kallisto with pseudobam option..."
    params:
        raw_data = config["raw_dir"],
        outdir = config["processed_dir"],
        bootstraps = config["kallisto"]["bootstraps"],
        ki=lambda wildcards: config["kallisto_index"][wildcards.ref],
        sample=lambda wildcards: wildcards.unit
    input:
        "fastq/{unit}_R1_001.fastq.gz",
        "fastq/{unit}_R2_001.fastq.gz"
    output:
        "processed_data/{ref}/{unit}"
    shell:
        """
            kallisto quant --index={params.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                           --pseudobam \
                           {input[0]} {input[1]} | \
            samtools view -Sb - > {params.sample}.bam
        """

rule all:
    input:
        #expand("{output}/{sample}" , output = config["processed_dir"], sample = config["units"])
        expand("{outdir}/{ref}/{unit}", outdir = config["processed_dir"], unit = config["units"], ref = config["kallisto_index"])
        # expand("{outdir}/L12362715qia_S9", outdir = config["processed_dir"]),
        # expand("{outdir}/L12363-6-15_S7", outdir = config["processed_dir"])
