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
            mkdir -p {output}/pseudobam;
            kallisto quant --index={params.ki} \
                           --output-dir={output} \
                           --threads=1 \
                           --bootstrap-samples=1 \
                           --pseudobam \
                           {input[0]} {input[1]} | \
            samtools view -Sb - > {output}/pseudobam/{params.sample}.bam
        """

rule bam_sort:
    params:
        sample=lambda wildcards: wildcards.unit
    input:
        "{kallisto_quant_pseudobam.output}/pseudobam/{params.sample}.bam"
    output:
        "{kallisto_quant_pseudobam.output}/pseudobam/{params.sample}.sorted.bam"
    shell:
        "samtools sort {input} -T {wildcards.unit}.sorted -o {output}"

rule bam_index:
    input:
        "{bam_sort.output}"
    output:
        "{bam_sort.output}.bai"
    shell:
        "samtools index {input} {output}"

rule all:
    input:
        #expand("{output}/{sample}" , output = config["processed_dir"], sample = config["units"])
        #expand("{outdir}/{ref}/{unit}", outdir = config["processed_dir"], unit = config["units"], ref = config["kallisto_index"]),
        expand("{outdir}/{ref}/{unit}/pseudobam/{sample}.bai", outdir = config["processed_dir"], unit = config["units"], ref = config["kallisto_index"], sample = config["units"][0])
        # expand("{outdir}/L12362715qia_S9", outdir = config["processed_dir"]),
        # expand("{outdir}/L12363-6-15_S7", outdir = config["processed_dir"])
