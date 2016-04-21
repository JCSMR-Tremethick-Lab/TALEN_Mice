__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

#configfile: "~/Development/JCSMR-Tremethick-Lab/Hodgkins-Lymphoma/snakemake/configs/config.json"

# REF_TO_PATH = {
#     path: paths in config["kallisto_index"].items()
#     for path in paths}

def getIDs( file ):
    fo = open(file, "r")
    line = [x.strip() for x in fo.readlines()]
    line = ' '.join(line)
    return line

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
        "processed_data/{ref}/{unit}/"
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
        "{outdir}/{ref}/{unit}/pseudobam/{unit}.bam"
    output:
        "{outdir}/{ref}/{unit}/pseudobam/{unit}.sorted.bam"
    shell:
        "samtools sort {input} -T {wildcards.unit}.sorted -o {output}"

rule bam_index:
    input:
        "{outdir}/{ref}/{unit}/pseudobam/{unit}.sorted.bam"
    output:
        "{outdir}/{ref}/{unit}/pseudobam/{unit}.sorted.bai"
    shell:
        "samtools index {input} {output}"

rule extract_fastq_data:
    params:
        id_dir = config["id_dir"],
        id_file = config["id_file"],
        ids = getIDs(config["id_dir"] + "/" + config["id_file"])
    input:
        "{outdir}/{ref}/{unit}/pseudobam/{unit}.sorted.bam",
        "fastq/{unit}_R1_001.fastq.gz",
        "fastq/{unit}_R2_001.fastq.gz"
    output:
        "fastq/subsets/{unit}_subset_IDs.txt",
        "fastq/subsets/{unit}_subset_R1_001.fastq.gz",
        "fastq/subsets/{unit}_subset_R2_001.fastq.gz"
    shell:
        """
            samtools view {input[0]} {params.ids} | cut -f 1 | sort | uniq > {output[0]}; \
            seqtk subseq {input[1]} {output[0]} > {output[1]}; \
            seqtk subseq {input[2]} {output[0]} > {output[2]}
        """




rule all:
    input:
        expand("{outdir}/mm10.ens74.cdna.all_incl_h2a.Lap1_mutants/NMG3-60hemi_S1/pseudobam/NMG3-60hemi_S1.sorted.bai", outdir = config["processed_dir"]),
        expand("{outdir}/mm10.ens74.cdna.all_incl_h2a.Lap1_mutants/NMG3-62wt_S2/pseudobam/NMG3-62wt_S2.sorted.bai", outdir = config["processed_dir"]),
        expand("{outdir}/mm10.ens74.cdna.all_incl_h2a.Lap1_mutants/NMG3-74wt_S3/pseudobam/NMG3-74wt_S3.sorted.bai", outdir = config["processed_dir"]),
        expand("{outdir}/mm10.ens74.cdna.all_incl_h2a.Lap1_mutants/NMG3-75hemi_S4/pseudobam/NMG3-75hemi_S4.sorted.bai", outdir = config["processed_dir"]),
        expand("{outdir}/mm10.ens74.cdna.all_incl_h2a.Lap1_mutants/NMG3-76wt_S5/pseudobam/NMG3-76wt_S5.sorted.bai", outdir = config["processed_dir"]),
        expand("{outdir}/mm10.ens74.cdna.all_incl_h2a.Lap1_mutants/NMG3-77hemi_S6/pseudobam/NMG3-77hemi_S6.sorted.bai", outdir = config["processed_dir"]),
        expand("fastq/subsets/{unit}_subset_IDs.txt", unit = config["units"]),
        expand("fastq/subsets/{unit}_subset_R1_001.fastq.gz", unit = config["units"]),
        expand("fastq/subsets/{unit}_subset_R2_001.fastq.gz", unit = config["units"])
