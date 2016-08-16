__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

rule kallisto_quant_se:
    message:
        "Running kallisto with single end data..."
    params:
        raw_data = config["raw_dir"],
        outdir = config["processed_dir"],
        bootstraps = config["kallisto"]["bootstraps"],
    input:
        unit = lambda wildcards: config["units"][wildcards.unit],
        read1 = "./fastq/" + lambda wildcards: config["units"][wildcards.unit],
        ki = lambda wildcards: config["kallisto_index"][wildcards.reference_version]
    output:
        "{outdir}/{reference_version}/kallisto_se/{unit}"
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                           "{params.raw_data}/{input.read1}"
        """

rule kallisto_quant:
    message:
        "Running kallisto with paired end data..."
    params:
        raw_data = config["raw_dir"],
        outdir = config["processed_dir"],
        bootstraps = config["kallisto"]["bootstraps"],
    input:
        rules.cutadapt_pe.output,
        ki = lambda wildcards: config["kallisto_index"][wildcards.reference_version]
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

# rule kallisto_quant_pseudobam:
#     message:
#         "Running kallisto with pseudobam option..."
#     params:
#         raw_data = config["raw_dir"],
#         outdir = config["processed_dir"] + "\kallisto",
#         bootstraps = config["kallisto"]["bootstraps"],
#         ki=lambda wildcards: config["kallisto_index"][wildcards.ref],
#         sample=lambda wildcards: wildcards.unit
#     input:
#         "fastq/{unit}_R1_001.fastq.gz",
#         "fastq/{unit}_R2_001.fastq.gz"
#     output:
#         "{outdir}/{ref}/{unit}/",
#         "{outdir}/{ref}/{unit}/pseudobam/{unit}.bam"
#     shell:
#         """
#             mkdir -p {output}/pseudobam;
#             kallisto quant --index={params.ki} \
#                            --output-dir={output[0]} \
#                            --threads=1 \
#                            --bootstrap-samples=1 \
#                            --pseudobam \
#                            {input[0]} {input[1]} | \
#             samtools view -Sb - > {output[1]}
#         """
#
# rule bam_sort:
#     message:
#         "Sorting BAM files..."
#     params:
#         sample=lambda wildcards: wildcards.unit
#     input:
#         rules.kallisto_quant_pseudobam.output
#     output:
#         "{outdir}/{ref}/{unit}/pseudobam/{unit}.sorted.bam"
#     shell:
#         "samtools sort {input[1]} -T {wildcards.unit}.sorted -o {output}"
#
# rule bam_index:
#     message:
#         "Indexing sorted BAM files..."
#     input:
#         rules.bam_sort.output
#     output:
#         "{outdir}/{ref}/{unit}/pseudobam/{unit}.sorted.bai"
#     shell:
#         "samtools index {input} {output}"
#
# rule extract_fastq_data:
#     message:
#         "Extracting FASTQ records which mapped to transcripts of interest..."
#     params:
#         id_dir = config["id_dir"],
#         id_file = config["id_file"],
#         ids = getIDs(config["id_dir"] + "/" + config["id_file"])
#     input:
#         "processed_data/mm10.ens74.cdna.all_incl_h2a.Lap1_mutants/{unit}/pseudobam/{unit}.sorted.bam",
#         "fastq/{unit}_R1_001.fastq.gz",
#         "fastq/{unit}_R2_001.fastq.gz",
#         "processed_data/mm10.ens74.cdna.all_incl_h2a.Lap1_mutants/{unit}/pseudobam/{unit}.sorted.bai"
#     output:
#         "fastq/subsets/{unit}_subset_IDs.txt",
#         "fastq/subsets/{unit}_subset_R1_001.fastq.gz",
#         "fastq/subsets/{unit}_subset_R2_001.fastq.gz"
#     shell:
#         """
#             samtools view {input[0]} {params.ids} | cut -f 1 | sort | uniq > {output[0]}; \
#             seqtk subseq {input[1]} {output[0]} > {output[1]}; \
#             seqtk subseq {input[2]} {output[0]} > {output[2]}
#         """

# rule run_oases:
#     message:
#         "Running Oases assembler on cDNA reads..."
#     input:
#         "fastq/subsets/{unit}_subset_R1_001.fastq.gz",
#         "fastq/subsets/{unit}_subset_R2_001.fastq.gz"
#     output:
#         "assembly/{unit}"
#     shell:
#         """
#             oases_pipeline.py -s 2 -o {output} \
#             -d '-fastq -shortPaired -separate {input[0]} {input[1]}' \
#             -p '-ins_length 200 -min_trans_lgth 100'
#         """
