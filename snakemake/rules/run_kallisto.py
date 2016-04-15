__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

#configfile: "~/Development/JCSMR-Tremethick-Lab/Hodgkins-Lymphoma/snakemake/configs/config.json"

rule kallisto_quant:
    message:
        "Running kallisto..."
    params:
        ki = config["references"]["mm10.ens74.cdna.all"],
        raw_data = config["raw_dir"],
        outdir = config["processed_dir"],
        bootstraps = config["kallisto"]["bootstraps"]
    input:
        "fastq/{unit}_R1_001.fastq.gz",
        "fastq/{unit}_R2_001.fastq.gz",
        "/home/skurscheid/Data/RefGenomes/Mus_musculus/mm10_GRCm38/{ref}""
    output:
        "processed_data/{unit}/{ref}"
    shell:
        """
            kallisto quant --index={input[3]} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                           {input[0]} {input[1]}
        """

rule all:
    input:
        #expand("{output}/{sample}" , output = config["processed_dir"], sample = config["units"])
        expand("{outdir}/{unit}/{ref}", outdir = config["processed_dir"], unit = config["units"], ref = config["references"])
        # expand("{outdir}/L12362715qia_S9", outdir = config["processed_dir"]),
        # expand("{outdir}/L12363-6-15_S7", outdir = config["processed_dir"])
