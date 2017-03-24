  __author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-22"

from snakemake.exceptions import MissingInputException

wrapper_dir = "/home/skurscheid/Development/snakemake-wrappers/bio"

def getGroups(wildcards):
    cond1 = []
    cond2 = []
    c1 = wildcards.condition.split("_vs_")[0]
    c2 = wildcards.condition.split("_vs_")[1]
    for i in config["Groups"][wildcards.tissue][c1]:
        cond1.append(wildcards.outdir + "/" + wildcards.reference_version + "/STAR/full/" + i + ".aligned.bam")
    for i in config["Groups"][wildcards.tissue][c2]:
        cond2.append(wildcards.outdir + "/" + wildcards.reference_version + "/STAR/full/" + i + ".aligned.bam")
    cond1 = ",".join(cond1)
    cond2 = ",".join(cond2)
    return(cond1, cond2)


rule star_align_full:
    version:
        0.5
    threads:
        lambda wildcards: int(str(config["STAR"]["runThreadN"]).strip("['']"))
    params:
            tmp_dir = home + "/tmp"
    input:
        rules.cutadapt_pe.output,
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

rule run_dexseq_count:
    version:
        0.1
    params:
        dexseq_dir = home + config["DEXSeq_dir"],
        dex_gtf = config["references"]["DEX_GTF"]
    input:
        bam = "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
        index = "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai"
    output:
        "{outdir}/{reference_version}/DEXSeq/count/{unit}.txt"
    shell:
        """
            python {params.dexseq_dir}/dexseq_count.py --format=bam \
                                                       --paired=yes \
                                                       --order=pos \
                                                       --stranded=reverse \
                                                       {params.dex_gtf} \
                                                       {input.bam} \
                                                       {output}
        """


rule collect_insert_size_metrics:
    version:
        0.1
    params:
        sampling = config["Picard"]["sampling"],
        jar_file = home + "/bin/picard.jar",
        tmp_dir = home + "/tmp"
    input:
        rules.star_align_full.output
    output:
        txt = "{outdir}/{reference_version}/PICARD/insert_size_metrics/{unit}.insert_size_metrics.txt",
        pdf = "{outdir}/{reference_version}/PICARD/insert_size_metrics/{unit}.insert_size_metrics.pdf"
    shell:
        """
            java -Djava.io.tmpdir= {params.tmp_dir} \
            -Xmx36G \
            -jar {params.jar_file} CollectInsertSizeMetrics \
            I={input} \
            O={output.txt} \
            H={output.pdf} \
            M=0.2
        """

rule run_rMats:
    version:
        0.1
    params:
        gtf = config["references"]["GTF"],
        bin = home + "/Bioinformatics/rMATS.3.2.5/RNASeq-MATS.py"
    input:
        getGroups
    output:
        "{outdir}/{reference_version}/rMATS/{tissue}/{condition}"
    shell:
        """
            python {params.bin} -b1 {input[0]} \
                                -b2 {input[1]} \
                                -gtf {params.gtf} \
                                -t paired \
                                -len 76 \
                                -analysis U \
                                -libType fr-firststrand \
                                -o {output}
        """
