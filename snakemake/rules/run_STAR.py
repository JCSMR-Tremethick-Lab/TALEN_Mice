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
        cond1.append("./" + wildcards.outdir + "/" + wildcards.reference_version + "/STAR/full/" + i + ".aligned.bam")
    for i in config["Groups"][wildcards.tissue][c2]:
        cond2.append("./" + wildcards.outdir + "/" + wildcards.reference_version + "/STAR/full/" + i + ".aligned.bam")
    #cond1 = ",".join(cond1)
    #cond2 = ",".join(cond2)
    return(cond1, cond2)


rule star_align_full:
    version:
        0.4
    params:
        runThreadN = config["STAR"]["runThreadN"]
    input:
        rules.cutadapt_pe.output,
        index = lambda wildcards: config["STAR"][wildcards.reference_version]
    output:
        "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {input.index} \
                 --readFilesIn {input[0]} {input[1]} \
                 --readFilesCommand zcat \
                 --outTmpDir /home/skurscheid/tmp/{wildcards.unit} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 > {output[0]}
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

rule run_htseq_count:
    version:
        0.3
    params:
        htseq_dir = config["HTSeq_dir"],
        gtf = config["references"]["GTF"]
    input:
        bam = "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
        index = "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai"
    output:
        "{outdir}/{reference_version}/HTSeq/count/{unit}.txt"
    shell:
        """
            {params.htseq_dir}/htseq-count --format=bam \
                                           --order=pos \
                                           --stranded=reverse \
                                           --type=exon \
                                           --idattr=gene_id \
                                           --order=pos \
                                           {input.bam} \
                                           {params.gtf} \
                                           > {output}
        """

rule collect_insert_size_metrics:
    version:
        0.1
    params:
        sampling = config["Picard"]["sampling"]
    input:
        rules.star_align_full.output
    output:
        txt = "{outdir}/{reference_version}/PICARD/insert_size_metrics/{unit}.insert_size_metrics.txt",
        pdf = "{outdir}/{reference_version}/PICARD/insert_size_metrics/{unit}.insert_size_metrics.pdf"
    shell:
        """
            java -Djava.io.tmpdir=/home/skurscheid/tmp \
            -Xmx36G \
            -jar /home/skurscheid/Bioinformatics/picard-tools-1.131/picard.jar CollectInsertSizeMetrics \
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
        bin = "/home/skurscheid/Bioinformatics/rMATS.3.2.2.beta/RNASeq-MATS.py"
    input:
        getGroups
    output:
        "{outdir}/{reference_version}/rMATS/{tissue}/{condition}"
    shell:
        """
            python {params.bin} -b1 ",".join({input[0]}) \
                                -b2 ",".join({input[1]}) \
                                -gtf {params.gtf} \
                                -t paired \
                                -len 76 \
                                -analysis P \
                                -libType fr-firststrand \
                                -o {output}
        """
