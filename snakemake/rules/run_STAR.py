  __author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-22"

from snakemake.exceptions import MissingInputException

rule star_align:
    params:
        genomeDir = config["STAR"]["genomeDir"],
        runThreadN = config["STAR"]["runThreadN"]
    input:
        "fastq/subsets/{unit}_subset_R1_001.fastq.gz",
        "fastq/subsets/{unit}_subset_R2_001.fastq.gz"
    output:
        "{outdir}/STAR/subsets/{unit}.aligned.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {params.genomeDir} \
                 --readFilesIn {input[0]} {input[1]} \
                 --outSAMmode Full \
                 --outStd SAM \
                 --outSAMattributes Standard\
            | samtools view -b > {output[0]}
        """

rule star_align_full:
    params:
        genomeDir = config["STAR"]["genomeDir"],
        runThreadN = config["STAR"]["runThreadN"]
    input:
        "fastq/{unit}_R1_001.fastq.gz",
        "fastq/{unit}_R2_001.fastq.gz"
    output:
        "{outdir}/STAR/full/{unit}.aligned.sam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {params.genomeDir} \
                 --readFilesIn {input[0]} {input[1]} \
                 --readFilesCommand zcat \
                 --outTmpDir /home/skurscheid/tmp/{wildcards.unit} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 > {output[0]}
        """

rule run_htseq_count:
    params:
        htseq_dir = config["HTSeq_dir"],
        gtf = config["references"]["GTF"]
    input:
        rules.star_align_full.output
    output:
        "{outdir}/HTSeq/count/{unit}.txt"
    shell:
        """
            {params.htseq_dir}/htseq-count --format=bam \
                                          --order=pos \
                                          --stranded=yes \
                                          --type=exon \
                                          --idattr=gene_id \
                                          {input} \
                                          {params.gtf} \
                                          > {output}
        """
