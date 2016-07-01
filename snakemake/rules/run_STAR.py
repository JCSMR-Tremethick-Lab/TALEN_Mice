  __author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-22"

from snakemake.exceptions import MissingInputException

wrapper_dir = "/home/skurscheid/Development/snakemake-wrappers/bio"

rule star_align:
    version:
        0.1
    params:
        runThreadN = config["STAR"]["runThreadN"]
    input:
        genomeDir = lambda wildcards: config["STAR"][wildcards.reference_version],
        "fastq/subsets/{unit}_subset_R1_001.fastq.gz",
        "fastq/subsets/{unit}_subset_R2_001.fastq.gz"
    output:
        "{outdir}/{reference_version}/STAR/subsets/{unit}.aligned.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {input.genomeDir} \
                 --readFilesIn {input[0]} {input[1]} \
                 --outSAMmode Full \
                 --outStd SAM \
                 --outSAMattributes Standard\
            | samtools view -b > {output[0]}
        """

rule star_align_full:
    version:
        0.3
    params:
        index = config["STAR"]["genomeDir"],
        runThreadN = config["STAR"]["runThreadN"]
    input:
        rules.cutadapt_pe.output
    output:
        "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {params.index} \
                 --readFilesIn {input[0]} {input[1]} \
                 --readFilesCommand zcat \
                 --outTmpDir /home/skurscheid/tmp/{wildcards.unit} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
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
        0.2
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
                                          --stranded=yes \
                                          --type=exon \
                                          --idattr=gene_id \
                                          {input.bam} \
                                          {params.gtf} \
                                          > {output}
        """
