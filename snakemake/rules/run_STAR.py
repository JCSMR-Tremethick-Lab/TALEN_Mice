  __author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-22"

from snakemake.exceptions import MissingInputException

wrapper_dir = "/home/skurscheid/Development/snakemake-wrappers/bio"

rule star_align_full:
    version:
        0.3
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
