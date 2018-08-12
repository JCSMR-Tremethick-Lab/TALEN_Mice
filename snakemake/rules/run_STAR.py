  __author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-22"

from snakemake.exceptions import MissingInputException

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

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
        "1"
    threads:
        8
    params:
        tmp_dir = home + "/tmp"
    input:
        trimmed_read1 = rules.run_fastp.output.trimmed_read1,
        trimmed_read2 = rules.run_fastp.output.trimmed_read2,
        index = lambda wildcards: config["STAR"][wildcards["reference_version"]]["index"],
        gtf = lambda wildcards: config["STAR"][wildcards["reference_version"]]["GTF"]
    output:
        "{assayType}/STAR/full/{reference_version}/{runID}/{library}.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {threads} \
                 --genomeDir {input.index} \
                 --readFilesIn {input.trimmed_read1} {input.trimmed_read1} \
                 --readFilesCommand zcat \
                 --outTmpDir {params.tmp_dir}/{wildcards.library} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 > {output}
        """

rule star_align_rMATs:
    version:
        "1"
    threads:
        16
    params:
        tempDir = home + "/tmp/",
        tophatAnchor = "2"
    input:
        trimmed_read1 = rules.run_fastp.output.trimmed_read1,
        trimmed_read2 = rules.run_fastp.output.trimmed_read2,
        index = lambda wildcards: config["STAR"][wildcards["reference_version"]]["index"],
        gtf = lambda wildcards: config["STAR"][wildcards["reference_version"]]["GTF"]
    output:
        "{assayType}/rMATs/BAMs/{reference_version}/{runID}/{library}.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {threads} \
                 --genomeDir {input.index} \
                 --readFilesIn {input.trimmed_read1} {input.trimmed_read2} \
                 --readFilesCommand zcat \
                 --outTmpDir {params.tempDir}{wildcards.library} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMstrandField intronMotif\
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 --chimSegmentMin 2\
                 --outFilterMismatchNmax 3\
                 --alignSJDBoverhangMin {params.tophatAnchor}\
                 --alignIntronMax 299999\
                 --sjdbGTFfile {input.gtf}\
                 > {output[0]}
        """


rule bam_index_STAR_output:
    version:
        "1"
    input:
        "{assayType}/{tool}/{subcommand}/{reference_version}/{runID}/{library}.bam"
    output:
        "{assayType}/{tool}/{subcommand}/{reference_version}/{runID}/{library}.bam.bai"
    wrapper:
        "file://" + wrapper_dir + "/samtools/index/wrapper.py"


rule create_bigwig_from_bam_1xgenome:
    version:
        "1"
    threads:
        16
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        outFileFormat = "bigwig",
        binSize = 10,
        normalizeUsing = "RPGC",
        effectiveGenomeSize = 2150570000
    input:
        bam = "{assayType}/{tool}/{subcommand}/{reference_version}/{runID}/{library}.bam",
        index = "{assayType}/{tool}/{subcommand}/{reference_version}/{runID}/{library}.bam.bai"
    output:
        bigwig = "{assayType}/deepTools/bamCoverage/{reference_version}/{runID}/{tool}/{subcommand}/{library}.bw"
    shell:
        """
            {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                               --outFileName {output} \
                                               --outFileFormat {params.outFileFormat} \
                                               --numberOfProcessors {threads} \
                                               --normalizeUsing {params.normalizeUsing} \
                                               --effectiveGenomeSize {params.effectiveGenomeSize} \
                                               --binSize {params.binSize}
        """
