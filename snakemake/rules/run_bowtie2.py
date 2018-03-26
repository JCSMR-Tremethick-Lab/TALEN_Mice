__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-03-26"

from snakemake.exceptions import MissingInputException

rule tophat_align:
    params:
        bt2_dir = home + "/miniconda3/bin"
    threads:
        8
    input:
        left = "{wildcards.outdir}/{wildcards.reference_version}/KMA_analysis/experiment/{wildcards.condition}/{wildcards.unit}/{wildcards.unit}_R1.fastq.gz",
        right = "{wildcards.outdir}/{wildcards.reference_version}/KMA_analysis/experiment/{wildcards.condition}/{wildcards.unit}/{wildcards.unit}_R2.fastq.gz",
        bt_index = config["references"]["trans_and_introns"]
    output:
        "{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{unit}/hits.bam"
    shell:
        """
            {params.bt2_dir}/bowtie2 -k 200 \
                                     --threads {threads}\
                                     --rdg 6,5\
                                     --rfg 6,5\
                                     --score-min L,-.6,-.4\
                                     -X {input.bt_index}\
                                     -1 {input.left} -2 {input.right} | samtools view -Sb - > {output}
        """
