__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-03-26"

from snakemake.exceptions import MissingInputException

rule bowtie2_align:
    params:
        bt2_dir = home + "/miniconda3/bin",
        bt_index = config["references"]["trans_and_introns"]
    threads:
        8
    input:
        left = "{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{unit}/{unit}_R1.fastq.gz",
        right = "{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{unit}/{unit}_R2.fastq.gz"
    output:
        "{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{unit}/hits.bam"
    shell:
        """
            {params.bt2_dir}/bowtie2 -k 200 \
                                     --threads {threads}\
                                     --rdg 6,5\
                                     --rfg 6,5\
                                     --score-min L,-.6,-.4\
                                     -X {params.bt_index}\
                                     -1 {input.left} -2 {input.right} | samtools view -Sb - > {output}
        """
