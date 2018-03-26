__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-03-26"

from snakemake.exceptions import MissingInputException

rule tophat_align:
    params:
        bt2_dir = home + "/miniconda3/bin"
    input:
        left = "{wildcards.outdir}/{wildcards.reference_version}/KMA_analysis/experiment/{wildcards.condition}/{wildcards.unit}",
        bt_index = config["references"]["trans_and_introns"]
    output:
        "{outdir}/{reference_version}/bowtie2/{unit}"
    shell:
        """
            {params.bt2_dir}/bowtie2 -k 200 \
                                     --rdg 6,5\
                                     --rfg 6,5\
                                     --score-min L,-.6,-.4\
                                     -X {input.bt_index}\
                                     -1 {input.left} -2 {input.right} | samtools view -Sb - > {output}
        """
