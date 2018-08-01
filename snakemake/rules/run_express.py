__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-03-26"

from snakemake.exceptions import MissingInputException

rule express_quantify:
    params:
        bt2_dir = home + "/miniconda3/bin",
        bt_index = config["references"]["trans_and_introns"]
    threads:
        1
    input:
        bam = rules.bowtie2_align.output,
        trans_and_introns = config["references"]["trans_and_introns_fasta"]
    output:
        file = "{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{unit}/express/results.xprs",
        dir = "{outdir}/{reference_version}/KMA_analysis/experiment/{condition}/{unit}/express/"
    shell:
        """
            {params.bt2_dir}/express {input.trans_and_introns} {input.bam} --output-dir {output.dir} --rf-stranded; touch {output.file}
        """
