__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-08-19"

from snakemake.exceptions import MissingInputException

rule run_rnaseq_qc:
    message:
        "Running RNA-Seq QC"
    params:
    input:
        rules.star_align_full.output,
        genome = config["references"]["genome"],
        gtf = config["references"]["GTF"],
        gc = config["references"]["transcript_GC"]
    output:
        "{outdir}/{reference_version}/RNASeqQC/{unit}/"
    shell:
        """
            java -Djava.io.tmpdir=/home/skurscheid/tmp \
            -Xmx36G \
            -jar /home/skurscheid/Bioinformatics/RNA-SeQC_v1.1.8.jar \
            -t {input.gtf} \
            -r {input.genome} \
            -strat gc \
            -gc {input.gc} \
            -s {input[0]}
        """
