__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-18"

from snakemake.exceptions import MissingInputException

rule tophat_align:
    params:
        bt_index = config["references"]["genome"],
        ref_trome = config["tophat"]["GTF"],
        del_length = config["tophat"]["max-deletion-length"],
        intron_length = config["tophat"]["min-intron-length"],
        lib_type = config["tophat"]["library-type"]
    input:
        "fastq/subsets/{unit}_subset_R1_001.fastq.gz",
        "fastq/subsets/{unit}_subset_R2_001.fastq.gz"
    output:
        "processed_data/tophat2/{unit}"
    shell:
        """
            tophat2 --output-dir {output} \
                    --GTF {params.ref_trome} \
                    --max-deletion-length {params.del_length} \
                    --library-type {params.lib_type} \
                    --min-intron-length {params.intron_length} \
                    {params.bt_index} \
                    {input[0]} {input[1]}
        """

rule all:
    input:
        expand("{outdir}/tophat2/{unit}", outdir = config["processed_dir"], unit = config["units"])
        # expand("{outdir}/tophat2/NMG3-60hemi_S1", outdir = config["processed_dir"]),
        # expand("{outdir}/tophat2/NMG3-62wt_S2", outdir = config["processed_dir"]),
        # expand("{outdir}/tophat2/NMG3-74wt_S3", outdir = config["processed_dir"]),
        # expand("{outdir}/tophat2/NMG3-75hemi_S4", outdir = config["processed_dir"]),
        # expand("{outdir}/tophat2/NMG3-76wt_S5", outdir = config["processed_dir"]),
        # expand("{outdir}/tophat2/NMG3-77hemi_S6", outdir = config["processed_dir"])