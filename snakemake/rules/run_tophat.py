__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-18"

from snakemake.exceptions import MissingInputException

rule tophat_align:
    params:
        bt_index = config["references"]["genome"],
        
    input:
    output:
    shell:

rule all:
    input:
