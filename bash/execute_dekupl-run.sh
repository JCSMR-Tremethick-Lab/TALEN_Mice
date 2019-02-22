docker run --rm -v /home/skurscheid/Analysis/TALEN_Mice/snakemake/configs/config_H2AB3DSVM_dekupl.json:/dekupl/my-config.json \
       -v /mnt/raw/NB501086_0219_TSoboleva_JCSMR_RNAseq:/dekupl/data  -v /mnt/processed:/dekupl/results \
       transipedia/dekupl-run --configfile my-config.json  \
       -j8 -p
