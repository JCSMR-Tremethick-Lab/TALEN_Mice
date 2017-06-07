snakemake --configfile config.json \
          --snakefile ~/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/workflows/TALEN_RNA-Seq.py all \
          --jobs 6 \
          --cluster "qsub -pe threads {cluster.threads} -q {cluster.queue} -l virtual_free={cluster.virtual_free} -l h_vmem={cluster.h_vmem}"\
          --cluster-config ~/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/configs/cluster.json

snakemake --configfile config_RNA-Seq_Brain.json \
          --snakefile ~/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/rules/run_RSeQC.py read_distribution\
          --jobs 6\
          -prn

snakemake --configfile config_RNA-Seq_Brain.json \
          --snakefile ~/Development/JCSMR-Tremethick-Lab/TALEN_Mice/snakemake/rules/run_RSeQC.py tin\
          --jobs 6\
          -prn
