java -Djava.io.tmpdir=/home/skurscheid/tmp \
     -Xmx40G \
     -jar /home/skurscheid/Bioinformatics/RNA-SeQC_v1.1.8.jar \
     -t /home/skurscheid/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf \
     -r /home/skurscheid/Data/References/Genomes/Mus_musculus/mm10_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa \
     -strat gc \
     -gc /home/skurscheid/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gc_content.txt \
     -s "10_28_hemi_PFC_ATTCCT|processed_data/GRCm38_ensembl84_ERCC/STAR/full/10_28_hemi_PFC_ATTCCT.aligned.bam|NA" \
     -gatkFlags --num_threads 16 \
     -o processed_data/GRCm38_ensembl84_ERCC/RNASeqQC/10_28_hemi_PFC_ATTCCT


java -Djava.io.tmpdir=/home/skurscheid/tmp \
      -Xmx36G \
      -jar /home/skurscheid/Bioinformatics/picard-tools-1.131/picard.jar CreateSequenceDictionary \
      REFERENCE=/home/skurscheid/Data/References/Genomes/Mus_musculus/mm10_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa \
      OUTPUT=/home/skurscheid/Data/References/Genomes/Mus_musculus/mm10_GRCm38/Mus_musculus.GRCm38.dna.toplevel.dict \
      GENOME_ASSEMBLY=mm10_GRCm38 \
      SPECIES=Mus_musculus


bowtie2 -k 200 --threads 8 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x /home/sebastian/Data/References/Transcriptomes/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.introns/trans_and_introns -1 processed_data/GRCm38_ensembl84/KMA_analysis/experiment/wt/wt1/wt1_R1.fastq.gz -2 processed_data/GRCm38_ensembl84/KMA_analysis/experiment/wt/wt1/wt1_R2.fastq.gz > ./test.bam
