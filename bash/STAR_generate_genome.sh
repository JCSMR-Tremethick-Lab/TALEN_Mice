STAR --runMode genomeGenerate\
     --runThreadN 96\
     --genomeDir ./STAR_Index_ERCC/\
     --genomeFastaFiles /home/skurscheid/Data/References/Genomes/Mus_musculus/mm10_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa\
     --sjdbGTFfile /home/skurscheid/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf\
     --sjdbOverhang 75

STAR --runThreadN 30 \
     --runMode genomeGenerate \
     --genomeDir ./STAR_Index_ERCC/ \
     --genomeFastaFiles /Data/References/Genomes/Mus_musculus/mm10_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa \
     --sjdbGTFfile /Data/References/Annotations/Mus_musculus/GRCm38_ensembl93/Mus_musculus.GRCm38.93.ERCC.gtf \
     --sjdbOverhang 75
