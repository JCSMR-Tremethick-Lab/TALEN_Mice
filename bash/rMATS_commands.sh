# starting with pre-mapped BAMs
source activate py27
python ~/Bioinformatics/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --b1 wt.txt --b2 hemi.txt --gtf /home/sebastian/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf --od rMATS -t paired --readLength 76 --cstat 0.0001 --libType fr-firststrand --nthread 32 --tstat 32

# starting from scratch - trimmed FASTQ
python ~/Bioinformatics/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --s1 wt_fastq.txt --s2 hemi_fastq.txt\
                                                                    --gtf /home/sebastian/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf\
                                                                    --bi /home/sebastian/Data/References/Transcriptomes/Mus_musculus/GRCm38_ensembl84/STAR_Index \
                                                                    --od rMATS_fastq\
                                                                    -t paired\
                                                                    --readLength 76\
                                                                    --cstat 0.0001\
                                                                    --libType fr-firststrand\
                                                                    --nthread 32\
                                                                    --tstat 32

# brain - HIPPO
python ~/Bioinformatics/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --s1 wt_HIPPO.txt --s2 hemi_HIPPO.txt\
                                                                    --gtf /home/sebastian/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf\
                                                                    --bi /home/sebastian/Data/References/Transcriptomes/Mus_musculus/GRCm38_ensembl84/STAR_Index \
                                                                    --od rMATS_HIPPO_wt_vs_hemi\
                                                                    -t paired\
                                                                    --readLength 76\
                                                                    --cstat 0.0001\
                                                                    --libType fr-firststrand\
                                                                    --nthread 32\
                                                                    --tstat 32

# brain - OB
python ~/Bioinformatics/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --s1 wt_OB.txt --s2 hemi_OB.txt\
                                                                    --gtf /home/sebastian/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf\
                                                                    --bi /home/sebastian/Data/References/Transcriptomes/Mus_musculus/GRCm38_ensembl84/STAR_Index \
                                                                    --od rMATS_OB_wt_vs_hemi\
                                                                    -t paired\
                                                                    --readLength 76\
                                                                    --cstat 0.0001\
                                                                    --libType fr-firststrand\
                                                                    --nthread 32\
                                                                    --tstat 32

# brain - PFC
python ~/Bioinformatics/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --b1 wt_PFC.txt --b2 hemi_PFC.txt\
                                                                    --gtf ~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf\
                                                                    --od rMATS_PFC_wt_vs_hemi\
                                                                    -t paired\
                                                                    --readLength 76\
                                                                    --cstat 0.0001\
                                                                    --libType fr-firststrand\
                                                                    --nthread 48\
                                                                    --tstat 48

# Fear conditioning experiment
#
cd /home/sebastian/Data/Tremethick/TALENs/RNA-Seq/Mus_musculus_brain_experiment_2/processed_data/GRCm38_ensembl84/rMATS

python ~/Bioinformatics/rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --b1 pft_wt_naive.txt --b2 pft_wt_naive.txt\
                                                                    --gtf ~/Data/References/Annotations/Mus_musculus/GRCm38_ensembl84/Mus_musculus.GRCm38.84.gtf\
                                                                    --od pfc_naive_wt_vs_hemi\
                                                                    -t paired\
                                                                    --readLength 76\
                                                                    --cstat 0.0001\
                                                                    --libType fr-firststrand\
                                                                    --nthread 32\
                                                                    --tstat 32
