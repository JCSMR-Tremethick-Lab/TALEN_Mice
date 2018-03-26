#!/bin/bash

# diff test on pooled replicate data
/home/sebastian/miniconda3/envs/irfinder/opt/irfinder-1.2.3/bin/analysisWithLowReplicates.pl \
  -A wt/pooled/IRFinder-IR-dir.txt wt/wt1/IRFinder-IR-dir.txt wt/wt2/IRFinder-IR-dir.txt wt/wt3/IRFinder-IR-dir.txt \
  -B hemi/pooled/IRFinder-IR-dir.txt hemi/hemi1/IRFinder-IR-dir.txt hemi/hemi2/IRFinder-IR-dir.txt hemi/hemi3/IRFinder-IR-dir.txt \
  > wt_vs_hemi.tab
