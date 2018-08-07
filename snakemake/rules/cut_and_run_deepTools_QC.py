__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-07"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools QC/QC on ChIP-Seq data
For usage, include this in your workflow.
"""

# global functions
def get_sample_labels(wildcards):
    sl = []
    runIDs = config["samples"][wildcards.assayType].keys()
    for i in runIDs:
        for k in config["samples"][wildcards.assayType][i].keys():
            sl.append(k)
    return(sl)

REF_VERSION = "GRCm38_ensembl93"
RUN_ID = "180731_NB501086_0217_CutandRun_Tanya"

rule multiBamSummary:
    version:
        "2"
    params:
        deepTools_dir = home + config["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        labels = get_sample_labels
    threads:
        32
    input:
        expand("{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.{suffix}",
               assayType = "CutRun",
               runID = RUN_ID,
               reference_version = REF_VERSION,
               library = config["CutRun"]["samples"][RUN_ID],
               suffix = ["bam"])
    output:
        npz = "{assayType}/deepTools/multiBamSummary/{reference_version}/{runID}/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input} \
                                                        --labels {params.labels} \
                                                        --numberOfProcessors {threads} \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz}
        """


rule plotCorrelation_heatmap:
    version:
        "2"
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = "Correlation heatmap - read counts"
    input:
        npz = rules.multiBamSummary.output.npz
    output:
        png = "{assayType}/deepTools/plotCorrelation/{reference_version}/{runID}/heatmap_SpearmanCorr_readCounts.png",
        tab = "{assayType}/deepTools/plotCorrelation/{reference_version}/{runID}/heatmap_SpearmanCorr_readCounts.tab"
    shell:
        """
            {params.deepTools_dir}/plotCorrelation --corData {input.npz} \
                                                   --corMethod spearman \
                                                   --skipZeros \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --whatToPlot heatmap \
                                                   --colorMap RdYlBu \
                                                   --plotNumbers \
                                                   -o {output.png} \
                                                   --outFileCorMatrix {output.tab}
        """

rule plotPCA:
    version:
        "2"
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = "PCA - read counts"
    input:
        npz = rules.multiBamSummary.output.npz
    output:
        png = "{assayType}/deepTools/plotPCA/{reference_version}/{runID}/PCA_readCounts.png"
    shell:
        """
            {params.deepTools_dir}/plotPCA --corData {input.npz} \
                                           --plotFile {output.png} \
                                           --plotTitle "{params.plotTitle}"
        """

rule bamPEFragmentSize:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fragment size",
        labels = get_sample_labels
    threads:
        32
    input:
        expand("{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.{suffix}",
               assayType = "CutRun",
               runID = RUN_ID,
               reference_version = REF_VERSION,
               library = config["CutRun"]["samples"][RUN_ID],
               suffix = ["bam"])
    output:
        "{assayType}/deepTools/bamPEFragmentSize/{reference_version}/{runID}/histogram_duplicates_marked.png"
    shell:
        """
            {params.deepTools_dir}/bamPEFragmentSize --bamfiles {input} \
                                                     --samplesLabel {params.labels} \
                                                     --numberOfProcessors {threads} \
                                                     --histogram {output}
        """


rule plotFingerprint:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = "BAM PE fingerprint",
        labels = get_sample_labels
    threads:
        32
    input:
        expand("{assayType}/samtools/rmdup/{reference_version}/{runID}/{library}.{suffix}",
               assayType = "CutRun",
               runID = RUN_ID,
               reference_version = REF_VERSION,
               library = config["CutRun"]["samples"][RUN_ID],
               suffix = ["bam"])
    output:
        "{assayType}/deepTools/plotFingerprint/{reference_version}/{runID}/fingerprints_duplicates_marked.png"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                                   --numberOfProcessors {threads} \
                                                   --centerReads \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --labels {params.labels} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """

rule all:
    input:
        expand("{assayType}/deepTools/plotFingerprint/{reference_version}/{runID}/fingerprints_duplicates_marked.png",
                assayType = "CutRun",
                reference_version = REF_VERSION,
                runID = RUN_ID),
        expand("{assayType}/deepTools/bamPEFragmentSize/{reference_version}/{runID}/histogram_duplicates_marked.png",
                assayType = "CutRun",
                reference_version = REF_VERSION,
                runID = RUN_ID),
        expand("{assayType}/deepTools/plotPCA/{reference_version}/{runID}/PCA_readCounts.png",
                assayType = "CutRun",
                reference_version = REF_VERSION,
                runID = RUN_ID),
        expand("{assayType}/deepTools/plotCorrelation/{reference_version}/{runID}/heatmap_SpearmanCorr_readCounts.{suffix}",
                assayType = "CutRun",
                reference_version = REF_VERSION,
                runID = RUN_ID,
                suffix = ["png", "tab"])
