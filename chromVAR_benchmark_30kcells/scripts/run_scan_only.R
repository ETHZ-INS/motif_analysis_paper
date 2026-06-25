library(SummarizedExperiment)
library(motifmatchr)
library(BiocParallel)
library(chromVARmotifs)
library(TFBSTools)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)

genome_pkg <- snakemake@params[["genome"]]
n_threads <- snakemake@threads
n_bg <- snakemake@params[["n_bg"]]

# Load Genome and Motifs
genome_obj <- get(genome_pkg)
motifs <- chromVARmotifs::encode_pwms

library(chromVAR)
rse <- readRDS(snakemake@input[["rse"]])
rse <- addGCBias(rse, genome = genome_obj)
motif_ix <- matchMotifs(motifs, rse, genome = genome_obj)
saveRDS(list(), snakemake@output[["dev"]])
