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

library(betterChromVAR)
rse <- readRDS(snakemake@input[["rse"]])
rse <- addGCBias(rse, genome = genome_obj)
motif_ix <- matchMotifs(motifs, rse, genome = genome_obj)

bg <- getBackgroundKNN(rse)
dev <- computeDeviationsFromKNN(rse, cBg=bg, motif_ix, l=1L)
saveRDS(dev, snakemake@output[["dev"]])
