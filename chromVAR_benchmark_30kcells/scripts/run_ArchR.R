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

    library(ArchR)
    addArchRThreads(threads = n_threads)
    addArchRGenome("mm10")
    proj <- loadArchRProject(snakemake@input[["arrow_dir"]])    
    proj <- addMotifAnnotations(proj, motifPWMs=motifs, force = TRUE)
    proj <- addBgdPeaks(proj, nIterations=n_bg, force = TRUE)    
    proj <- addDeviationsMatrix(proj, force=TRUE)
    
    dev <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")

    saveRDS(dev, snakemake@output[["dev"]])

