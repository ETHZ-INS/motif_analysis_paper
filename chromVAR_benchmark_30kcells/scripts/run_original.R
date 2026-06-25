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
    if(n_threads==1){
      register(SerialParam())
    }else{
      register(MulticoreParam(n_threads))
    }
    
    rse <- addGCBias(rse, genome = genome_obj)
    motif_ix <- matchMotifs(motifs, rse, genome = genome_obj)
    
    bg <- getBackgroundPeaks(rse, niterations=n_bg)
    dev <- computeDeviations(object = rse, annotations = motif_ix, background_peaks=bg)
    saveRDS(dev, snakemake@output[["dev"]])

