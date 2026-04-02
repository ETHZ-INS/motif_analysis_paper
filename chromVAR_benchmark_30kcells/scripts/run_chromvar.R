library(SummarizedExperiment)
library(motifmatchr)
library(BiocParallel)
library(chromVARmotifs)
library(TFBSTools)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)

method <- snakemake@params[["method"]]
genome_pkg <- snakemake@params[["genome"]]
n_threads <- snakemake@threads
n_bg <- snakemake@params[["n_bg"]]

# Load Genome and Motifs
genome_obj <- get(genome_pkg)
motifs <- chromVARmotifs::encode_pwms

if (method == "original") {

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

} else if (method == "ArchR") {

    library(ArchR)
    addArchRThreads(threads = n_threads)
    addArchRGenome("mm10")
    proj <- loadArchRProject(snakemake@input[["arrow_dir"]])    
    proj <- addMotifAnnotations(proj, motifPWMs=motifs, force = TRUE)
    proj <- addBgdPeaks(proj, nIterations=n_bg, force = TRUE)    
    proj <- addDeviationsMatrix(proj, force=TRUE)
    
    dev <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")

    saveRDS(dev, snakemake@output[["dev"]])

} else if (method == "betterChromVAR") {

   library(betterChromVAR)
   rse <- readRDS(snakemake@input[["rse"]])
   rse <- addGCBias(rse, genome = genome_obj)
   motif_ix <- matchMotifs(motifs, rse, genome = genome_obj)

   dev <- betterChromVAR(rse, motif_ix, nthreads=as.integer(n_threads))
   saveRDS(dev, snakemake@output[["dev"]])

} else if (method == "scan_only") {

   library(chromVAR)
   rse <- readRDS(snakemake@input[["rse"]])
   rse <- addGCBias(rse, genome = genome_obj)
   motif_ix <- matchMotifs(motifs, rse, genome = genome_obj)
   saveRDS(list(), snakemake@output[["dev"]])
}
