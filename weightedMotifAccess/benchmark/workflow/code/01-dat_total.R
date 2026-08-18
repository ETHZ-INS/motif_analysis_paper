suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(SummarizedExperiment)
    #remotes::install_github("Jiayi-Wang-Joey/weightedMotifAccess")
    library(weightedMotifAccess)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(rtracklayer)
})

m <- wcs$mot
atacFrag <- readRDS(args$frag)
peaks <- import(args$peak)
spec <- args$genome
if (spec == "Hsapiens") {
    genome <- BSgenome.Hsapiens.UCSC.hg38
    species <- "Homo sapiens"
} else {
    genome <- BSgenome.Mmusculus.UCSC.mm10
    species <- "Mus_musculus"
}

se <- getWeightedCounts(atacFrag = atacFrag,
    ranges = peaks,
    genome = genome,
    species = species,
    resize = TRUE,
    width = 300)

saveRDS(se, args$res)
