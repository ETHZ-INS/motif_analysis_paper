suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(SummarizedExperiment)
    library(WeightInsertModels)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(rtracklayer)
    library(betterChromVAR)
    library(Rsamtools)
})

m <- wcs$mot
fs <- list.files(paste0("/mnt/plger/plger/DTFAB/fullFrags/",
    m, "/seq_files"), pattern = "\\.bam$|\\.bed$",
    full.names = TRUE)

# bgzip and tabix-index any plain BED files so tabixChrApply can read them
bed_files <- fs[grepl("\\.bed$", fs)]
if(length(bed_files) > 0){
    idx_dir <- file.path("/mnt/plger/jwang/data/bed_indexed", m)
    dir.create(idx_dir, recursive=TRUE, showWarnings=FALSE)
    fs <- c(
        fs[!grepl("\\.bed$", fs)],
        sapply(bed_files, function(f){
            out <- file.path(idx_dir, paste0(basename(f), ".gz"))
            if(!file.exists(out)) bgzip(f, dest=out, overwrite=FALSE)
            if(!file.exists(paste0(out, ".tbi"))) indexTabix(out, format="bed")
            out
        })
    )
}

lst <- as.list(fs)
names(lst) <- basename(fs)

spec <- args$genome
species_name <- ifelse(spec == "Hsapiens", "Homo_sapiens", "Mus_musculus")

peakRanges <- import(args$peak)
peakRanges <- keepStandardChromosomes(peakRanges, species=species_name, pruning.mode="coarse")
motifRanges <- readRDS(args$motfrange)
motifRanges <- keepStandardChromosomes(motifRanges, species=species_name, pruning.mode="coarse")
# match motifRanges style to peakRanges (which matches the fragment files)
seqlevelsStyle(motifRanges) <- seqlevelsStyle(peakRanges)[1]

if (spec == "Hsapiens") {
    genome <- BSgenome.Hsapiens.UCSC.hg38
    species <- "Homo sapiens"
} else {
    genome <- BSgenome.Mmusculus.UCSC.mm10
    species <- "Mus_musculus"
}
motifRanges <- subsetByOverlaps(motifRanges, peakRanges)

shift <- fifelse(m %in% c("GATA1", "GATA2", "MYC", "KLF1", "RUNX1", "RUNX2"),
    FALSE, TRUE)

res <- WeightInsertModels::weightedInsertions(lst, peakRanges,
    motifRanges, ncores=1, shift=shift)
pc <- res$peakCounts
seqlevelsStyle(pc) <- "UCSC"

w <- which(rowSums(assay(pc)) > 0)
pc <- addGCBias(pc[w,], genome=genome)
dev <- betterChromVAR::computeDeviationsWeighted(res$wInsCounts,
    pc, annotations=res$peakAnnotation[w,])


saveRDS(dev, args$res)
