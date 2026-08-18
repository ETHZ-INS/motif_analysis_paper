suppressPackageStartupMessages({
    library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    source("workflow/code/utils.R")
    library(MotifDb)
    library(BiocParallel)
})
register(MulticoreParam(4))
Sys.setenv(OMP_NUM_THREADS = "4")

m <- toupper(wcs$mot)
se <- readRDS(args$dat)

spec <- args$genome

if (spec == "Hsapiens") {
    genome <- BSgenome.Hsapiens.UCSC.hg38
    species <- "Homo sapiens"
    motif <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
} else {
    genome <- BSgenome.Mmusculus.UCSC.mm10
    species <- "Mus_musculus"
    #motif <- readRDS("data/Mmotifs.rds")
    motif <- getNonRedundantMotifs(format="PFMatrix", species="Mmusculus")
    banp <- readRDS("data/BANP.PFMatrix.rds")
    motif$BANP <- banp
    Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
    motif$NR1H3 <- Hmotifs$NR1H3
}
source(args$fun)

res <- fun(se,
           genome = genome,
           motif = motif)

if (!is.null(res$z)) {
    #saveRDS(res$z, paste0("outs/dat/chromVAR-z-", wcs$mot, ".rds"))
    res <- res$res
}


if (m=="MYC") {
    target <- paste("MYC", "MAX", sep=",")
} else if (m=="NR1H4"){
    target <- paste("NR1H4","RXRA","RXRB",sep=",")
} else if (m=="ESR1"){
    target <- paste("ESR1", "ESR2", sep=",")
} else if (m=="PPARA") {
    target <- paste("PPARA", "RXRA", "RXRB", "FXR", sep=",")
} else {
    target <- m
}


if (grepl("weight", args$dat)) {
    df <- data.frame(res, truth=target,
        dif=wcs$dif, method="weight", row.names=NULL)

} else if (grepl("insert", args$dat)) {
    df <- data.frame(res, truth=target,
        dif=wcs$dif, method="inserts", row.names=NULL)
} else {
    df <- data.frame(res, truth=target,
        dif=wcs$dif, method="origin", row.names=NULL)
}

saveRDS(df, args$res)
