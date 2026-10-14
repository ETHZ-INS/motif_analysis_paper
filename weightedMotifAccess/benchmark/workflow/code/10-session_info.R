x <- c(
    "dplyr",
    "tidyr",
    "edgeR",
    "limma",
    "chromVAR",
    "GenomicRanges",
    "R.utils",
    "matrixStats",
    "SummarizedExperiment",
    "affy",
    "ggplot2",
    "ggrastr",
    "patchwork",
    "data.table",
    "BSgenome.Mmusculus.UCSC.mm10",
    "BSgenome.Hsapiens.UCSC.hg38",
    "motifmatchr",
    "ggrepel",
    "reshape2",
    "MotifDb",
    "rtracklayer",
    "fields",
    "Seurat",
    "Signac",
    "ggh4x",
    "tidytext")

# install dependencies
if (!require(BiocManager))
    install.packages("BiocManager")
for (. in x)
    if (!require(basename(.), character.only = TRUE))
        BiocManager::install(., ask = FALSE, update = TRUE)

# capture session
for (. in x) {
    . <- gsub(".*/", "", .)
    suppressPackageStartupMessages(
        library(., character.only = TRUE))
}
si <- capture.output(sessionInfo())
writeLines(si, args[[1]])