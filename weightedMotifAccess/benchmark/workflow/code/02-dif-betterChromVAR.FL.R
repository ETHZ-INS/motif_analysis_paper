suppressPackageStartupMessages({
    library(betterChromVAR)
    library(limma)
    library(SummarizedExperiment)
    library(BiocParallel)
})

fun <- \(x, genome, motif) {
    x <- filterPeaks(x, non_overlapping = TRUE)
    x <- addGCBias(x,
        genome = genome)
    rowData(x)$flbias <- log10(rowMeans(assay(x, "median_width")))
    keep <- is.finite(rowData(x)$bias)
    x <- x[keep, ]
    motif_ix <- matchMotifs(motif,
        x,
        genome = genome)
    dev <- betterChromVAR(x, motif_ix)
    group_id <- ifelse(grepl("CTRL", colnames(x), ignore.case = TRUE), "A", "B")
    group_id <- relevel(factor(group_id), ref = "A")
    design <- model.matrix(~ group_id)
    fit <- eBayes(lmFit(assay(dev, "z"), design))
    res <- topTable(fit, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$name <- rownames(res)
    res$rank <- seq_len(nrow(res))
    list(res=res, z=assay(dev, "z"))
}
