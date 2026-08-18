suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(chromVAR)
    library(limma)
    library(edgeR)
    library(pryr)
})

fun <- function(se, genome=NULL, motif=NULL) {
    x <- assay(se, "z")
    group_id <- ifelse(grepl("CTRL", colnames(x), ignore.case = TRUE), "A", "B")
    group_id <- relevel(factor(group_id), ref = "A")
    design <- model.matrix(~ group_id)

    fit <- lmFit(x, design)
    fit <- eBayes(fit, trend = sqrt(rowData(se)$N))
    res <- topTable(fit, n = Inf, sort.by = "none")
    res$name <- rownames(res)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    res
}
