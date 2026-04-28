suppressPackageStartupMessages({
    library(chromVAR)
    library(limma)
    library(SummarizedExperiment)
    library(BiocParallel)
})
register(MulticoreParam(4))
Sys.setenv(OMP_NUM_THREADS = "4")
run_once <- function(x, genome, motif, seed) {
    set.seed(seed)

    x <- filterPeaks(x, non_overlapping = TRUE)
    x <- addGCBias(x, genome = genome)

    bg <- getBackgroundPeaks(object = x, niterations = 100)
    motif_ix <- matchMotifs(motif, x, genome = genome)

    dev <- chromVAR::computeDeviations(
        object = x,
        annotations = motif_ix,
        expectation = computeExpectations(x),
        background_peaks = bg
    )

    group_id <- ifelse(grepl("CTRL", colnames(x), ignore.case = TRUE), "A", "B")
    group_id <- relevel(factor(group_id), ref = "A")
    design <- model.matrix(~ group_id)

    fit <- eBayes(lmFit(assay(dev, "z"), design))
    res <- topTable(fit, n = Inf, sort.by = "none")
    res$name <- rownames(res)

    list(
        res = res,
        z = assay(dev, "z")
    )
}

fun <- function(x, genome, motif, seeds = c(42, 2026)) {
    runs <- lapply(seeds, function(s) run_once(x, genome, motif, seed = s))

    res_list <- lapply(runs, `[[`, "res")
    z_list <- lapply(runs, `[[`, "z")

    common_names <- Reduce(intersect, lapply(res_list, rownames))
    res_list <- lapply(res_list, function(r) r[common_names, , drop = FALSE])

    avg_res <- res_list[[1]]
    num_cols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
    num_cols <- intersect(num_cols, colnames(avg_res))

    for (cc in num_cols) {
        avg_res[[cc]] <- Reduce(`+`, lapply(res_list, `[[`, cc)) / length(res_list)
    }

    avg_res <- avg_res[order(avg_res$P.Value, -abs(avg_res$logFC)), , drop = FALSE]
    avg_res$name <- rownames(avg_res)
    avg_res$rank <- seq_len(nrow(avg_res))

    common_z <- Reduce(intersect, lapply(z_list, rownames))
    z_list <- lapply(z_list, function(z) z[common_z, , drop = FALSE])
    avg_z <- Reduce(`+`, z_list) / length(z_list)

    list(res = avg_res, z = avg_z, runs = runs)
}

