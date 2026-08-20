# Code accompanying the motif accessibility analysis paper by Wang, Sonder et al.

For the actual packages, see [betterChromVAR](https://bioconductor.org/packages/devel/bioc/html/betterChromVAR.html), [epiwraps](https://github.com/ETHZ-INS/epiwraps), and [weightedMotifAccess](https://github.com/Jiayi-Wang-Joey/weightedMotifAccess).

The code to reproduce the following figures:

- Figure 1 (technical variations in ATAC-seq) in [normalization/bias_examples.Rmd](normalization/bias_examples.Rmd)
- Figure 2 (TF activity inference benchmark) in [weightedMotifAccess/](weightedMotifAccess/)
- Figure 3 (comparison of betterChromVAR to chromVAR, CVnorm, etc.): figure in [figures/betterChromVAR.Rmd](figures/betterChromVAR.Rmd), snakemake for the actual benchmark in [chromVAR_benchmark_30kcells/](chromVAR_benchmark_30kcells/).
- Figure 4 (motif interactions analysis) in [motifInteractions/motifInter.Rmd](motifInteractions/motifInter.Rmd)
- Supplementary Figure 1 (ATAC duplicates) in [ATAC_duplicates/](ATAC_duplicates/)
- Supplementary Figure 4 (benchmark knn variants at single-cell level) in [chromVAR_benchmark_30kcells/bin_vs_knn.Rmd](chromVAR_benchmark_30kcells/bin_vs_knn.Rmd).
