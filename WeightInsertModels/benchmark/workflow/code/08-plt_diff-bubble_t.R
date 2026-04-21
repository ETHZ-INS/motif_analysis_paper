# setwd("~/motif_analysis_paper/WeightInsertModels/benchmark")
# args <- list(list.files("outs/dat", "^dif-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(viridis)
    library(patchwork)
    library(dplyr)
    library(tibble)
})

res <- lapply(args[[1]], readRDS)
d <- lapply(res, select, rank, t, method,
    truth, name, dif, adj.P.Val) |>
    do.call(what=rbind) |>
    mutate(method=paste(method,dif,sep=".")) |>
    group_by(method, truth) |>
    filter(name %in% unlist(strsplit(truth, ","))) |>
    mutate(min_rank = min(rank)) |>
    filter(rank == min_rank) |>
    ungroup() |>
    select(-min_rank) |>
    rename("label_rank"="rank") |>
    mutate(significance=ifelse(adj.P.Val < 0.05, "TRUE", "FALSE"),
        rank=case_when(label_rank > 100 ~ 100, TRUE ~ label_rank))

d$col1 = ifelse(d$rank>=100,"black","white")
d$col2 = ifelse(abs(d$t)<15,"black","white")

aes_r <- list(geom_point(aes(
    color=rank, size=significance)),
    scale_size_manual(values = c("TRUE"=9, "FALSE"=6.5, 1)),
    scale_color_viridis(discrete=FALSE, option = "D",
        breaks=c(5, 50, 100),
        labels=c("5", "50", ">= 100")),
    theme_bw(),
    labs(x="perturbed TFs", y="method"),
    theme(legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        legend.position="bottom"),
    guides())


aes_t <- list(geom_point(aes(color=abs(t), size=significance)),
    scale_size_manual(values = c("TRUE"=9, "FALSE"=6.5, 1)),
    scale_color_viridis(discrete=FALSE,
        option = "magma",
        direction=-1,
        guide = guide_legend(nrow = 2)
        ),
    theme_bw(),
    labs(x="perturbed TFs", y="method"),
    theme(legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        legend.position="bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
    guides())


d$method <- factor(d$method, c("origin.limma-voom", "origin.chromVAR.limma", "origin.betterChromVAR", "weight.limma-voom", "origin.betterChromVAR.FL"),
    c("raw crossprod", "chromVAR", "betterChromVAR", "weights", "betterChromVAR.FL"))
o <- d
signedsqrt <- scales::new_transform("signedsqrt", \(x) sign(x)*sqrt(abs(x)), inverse = \(x) sign(x)*x^2)
o$t2 <- o$t * (2*(o$truth %in% c("NR1H3", "NR1H4"))-1)
o$dataset <- gsub(",","/\n",o$truth)
o$truth <- gsub(",","/\n",o$truth)

cols <- c(
    "raw crossprod"      = "#BDBDBD",
    "chromVAR"           = "#9ECAE1",
    "betterChromVAR"     = "#4292C6",
    "betterChromVAR.FL"  = "#08519C",
    "weights"            = "#E7298A"
)

labels <- c(
    "raw crossprod"     = "Raw crossprod",
    "chromVAR"          = "chromVAR (GC + enrichment)",
    "betterChromVAR"    = "betterChromVAR (GC + enrichment)",
    "betterChromVAR.FL" = "betterChromVAR.FL (GC + FL + enrichment)",
    "weights"           = "weights (GC + FL + enrichment)"
)

gt <- ggplot(o, aes(reorder(dataset,t2), t2, fill=method)) +
    geom_col(position=position_dodge2()) +
    scale_y_continuous(trans=signedsqrt) +
    theme_bw() +
    #scale_fill_brewer(palette = "Paired") +
    scale_fill_manual(values = cols, labels = labels,
        guide = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom",
        legend.title = element_text(face = "bold")
        ) +
    labs(x="Dataset", y=NULL) +
    ggtitle("t-value (higher absolute value = better)")




gg <- ggplot(o, aes(reorder(truth,-sqrt(rank)),
    reorder(method,-sqrt(rank)))) + aes_r +
    geom_text(aes(label=label_rank), size=2.5, color=o$col1)+
    ggtitle("rank (lower = better)") + labs(y = NULL) +
    theme(legend.title = element_text(face = "bold")) +
    gt +
    plot_layout(nrow=1) &
    plot_annotation(tag_levels = list(c("B", "C"))) &
    theme(
        plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0, 1)
    )

ggsave(args[[2]], gg, width=35, height=12, units="cm")
