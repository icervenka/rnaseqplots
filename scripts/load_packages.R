# package list used to load packages and fix namespace qualification
# using origin r package
pack <- c(
    "readxl",
    "readr",
    "colorspace",
    "RColorBrewer",
    "viridis",
    "ggplot2",
    "ggrepel",
    "VennDiagram",
    "futile.logger",
    "ComplexHeatmap",
    "pheatmap",
    "cowplot",
    "circlize",
    "egg",
#    "SummarizedExperiment",
    "cummeRbund",
    "DESeq2",
#    "limma",
    "matrixStats",
#    "broom",
    "stringr",
    "tibble",
    "purrr",
    "dplyr"
)

lapply(pack, function(x) {
    suppressPackageStartupMessages(library(x, character.only = TRUE))
})

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")