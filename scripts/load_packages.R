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
    "ComplexHeatmap",
    "pheatmap",
    "cowplot",
    "circlize",
    "egg",
    "SummarizedExperiment",
    "cummeRbund",
    "DESeq2",
    "limma",
    "matrixStats",
    "broom",
    "stringr",
    "tibble",
    "purrr",
    "dplyr"
)

# TODO add suppress startup messages
lapply(pack, library, character.only = TRUE)