#' A tibble with differentially regulated pathways identified by clusterProfiler
#'
#' @format tibble with 116 rows and 5 columns
#' \describe{
#'   \item{ID}{Pathway ID, chr}
#'   \item{Description}{Pathway name, chr}
#'   \item{GeneRatio/NES}{GeneRatio for ORA or NES for GSEA, num}
#'   \item{pvalue}{P-value for being differentially regulated, num}
#'   \item{p.adjust}{Adjusted p-value for being differentially regulated, num}
#' }
"cp_pathways"

#' A list of tibbles with differentially expressed genes identified by DESeq2.
#' Each tibble corresponds to one comparison between groups.
#'
#' @format a list of tibbles of length 2 each with 19229 rows and 10 columns
#' \describe{
#'   \item{ex_1}{group comparison 1}
#'     \item{ENSEMBL}{ensembl gene id, chr}
#'     \item{ENTREZID}{entrez id, num}
#'     \item{SYMBOL}{gene symbol, chr}
#'     \item{GENENAME}{long gene name, chr}
#'     \item{baseMean}{mean of normalized counts for all samples, num}
#'     \item{log2FoldChange}{log2 fold-change difference of expression, num}
#'     \item{lfcSE}{standard error, num}
#'     \item{stat}{Wald statistic, num}
#'     \item{pvalue}{P-value for differential expression, num}
#'     \item{padj}{Adjusted p-value for differential expression, num}
#'   \item{ex_2}{group comparison 2}
#'     \item{ENSEMBL}{ensembl gene id, chr}
#'     \item{ENTREZID}{entrez id, num}
#'     \item{SYMBOL}{gene symbol, chr}
#'     \item{GENENAME}{long gene name, chr}
#'     \item{baseMean}{mean of normalized counts for all samples, num}
#'     \item{log2FoldChange}{log2 fold-change difference of expression, num}
#'     \item{lfcSE}{standard error, num}
#'     \item{stat}{Wald statistic, num}
#'     \item{pvalue}{P-value for differential expression, num}
#'     \item{padj}{Adjusted p-value for differential expression, num}
#' }
"diffexp_data"

#' A tibble with upstream transcription factors (TFs) identified by dire.dcode.org
#'
#' @format tibble with 241 rows and 4 columns
#' \describe{
#'   \item{Transcription Factor}{transcription factor name, chr}
#'   \item{Occurrence}{TF occurrence, num}
#'   \item{Importance}{TF importance, num}
#'   \item{file}{name of the file the data was parsed from, chr}
#' }
"dire"

#' A tibble with normalized gene counts for all samples obtained from DESeq2
#'
#' @format tibble with 19229 rows and 19 columns
#' \describe{
#'   \item{ENSEMBL}{ensembl gene id, chr}
#'   \item{ENTREZID}{entrez id, num}
#'   \item{SYMBOL}{gene symbol, chr}
#'   \item{GENENAME}{long gene name, chr}
#'   \item{cgrp1}{normalized counts for sample cgrp1, num}
#'   \item{cgrp2}{normalized counts for sample cgrp2, num}
#'   \item{cgrp3}{normalized counts for sample cgrp3, num}
#'   \item{ctrl1}{normalized counts for sample ctrl1, num}
#'   \item{ctrl2}{normalized counts for sample ctrl2, num}
#'   \item{ctrl3}{normalized counts for sample ctrl3, num}
#'   \item{ne1}{normalized counts for sample ne1, num}
#'   \item{ne2}{normalized counts for sample ne2, num}
#'   \item{ne3}{normalized counts for sample ne3, num}
#'   \item{npy1}{normalized counts for sample npy1, num}
#'   \item{npy2}{normalized counts for sample npy2, num}
#'   \item{npy3}{normalized counts for sample npy3, num}
#'   \item{sp1}{normalized counts for sample sp1, num}
#'   \item{sp2}{normalized counts for sample sp2, num}
#'   \item{sp3}{normalized counts for sample sp3, num}
#' }
"expression_data"

#' A named list of character vectors containing user-defined gene lists parsed
#' from json file
#'
#' @format named list with 6 items
#' \describe{
#'   \item{mito_1}{mitchondrial genes, chr vector}
#'   \item{metabolism}{metabolism genes, chr vector}
#'   \item{volcano_example}{gene lists for example volcano plot, chr vector}
#'   \item{heatmap_fc_example_1}{gene lists for example heatmap, chr vector}
#'   \item{heatmap_fc_example_2}gene lists for example heatmap, chr vector}
#'   \item{oxhphos}{oxphos genes, chr vector}
#' }
"gene_lists"

#' A tibble with differentialy regulated pathways identified by clusterProfiler
#'
#' @format tibble with 116 rows and 5 columns
#' \describe{
#'   \item{ID}{Pathway ID, chr}
#'   \item{Description}{Pathway name, chr}
#'   \item{GeneRatio/NES}{GeneRatio for ORA or NES for GSEA, num}
#'   \item{pvalue}{P-value for being differentially regulated, num}
#'   \item{p.adjust}{Adjusted p-value for being differentially regulated, num}
#' }
"gsea_pathways"

#' A tibble with sample-group mapping
#'
#' @format tibble with 15 rows and 2 columns
#' \describe{
#'   \item{sample}{sample name/id, chr}
#'   \item{group}{group name, chr}
#' }
"metadata"

#' A named list used for replacing gene names. Names contain gene names to be
#' replaced and values contain new gene names
#'
#' @format a named list of length 13
#' }
"mito_atp_replacer"

#' A named list of either numeric or character vectors containing pathway IDs or
#' names to filter from a tibble containing differentially regulated pathways
#'
#' @format a named list of length 2
#' \describe{
#'   \item{example_1}{example pathways to filter, num vector}
#'   \item{example_volcano}{example pathways to filter for volcano plot, num vector}
#' }
"pathway_lists"
