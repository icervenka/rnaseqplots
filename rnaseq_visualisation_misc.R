library(readxl)
library(readr)
library(colorspace)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(ComplexHeatmap)
library(pheatmap)
library(cowplot)
library(circlize)
library(cummeRbund)
library(DESeq2)
library(limma)
library(matrixStats)
library(broom)
library(stringr)
library(tibble)
library(purrr)
library(dplyr)

gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# save_plot = function(plot, savename, width, height, unit, dpi) {
#     ggplot2::ggsave(savename,
#                     p,
#                     units = unit,
#                     dpi = dpi,
#                     width = width,
#                     height = height)
# }

# upgrade of addFeatures from cummeRbund package, was using deprecated functions
.addFeatures = function(object,features,level="genes",...){
  if(!is.data.frame(features)){
    stop("features must be a data.frame")
  }
  colnames(features)[1] = slot(object,level)@idField
  colnames(features) = make.db.names(object@DB,colnames(features),unique=T)
  DBI::dbWriteTable(object@DB,slot(object,level)@tables$featureTable,features,row.names=F,overwrite=T)
  indexQuery = paste("CREATE INDEX ",slot(object,level)@idField," ON ",
                    slot(object,level)@tables$featureTable," (",slot(object,level)@idField,")",sep="")
  res = DBI::dbExecute(object@DB,indexQuery)
}
setMethod("addFeatures",signature(object="CuffSet"),.addFeatures)

# user defined reading functions -----------------------------------------------
read_cuffdiff_diff = function(dir) {
  cuff = cummeRbund::readCufflinks(dir)
  annot = read.delim(paste0(dir,"/gene_exp.diff"),
                     sep = "\t",
                     header=T,
                     na.string="-") %>%
    dplyr::select(gene_id, gene)
  cummeRbund::addFeatures(cuff,annot,level="genes")

  # munging of differential expression data
  diff = cummeRbund::diffData(genes(cuff)) %>%
    tibble::as_tibble()
  diff = diff %>%
    dplyr::filter(value_1 > 0 & value_2 > 0) %>%
    dplyr::filter(status == "OK")
  diff = diff %>%
    dplyr::mutate(comparison = paste0(sample_1, "_", sample_2)) %>%
    dplyr::left_join(annot %>% unique())
  return(diff)
}

merge_cuffdiff_fc = function(x, y, annot, by = "gene") {
  lfcs = merge(x %>%
                 dplyr::select(gene, gene_id, log2_fold_change, q_value),
               y %>%
                 dplyr::select(gene, gene_id, log2_fold_change, q_value),
               by = by)
  lfcs = lfcs %>%
    # Fisher method of combining p-values
    dplyr::mutate(chi_pcomb = -2*(log(q_value.x)+log(q_value.y))) %>%
    # calculate significance of combined p-values
    dplyr::mutate(p_chi = pchisq(chi_pcomb, df=4, lower.tail=FALSE)) %>%
    # calculate squared error of genes for ranking
    dplyr::mutate(err_sq = (log2_fold_change.y - log2_fold_change.x)^2)

  lfcs = lfcs %>% dplyr::mutate(err_sq = (log2_fold_change.y - log2_fold_change.x)^2)
  return(lfcs)
}

read_dire = function(filename, sheet_name = "Sheet1") {

  # dire_df = function(fl) {
  #   conn = file(fl,open="r")
  #   linn = readLines(conn)
  #   linn =  matrix(linn, ncol = 4, byrow = T)
  #   print(linn)
  #   linn = data.frame(linn, stringsAsFactors = F)
  #   linn[,1] = NULL
  #   names(linn) = c("tf", "occurence", "importance")
  #   linn$occurence = as.numeric(gsub("%", "", linn$occurence))
  #   linn$importance = as.numeric(linn$importance)
  #   close(conn)
  #   return(linn)
  # }


  if(grepl(".xlsx", filename)) {
    data = readxl::read_excel(filename, sheet = sheet_name) %>%
      select(-`#`)
  } else {
    data = readr::read_delim(filename) %>%
      select(-`#`)
  }
  return(data)
}

# user defined plotting functions ----------------------------------------------
plot_pca_deseq = function(dds, group = "group") {
  vsd = DESeq2::vst(dds, blind = F)

  pcaData = DESeq2::plotPCA(vsd, intgroup=c(group), returnData=TRUE)
  percentVar = round(100 * attr(pcaData, "percentVar"))
  segments = pcaData %>%
    dplyr::group_by(!!as.symbol(group)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pcaData = merge(pcaData, segments, by = group)

  no_colors = pcaData[,group] %>% unique() %>% length()

  p = pcaData %>%
    ggplot2::ggplot(aes(PC1,
                        PC2,
                        color=!!as.symbol(group),
                        shape=!!as.symbol(group))) +
    geom_point(size=3) +
    geom_segment(aes(x = PC1,
                     y = PC2,
                     xend = xend,
                     yend = yend),
                 size = 0.5,
                 linetype = "dashed") +
    geom_point(data = segments, aes(x = xend, y = yend), size = 2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_bw() +
    scale_color_manual(values=viridis::viridis(no_colors+1)) +
    theme(legend.position = "top")

  return(p)
}

plot_pca_common = function(expression_data, metadata, color_by = "group") {
  pca = prcomp(t(expression_data %>% select(matches(metadata$sample))), scale = T)
  pcaData = as.data.frame(pca$x) %>% rownames_to_column("id")
  pcaData = merge(pcaData, metadata, by.x = "id", by.y = "sample")
  segments = pcaData %>%
    dplyr::group_by(!!as.symbol(color_by)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pcaData = merge(pcaData, segments, by = color_by)

  percentage = round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage = paste0(colnames(pca), "(", paste(as.character(percentage), "%", ")", sep=""))

  no_colors = pcaData[,color_by] %>% unique() %>% length()

  p = pcaData %>%
    ggplot2::ggplot(aes(PC1,
                        PC2,
                        color=!!as.symbol(color_by),
                        shape=!!as.symbol(color_by))) +
    geom_point(size=3) +
    geom_segment(aes(x = PC1,
                     y = PC2,
                     xend = xend,
                     yend = yend),
                 size = 0.5,
                 linetype = "dashed") +
    geom_point(data = segments, aes(x = xend, y = yend), size = 2) +
    xlab(paste0("PC1: ",percentage[1]," variance")) +
    ylab(paste0("PC2: ",percentage[2]," variance")) +
    theme_bw() +
    scale_color_manual(values=viridis::viridis(no_colors+1)[-length(viridis::viridis(no_colors+1))]) +
    theme(legend.position = "top")

  return(p)
}

plot_corr = function(gene_expr, param_values) {
 p = cbind.data.frame(expr = gene_expr,
                       param = param_values,
                       color_by = param_values) %>%
    ggplot2::ggplot(aes(x = expr, y = param, color = color_by)) +
    geom_point()
 return(p)
}

plot_lfc_scatter = function(lfc_data) {
  p = lfc_data %>%
    dplyr::filter(p_chi < 0.1) %>%
    ggplot2::ggplot(aes(x = log2_fold_change.x,
                        y = log2_fold_change.y,
                        text = gene)) +
    geom_point(aes(size = -log10(p_chi), color = err_sq)) +
    colorspace::scale_colour_continuous_sequential("viridis", rev = F) +
    theme_bw() +
    labs(size="-log10(p-value)", colour="squared\nerror")
  return(p)
}

plot_venn2 = function(x, y, names, savename) {
  VennDiagram::venn.diagram(
    x = list(tolower(x), tolower(y)),
    category.names = names,
    filename = savename,
    output=TRUE,

    # Output features
    imagetype="tiff" ,
    units = 'mm',
    height = 35 ,
    width = 35 ,
    resolution = 600,
    compression = "lzw",

    # Circles
    lwd = 2,
    lty = 'blank',
    fill = viridis(2),

    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-25, 25),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
  )
}

plot_dire = function(df, color = "steelblue") {
  p = df %>%
    ggplot2::ggplot(aes(x = Occurrence, y = Importance)) +
    geom_point(color = color) +
    theme_bw()
  return(p)
}

plot_dire_labeled = function(df,
                             occurrence_threshold = 0.05,
                             importance_threshold = 0.05,
                             dot_size = 3.5,
                             color = "steelblue") {

  p = df %>%
    plot_dire(color) +
    geom_label_repel(data = . %>%
                       filter(Occurrence > occurrence_threshold | Importance > importance_threshold),
                     aes(x = Occurrence, y = Importance, label = `Transcription Factor`),
                     size = dot_size)
  return(p)
}

# heatmap all
plot_heatmap = function(expression_df, metadata,
                             gene_list = NULL,
                             geneid_colname = SYMBOL,
                             metadata_sample_colname = sample,
                             gene_ranking_fun = matrixStats::rowVars,
                             cell_dims = c(10, 1),
                             palette = "RdBu",
                             ...) {

  sc = deparse(substitute(metadata_sample_colname))

  if(!is.null(gene_list)) {
    if(is.character(gene_list)) {
      expression_df = expression_df %>%
        dplyr::filter(!!as.symbol(geneid_colname) %in% gene_list)
    } else if(is.numeric(gene_list)) {
      expr_matrix = expression_df %>%
        dplyr::select(dplyr::matches(metadata[[sc]])) %>%
        as.matrix()
      expression_df = expression_df[order(gene_ranking_fun(expr_matrix), decreasing = T),][1:gene_list,]
    } else {
      stop("Unrecognized type of argument for gene list.")
    }

  }

  expression_df = expression_df %>%
    dplyr::select(dplyr::matches(metadata[[sc]]))
  
  p = expression_df %>%
  pheatmap::pheatmap(scale = "row",
                     color = colorRampPalette(rev(RColorBrewer::brewer.pal(
                       n = 7, name = palette)))(100),
                     cellwidth = cell_dims[1],
                     cellheight = cell_dims[2],
                     border_color = "white",
                     treeheight_row = 15,
                     treeheight_col = 20,
                     annotation_col = metadata %>%
                       tibble::column_to_rownames(sc),
                     ...)
  return(p)
}

# # ma plot
# ko_vs_wt_all %>%
#   mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
#   mutate(significant = as.factor(case_when(padj < 0.05 & log2FoldChange > 0 ~ 1,
#                                            padj < 0.05 & log2FoldChange < 0 ~ 2,
#                                            TRUE ~ 0))) %>%
#   ggplot(aes(x = baseMean, y = log2FoldChange, color = significant)) +
#   geom_point() +
#   scale_color_discrete(type = c("gray60", "darkred", "steelblue")) +
#   theme_bw() +
#   theme(legend.position = "none")


# volcano plot all
volcano_plot = function(results,
                        label = SYMBOL,
                        x = log2FoldChange,
                        y = pvalue,
                        sig_threshold = 0.05,
                        log2fc_threshold = 0.585, # rougly 1.5x
                        filter_sig_on = padj,
                        label_bottom_n = 5,
                        label_top_n = 5,
                        label_genes = NULL, # if specified overrides top_n/bottom_n arguments,
                        add_vhlines = TRUE,
                        vhline_color = "gray40",
                        vhline_type = "dashed",
                        xlab_label = expression(paste("lo", g[2], "FC")),
                        ylab_label = expression(paste("lo", g[10], "FDR")),
                        color_palette = c("gray70", "steelblue", "darkred")) {

  require(ggrepel)

  # validate function arguments
  if(length(color_palette) != 3) {
    stop("Color palette requires a vector length 3.")
  }

  results_fil = results %>%
    dplyr::select({{ label }}, {{ x }}, {{ y }}, {{ filter_sig_on }}) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(significant = dplyr::case_when({{ filter_sig_on }} <= sig_threshold & {{ x }} > log2fc_threshold ~ "up",
                                                 {{ filter_sig_on }} <= sig_threshold & {{ x }} <= (-1)*log2fc_threshold ~ "down",
                                                 TRUE ~ "unchanged")) %>%
    dplyr::mutate(significant = factor(significant, levels = c("unchanged", "down", "up")))

  p = results_fil %>%
    ggplot2::ggplot(aes(x = {{ x }}, y = -log10({{ y }}))) +
    geom_point(aes(color = significant), size = 1)

  if(is.null(label_genes)) {
    p = p +
      geom_label_repel(data = results_fil %>%
                         filter(significant != "unchanged") %>%
                         top_n(-label_bottom_n, {{ x }}),
                       aes(label = {{ label }}),
                       min.segment.length = unit(0, 'lines')) +
      geom_label_repel(data = results_fil %>%
                         filter(significant != "unchanged") %>%
                         top_n(label_top_n, {{ x }}) ,
                       aes(label = {{ label }}),
                       min.segment.length = unit(0, 'lines'))
  } else if(is.vector(label_genes, mode = "character")) {
    p = p +
      geom_label(data = results_fil %>% dplyr::filter({{ label }} %in% label_genes),
                 aes(label = {{ label }}),
                 min.segment.length = unit(0, 'lines'))
  }

  if(add_vhlines) {
    p = p +
      geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold),
                 linetype = vhline_type,
                 color = vhline_color) +
      geom_hline(yintercept = -log10(sig_threshold),
                 linetype = vhline_type,
                 color=vhline_color)
  }

  p = p +
    theme_bw() +
    xlab(xlab_label) +
    ylab(ylab_label) +
    theme(legend.position = "none") +
    scale_color_manual(values = color_palette)
  p
}

volcano_plot_cuffdiff = function(results, sample1, sample2, ...) {
  results = results %>%
    dplyr::filter(sample_1 == sample1 & sample_2 == sample2)
  volcano_plot(results, ...)
}


plot_heatmap_fc = function(expression_data, diffexp_data, metadata, gene_list,
                           id_colname = "SYMBOL", fc_colname = "log2FoldChange",
                           padj_colname = "padj") {

  expression_data_fil = expression_data %>%
    dplyr::arrange(!!as.symbol(id_colname)) %>%
    dplyr::filter(!!as.symbol(id_colname) %in% gene_list) %>%
    dplyr::select(!!as.symbol(id_colname), matches(metadata$sample)) %>%
    tibble::column_to_rownames(id_colname)

  diffexp_data_fil = diffexp_data %>%
    dplyr::arrange(!!as.symbol(id_colname)) %>%
    dplyr::filter(!!as.symbol(id_colname) %in% gene_list)

  log2fc_vals = diffexp_data_fil %>%
    dplyr::pull(!!as.symbol(fc_colname))
  log2fc_colors = ifelse(log2fc_vals < 0, "steelblue", "darkred")

  pvals = diffexp_data_fil %>%
    dplyr::mutate(log_pval = -log10(!!as.symbol(padj_colname))) %>%
    dplyr::pull(log_pval)
  pvals_colors = circlize::colorRamp2(c(-log10(0.05), max(pvals)),
                                      c("white", "steelblue"))
  pvalues_legend = ComplexHeatmap::Legend(col_fun = pvals_colors,
                                          title = "-log10(p-value)")

  har = ComplexHeatmap::rowAnnotation(
    pvalue = ComplexHeatmap::anno_simple(pvals,
                                         col = pvals_colors,
                                         gp = gpar(col = "black", lwd = 1)),
    log2fc = ComplexHeatmap::anno_barplot(log2fc_vals,
                                          baseline = 0,
                                          bar_width = 0.9,
                                          gp = gpar(fill = log2fc_colors, col = "white")),
    simple_anno_size = unit(0.5, "cm"), width = unit(2, "cm"),
    gap = unit(2, "mm"))

  col_colors = gg_color_hue(length(unique(metadata$group)))
  names(col_colors) = unique(metadata$group)

  hat = ComplexHeatmap::columnAnnotation(
    genotype = metadata$group,
    col = list(genotype = col_colors),
    border = TRUE)

  ht = ComplexHeatmap::Heatmap(
    expression_data_fil %>%
      as.matrix %>%
      t %>%
      scale %>%
      t,
    name = "expression",
    right_annotation = har,
    top_annotation = hat,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    border_gp = gpar(col = "black", lty = 1),
    rect_gp = gpar(col = "black", lwd = 1),
    width = ncol(expression_data_fil)*unit(5, "mm"),
    height = nrow(expression_data_fil)*unit(5, "mm"))

  ComplexHeatmap::draw(
    ht,
    annotation_legend_list = list(pvalues_legend),
    merge_legends = TRUE)
}

collate_pathways = function(pathway_files_basepath, pattern, padj_threshold = 0.05) {
  if(stringr::str_sub(path_to_pathway_files, -1) != "/") {
    pathway_files_basepath = paste0(pathway_files_basepath, "/")
  }

  pathway_files = list.files(pathway_files_basepath,
                             pattern = pattern,
                             full.names = F)

  cp_unified_colnames = c("ID",	"Description",	"GeneRatio/NES",	"pvalue",
                          "p.adjust",	"SYMBOL",	"ENTREZID",	"log2FoldChange")

  diffexp_pathways = purrr::map_dfr(pathway_files, function(x) {
    readr::read_delim(paste0(pathway_files_basepath, x),
                      delim = "\t", escape_double = FALSE,
                      trim_ws = TRUE,
                      show_col_types = F)  %>%
      setNames(cp_unified_colnames) %>%
      # TODO needs to be changed
      dplyr::mutate(source = gsub("_contrast.txt", "", x),
                    ID = as.character(ID)) %>%
      # might not reflect actual columns
      dplyr::select(-SYMBOL, -ENTREZID, -log2FoldChange) %>%
      unique()
  }) %>%
    filter(p.adjust < padj_threshold)
  return(diffexp_pathways)
}


plot_pathways_meta = function(df, top_pathways = 30) {
  pathways_summary = df %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(n_pathways = n(),
                     mean_enrichment = mean(`GeneRatio/NES`))

  p1 = pathways_summary %>%
    ggplot2::ggplot() +
    geom_histogram(aes(x = n_pathways),
                   bins = 100,
                   fill = "steelblue") +
    theme_bw()

  # top pathway contributors
  p2 = pathways_summary %>%
    dplyr::top_n(30, n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(aes(x = n_pathways, y = source)) +
    geom_bar(stat="identity", fill = "steelblue") +
    theme_bw()

  # bottom pathway contributors
  p3 = pathways_summary %>%
    dplyr::top_n(30, -n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(aes(x = n_pathways, y = source)) +
    geom_bar(stat="identity", fill = "steelblue") +
    theme_bw()

  return(plot_grid(p1, p2, p3, nrow = 1))
}

plot_pathway_bargraph = function(df, pathway_source, top_n = 20, truncate_desc = 80) {
  df %>%
    dplyr::filter(grepl(pathway_source, source)) %>%
    dplyr::arrange(-abs(`GeneRatio/NES`)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::arrange(`GeneRatio/NES`) %>%
    dplyr::mutate(Description = stringr::str_trunc(Description, truncate_desc)) %>%
    dplyr::mutate(Description = factor(Description, levels = (.) %>% pull(Description))) %>%
    ggplot2::ggplot(aes(x = `GeneRatio/NES`,
               y = Description,
               color = -log10(p.adjust),
               fill = -log10(p.adjust))) +
    geom_bar(stat="identity") +
    ylab("") +
    theme_bw()
}
