library(magrittr, include.only = "%>%")

#' Save ggplot graph with fixed panel sizes
#'
#' Useful when the axis labels between group of graphs are of different sizes
#' and you want to keep the plotting areas of equal dimensions
#'
#' @param file character string, file name of the exported file
#' @param plot plot to save, defaults to last plot
#' @param units character string, plot size units, one of
#' ("in", "cm", "mm", or "px"), default: "mm"
#' @param margin double, plot panel size in specified units, default: 1
#' @param panel_width double, plot panel width in specified units, default: 4
#' @param panel_height double, plot panel height in specified units, default: 4
#' @param width double, plot width in specified units
#' @param height double, plot height in specified units
#'
#' @return ggplot graph with fixed panel sizes
#' @export
#'
#' @importFrom egg set_panel_size
#' @importFrom grid unit
#' @importFrom ggplot2 ggsave
#'
#' @examples
ggsave_fixed <- function(file, plot = last_plot(), units = "mm", margin = 1,
                         panel_width = 4, panel_height = 4,
                         width = round(dev.size()[1], digits = 1),
                         height = round(dev.size()[1], digits = 1)) {
  pf <- egg::set_panel_size(
    p = plot,
    file = NULL,
    margin = grid::unit(margin, units),
    width = grid::unit(panel_width, units),
    height = grid::unit(panel_height, units)
  )
  ggplot2::ggsave(file,
    plot = pf,
    units = units,
    width = width,
    height = height
  )
}

#' Create normalized gene expression file for GSEA analysis
#'
#' @param data data frame containing gene expression data with gene ID/symbol
#' column and sample names as columns, only samples found in metadata will be
#' used
#' @param metadata data frame containing sample and group information, requires
#' a column for samples and column for grouping
#' @param id_colname column name with gene IDs/symbols in the data file,
#' supplied as variable
#' @param data_description_col optional column name with gene descriptions
#' in the data file, supplied as variable. If NULL the resulting DESCRIPTION
#' column in the GSEA expression file will be empty. default: NULL
#' @param sample_colname columns column name of the sample column in the
#' metadata file, supplied as variable. default: sample
#' @param out character string, optional output filename. If NULL, only
#' formatted data frame is returned. default: NULL
#'
#' @return formatted data frame with normalized gene expression for  GSEA
#' analysis, if out is supplied, it will besaved to the corresponding file
#' @export
#'
#' @importFrom dplyr select mutate relocate
#'
#' @examples
create_gsea_normalized <- function(data,
                                   metadata,
                                   id_colname,
                                   data_description_col = NULL,
                                   sample_colname = sample,
                                   out = NULL) {
  # create string from supplied variables, needed for certain functions
  samples <- metadata[[deparse(substitute((sample_colname)))]]
  data_cols <- names(data)[names(data) %in% samples]

  # create the data frame to output
  # see https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
  out_df <- data %>%
    dplyr::select(NAME = {{ id_colname }}, dplyr::any_of(data_cols))

  # add description column based on user settings
  if (is.null(data_description_col)) {
    out_df <- out_df %>%
      dplyr::mutate(DESCRIPTION = "") %>%
      dplyr::relocate(DESCRIPTION, .after = NAME)
  } else {
    out_df <- out_df %>%
      dplyr::mutate(DESCRIPTION = {{ data_description_col }}) %>%
      dplyr::relocate(DESCRIPTION, .after = NAME)
  }

  # write to file if specified
  if (!is.null(out)) {
    write.table(out_df,
      out,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  return(out_df)
}

#' Create class file for GSEA analysis
#'
#' @param data data frame containing gene expression data with gene ID/symbol
#' column and sample names as columns, only samples found in metadata will be
#' used
#' @param metadata data frame containing sample and group information, requires
#' a column for samples and column for grouping
#' @param out character string, optional output filename. If NULL, only
#' formatted text is returned. default: NULL
#' @param sample_colname column name of the sample column, supplied as variable.
#' default: sample
#' @param group_colname column name of the sample column, supplied as variable.
#' default: group
#'
#' @return formatted text GSEA class file, if out is supplied, text will be
#' saved to the corresponding file
#' @export
#'
#' @importFrom dplyr select
#' @importFrom stringr str_replace_all
#'
#' @examples
create_gsea_cls <- function(data,
                            metadata,
                            out = NULL,
                            sample_colname = sample,
                            group_colname = group) {
  metadata_filtered <- metadata %>%
    dplyr::select({{ sample_colname }}, {{ group_colname }})

  # create string from supplied variables, needed for certain functions
  samples <- metadata_filtered[[deparse(substitute((sample_colname)))]]
  groups <- metadata_filtered[[deparse(substitute((group_colname)))]] %>%
    stringr::str_replace_all("[:space:]+", "")

  # which data columns and sample indices to include
  data_cols <- names(data)[names(data) %in% samples]
  sample_indices <- match(samples, data_cols)

  # compse the text of the final file
  # see https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
  text_lines <- c(
    paste(
      samples %>% unique() %>% length(),
      groups %>% unique() %>% length(),
      "1"
    ),
    paste("#", paste(groups %>% unique(), collapse = " ")),
    paste(groups[order(sample_indices)], collapse = " ")
  )

  if (!is.null(out)) {
    file_conn <- file(out)
    writeLines(text_lines, file_conn)
    close(file_conn)
  }
  return(text_lines)
}

#' Create gene ranking file for GSEA Preranked analysis
#'
#' @param data data frame containing gene expression data with a gene ID/symbol
#' column
#' @param out character string, optional output filename. If NULL, only
#' formatted text is returned. default: NULL
#' @param ranking_formula formula, way to calculate ranking of genes,
#' downregulated genes should have negative ranking, formula shouldn't produce
#' duplicates
#' @param id_colname column name with gene IDs/symbols in the data file,
#' supplied as variable
#' @param .inf one of ("drop", "replace"). How Inf values should be treated in
#' the data. "replace" will replace the Inf values with values 1% higher than
#' the current highest value
#'
#' @return formatted data frame with ranked genes for GSEA Preranked analysis,
#' if out is supplied, it will besaved to the corresponding file
#' @export
#'
#' @importFrom dplyr mutate select arrange filter desc pull
#' @importFrom lazyeval f_rhs
#' @importFrom tidyr drop_na
#'
#' @examples
create_gsea_rank <- function(data,
                             out = NULL,
                             ranking_formula = ~ -log10(pvalue) * sign(log2FoldChange),
                             id_colname = ensembl_gene_id,
                             .inf = "replace") {
  out_df <- data %>%
    dplyr::mutate(rank = !!lazyeval::f_rhs(ranking_formula)) %>%
    dplyr::select({{ id_colname }}, rank) %>%
    dplyr::arrange(dplyr::desc(rank)) %>%
    tidyr::drop_na()

  # process the +/- Inf values
  if (.inf == "drop") {
    out_df <- out_df %>%
      dplyr::filter(is.finite(rank))
  } else if (.inf == "replace") {
    # find min and max values that are not infinite
    min_val <- out_df %>%
      dplyr::filter(is.finite(rank)) %>%
      dplyr::pull(rank) %>%
      min()
    max_val <- out_df %>%
      dplyr::filter(is.finite(rank)) %>%
      dplyr::pull(rank) %>%
      max()

    # add 1% to the max value and replace the +/- Inf with new value
    out_df$rank[which(out_df$rank == Inf)] <- max_val + 0.01 * max_val
    out_df$rank[which(out_df$rank == -Inf)] <- min_val + 0.01 * min_val
  }

  # write to file if specified
  if (!is.null(out)) {
    write.table(out_df,
      out,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  return(out_df)
}

#' Return export parameters for a plot
#'
#' Paramters will contain properties such as filename, export device, dpi,
#' units, width and height of plot
#'
#' @param type character string, type of plots, corresponds to the graph_type
#' property in the config file
#' @param plot_params data frame or a character string pointing to a json file
#' with plot export parameters
#'
#' @return named list of plot export parameters
#' @export
#'
#' @importFrom dplyr filter select
#' @importFrom purrr map_dfr
#' @importFrom rjson fromJSON
#'
#' @examples
get_plot_params <- function(type, plot_params) {
  # check whether plot params are a data frame of path to config file
  if (class(plot_params) == "data.frame") {
    config_df <- plot_params
  } else if (class(plot_params) == "character" &&
    endsWith(plot_params, "json")) {
    config_df <- rjson::fromJSON(file = plot_params)$plot_export_params %>%
      purrr::map_dfr(data.frame)
  }

  # if the graph type is not present in config file R gives unhelpful error
  if (!(type %in% config_df$graph_type)) {
    stop("Graph type not present in parameter list for exporting plots.")
  }

  config_list <- config_df %>%
    dplyr::filter(graph_type == type) %>%
    dplyr::select(-graph_type) %>%
    as.list()
  return(config_list)
}

#' Save ggplot graph with specified parameters
#'
#' Parameters would be ideally supplied by a config json file
#'
#' @param output_dir character string, path to ouput directory, filename is
#' specified in plot_params
#' @param plot_params named list of export parameters for ggsave function such
#' as width, height, dpi, device etc. Includes the file name of the plot
#' @param plot plot to export. default: last plot
#' @param filename_prefix character string, add this as a prefix to the file
#' name
#' @param filename_suffix character string, add this as a suffix to the file
#' name
#' @param date_prefix logical, whether to include current date at the beginning
#' of the plot name
#'
#' @return NULL
#' @export
#'
#' @importFrom ComplexHeatmap draw
#' @importFrom rlang exec
#' @importFrom ggplot2 ggsave
#'
#' @examples
ggsave_param <- function(output_dir,
                         plot_params,
                         plot = last_plot(),
                         filename_prefix = "",
                         filename_suffix = "",
                         date_prefix = FALSE) {
  # assemble together the final file name
  if (date_prefix) {
    tp <- Sys.Date()
  } else {
    tp <- ""
  }
  out_filename <- paste0(
    output_directory,
    tp,
    filename_prefix,
    plot_params$filename,
    filename_suffix,
    ".",
    plot_params$device
  )

  # actual filename needs to be removed from the list otherwise its spliced into
  # the ggsave call and will cause an error
  plot_params$filename <- NULL

  # Complex Heatmap needs to be saved in specific manner
  if (!is.null(attr(class(plot), "package")) &&
      attr(class(plot), "package") == "ComplexHeatmap") {
    get(plot_params$device)(
      out_filename,
      width = plot_params$width,
      height = plot_params$height,
      units = plot_params$units,
      res = plot_params$dpi
    )
    ComplexHeatmap::draw(plot)
    dev.off()
  } else {
    rlang::exec(ggplot2::ggsave, out_filename, plot = plot, !!!plot_params)
  }
}

#' Wrapper around ggsave_param function
#'
#' To make saving graphs easier, only graph_type property from the config file
#' needs to be supplied. Functions assumes there is a file called 'config.json'
#' in the current directory and a output_directory global variable that points
#' to the locations graphs will be saved in.
#'
#' @param graph_type character string. Type of graph to be saved, needs to be a
#' property in the config file
#' @param ... other parameters passed to the ggsave_param function
#'
#' @return
#' @export NULL
#'
#' @examples
ggsave_param_wrapper <- function(graph_type, ...) {
  ggsave_param(
    output_directory,
    get_plot_params(graph_type, "config.json"),
    ...
  )
}
