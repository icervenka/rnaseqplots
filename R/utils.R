#' Checks is item is present in a configuration dictionary
#'
#' @param entry character string to check
#' @param config_dict named list usually obtained by parsing json configuration
#' file
#'
#' @return TRUE or FALSE depending whether the entry is present as a name in the
#' list
#' @export
#'
#' @examples
validate_config <- function(entry, config_dict) {
  # empty strings are sometimes needed for patterns for file listing
  if (entry %in% names(config_dict)) { # && nchar(config[[item]]) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Adds a slash at the end of path
#'
#' Slash is added only if there is none present.
#'
#' @param dir_string character string representing path to directory
#'
#' @return path too directory with slash as a last character
#'
#' @export
#'
#' @examples
append_dir_slash <- function(dir_string) {
  if (!base::endsWith(dir_string, "/")) {
    dir_string <- paste0(dir_string, "/")
  }
  return(dir_string)
}

#' Check if directory exists
#'
#' @param dir_string character string of path to check
#'
#' @return character string of path if it exists, stops executions otherwise
#'
#' @export
#'
#' @examples
dir_exists <- function(dir_string) {
  if (!file.exists(dir_string)) {
    stop(paste0("Specified directory does not exist: ", dir_string))
  } else {
    return(dir_string)
  }
}

#' Create default ggplot2 color hues
#'
#' @param i integer, number of colors to generate
#'
#' @return character vector with equally spaced color hues
#' @export
#'
#' @examples
gg_color_hue <- function(i) {
  hues <- seq(15, 375, length = i + 1)
  hcl(h = hues, l = 65, c = 100)[1:i]
}

#' Helper switch logical toggle for descending ranking function
#'
#' @param descending logical, whether to rank data frame columns in descending
#' order
#'
#' @return function  for ranking data frame entries
#' @export
#'
#' @examples
rank_how <- function(descending) {
  if (descending == TRUE) {
    return(dplyr::desc)
  } else {
    return(identity)
  }
}

#' Text processing of gene symbols
#'
#' @param gene_list character vector of gene symbols to process
#' @param text_replacements named vector, where names are patterns and values
#' replacement. Will be applied in order it is specified.
#' @param text_modifications list of single argument text processing functions
#' to apply to gene_list. Will be applied from right to left.
#'
#' @return character vector of modified gene names
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
modify_gene_symbols <- function(gene_list,
                                text_replacements,
                                text_modifications) {
  gene_list <- gene_list %>%
    stringr::str_replace_all(regex(text_replacements, ignore_case = TRUE)) %>%
    purrr::map_chr(purrr::compose(!!!text_modifications))
  return(gene_list)
}

#' Replace mouse gene symbols
#'
#' Applies simple text replacements for mitochondrial gene symbols  and converts
#' gene symbols to title case. This is a more specific version of
#' modify_gene_symbols function and a fast way to replace human gene symbols
#' with reasonable accuracy.
#'
#' @param gene_list character vector of gene symbols to process
#'
#' @return character vector of modified gene names
#' @export
#'
#' @examples
mouseify_gene_symbols <- function(gene_list) {
  modify_gene_symbols(gene_list,
    text_replacements = c(
      mito_atp_replacer,
      c("MT-CO" = "MT-COX", "MT-" = "")
    ),
    text_modifications = list(tools::toTitleCase, tolower)
  )
}
