#' Upgrade of addFeatures from cummeRbund package
#'
#' It was using some deprected function
#'
#' @param object aaa
#' @param features bbb
#' @param level ccc
#' @param ... ddd
#'
#' @return dbi
#'
#' @examples
.addFeatures <- function(object, features, level = "genes", ...) { # nolint
  if (!is.data.frame(features)) {
    stop("features must be a data.frame")
  }
  colnames(features)[1] <- slot(object, level)@idField
  colnames(features) <- make.db.names(object@DB,
    colnames(features),
    unique = TRUE
  )
  DBI::dbWriteTable(object@DB,
    slot(object, level)@tables$featureTable,
    features,
    row.names = FALSE,
    overwrite = TRUE
  )
  indexQuery <- paste( # nolint
    "CREATE INDEX ",
    slot(object, level)@idField,
    " ON ",
    slot(object, level)@tables$featureTable,
    " (",
    slot(object, level)@idField,
    ")",
    sep = ""
  )
  res <- DBI::dbExecute(object@DB, indexQuery)
}
# TODO fix
# setMethod("addFeatures", signature(object = "CuffSet"), .addFeatures)

# user defined reading functions -----------------------------------------------

#' Universal wrapper for reading tabular text data
#'
#' Determines the correct function for reading based on extension
#' (csv, tsv, txt).
#'
#' @param filename character string, filename with tabular data to read.
#'
#' @return data frame
#' @export
#'
#' @examples
read_data <- function(filename) {
  ext <- tools::file_ext(filename)
  if (ext %in% c("tsv", "txt")) {
    separator <- "\t"
  } else if (ext %in% c("csv")) {
    separator <- ","
  }
  df <- readr::read_delim(filename,
    delim = separator,
    col_names = TRUE,
    escape_double = FALSE,
    trim_ws = TRUE
  )
  return(df)
}

#' Read gene list from single column text file
#'
#' File shouldn't have a header.
#'
#' @param filename text file to read that contains single column with gene names
#' or IDs.
#'
#' @return character vector
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
read_gene_list_from_file <- function(filename) {
  gene_list <- readr::read_delim(filename,
    delim = "\t",
    col_names = c("SYMBOL")
  ) %>%
    dplyr::pull(SYMBOL)
  return(gene_list)
}

#' Read cufdiff data
#'
#' @param dir character string, directory path that contains cuffdiff data.
#'
#' @return data frame
#' @export
#'
#' @examples
read_cuffdiff <- function(dir) {
  cuff <- cummeRbund::readCufflinks(dir)
  annot <- read.delim(paste0(dir, "/gene_exp.diff"),
    sep = "\t",
    header = TRUE,
    na.string = "-"
  ) %>%
    dplyr::select(gene_id, gene)
  cummeRbund::addFeatures(cuff, annot, level = "genes")

  # munging of differential expression data
  diff <- cummeRbund::diffData(cummeRbund::genes(cuff)) %>%
    tibble::as_tibble()
  diff <- diff %>%
    dplyr::filter(value_1 > 0 & value_2 > 0) %>%
    dplyr::filter(status == "OK")
  diff <- diff %>%
    dplyr::mutate(comparison = paste0(sample_1, "_", sample_2)) %>%
    dplyr::left_join(annot %>% unique())
  return(diff)
}

#' Read dire data from excel file
#'
#' Result data from dire.dcode can be copied to excel while preserving tabular
#' formatting. The data should be kept as is from the webpage, the function
#' will take care of proper formatting.
#'
#' @param filename character string, file to read.
#' @param sheet_name character string, sheet name that contains the data.
#' default: Sheet1
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
read_dire_xlsx <- function(filename, sheet_name = "Sheet1") {
  if (grepl(".xlsx", filename)) {
    df <- readxl::read_excel(filename, sheet = sheet_name) %>%
      dplyr::select(-`#`)
  } else {
    df <- readr::read_delim(filename) %>%
      dplyr::select(-`#`)
  }
  return(df)
}

#' Batch reading of dire data from directory
#'
#' Function will attempt to read all text and excel files in directory that
#' match the specified pattern. If other types of files match, it will throw
#' an error.
#'
#' @param pathway_files_basepath character string, directory where dire data is
#' located.
#' @param pattern character string, regular expression to match the filenames
#' in specified directory, empty string will match all files.
#' default: empty string
#' @param sheet_name character string, for excel files the name of the sheet
#' that contains the dire data. default: "Sheet1"
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
collate_dire_pathways <- function(pathway_files_basepath,
                                  pattern = "",
                                  sheet_name = "Sheet1",
                                  recursive = TRUE) {
  empty_dire_df <- data.frame(
    `Transcription Factor` = character(0),
    Occurrence = numeric(0),
    Importance = numeric(0)
  )

  pathway_files <- list.files(
    pathway_files_basepath,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive
  )

  dire_pathways <- purrr::map_dfr(pathway_files, function(x) {
    if (endsWith(x, "xlsx")) {
      tryCatch(
        {
          dire_df <- read_dire_xlsx(x, dire_sheet_name)
        },
        error = function(cond) {
          message("There were problems dire xlsx file, please check that the
          specifications are correct, skipping.")
          message(cond)
          return(empty_dire_df)
        }
      )
    } else {
      tryCatch(
        {
          dire_df <- read_data(x) %>%
            dplyr::select(-`#`) %>%
            mutate(Occurrence = as.numeric(gsub("%", "", Occurrence)) / 100)
        },
        error = function(cond) {
          message("There were problems dire text file, please check that the
          specifications are correct, skipping.")
          message(cond)
          return(empty_dire_df)
        }
      )
    }
    dire_df <- dire_df %>%
      dplyr::mutate(file = x)
  })
  return(dire_pathways)
}

#' Parse json file generated by AmiGO analysis
#'
#' Jsonfile generated by http://amigo.geneontology.org/amigo preserves pathway
#' hierarchy. Resulting data frame will have a level column indicating
#' hierarchy level.
#'
#' @param jsonfile json like list, generated by rjson::fromJSON function.
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
read_amigo_json <- function(jsonfile) {

  # If the json file specification changes, this needs to be updated
  df <- purrr::map_dfr(jsonfile$overrepresentation$group, function(y) {
    if (length(y$result) == 3 &
      !is.list(y$result[[2]])) {
      data.frame(t(c(
        id = y$result$term$id,
        label = y$result$term$label,
        fold_enrichment = y$result$input_list$fold_enrichment,
        pvalue = y$result$input_list$pValue,
        level = y$result$term$level
      )))
    } else {
      purrr::map_dfr(y$result, function(x) {
        data.frame(t(c(
          id = x$term$id,
          label = x$term$label,
          fold_enrichment = x$input_list$fold_enrichment,
          pvalue = x$input_list$pValue,
          level = x$term$level
        )))
      })
    }
  }) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      fold_enrichment = as.double(fold_enrichment),
      pvalue = as.double(pvalue),
      level = as.integer(level)
    )
  return(df)
}

#' Read GSEA result file
#'
#' Adds a filename column containing the filename parameter.
#'
#' @param filename character string, text file with GSEA pathways to process
#' @param rank_by which column should the data be ranked on, supplied as
#' variable. default: NES
#' @param descending logical, whethere the renking should be from highest to
#' lowest
#' @param rename_columns character vector of size 11, new names for data frame
#' columns. Applied as a last operation. default: c("name", "link", "details",
#' "size", "es", "nes", "nom_pval", "fdr_qval", "fwer_pval", "rank_at_max",
#' "leading_edge")
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
read_gsea <- function(filename,
                      rank_by = NES,
                      descending = TRUE,
                      rename_columns = c(
                        "name", "link", "details", "size", "es", "nes",
                        "nom_pval", "fdr_qval", "fwer_pval", "rank_at_max",
                        "leading_edge"
                      )) {
  orientation <- rank_how(descending)

  df <- read.table(filename, header = TRUE, sep = "\t") %>%
    dplyr::select(-X) %>%
    dplyr::mutate(file = filename) %>%
    dplyr::arrange(orientation(abs({{ rank_by }}))) %>%
    dplyr::mutate(pathway_rank = dplyr::row_number()) %>%
    # small e is added to avoid problems in volcano plots with log transforms
    dplyr::mutate(`NOM p-val` = `NOM p-val` + 0.00001) %>%
    dplyr::mutate(`FDR q-val` = `FDR q-val` + 0.00001) %>%
    setNames(rename_columns)
  return(df)
}

#' Recursively find GSEA result files in a directory
#'
#' @param dirpath character string, directory path where the result files are
#' located.
#' @param pattern character string, regex pattern restricting time file names
#' that will be processed.
#'
#' @return data frame containg the locations of result files and their paths
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
get_gsea_files_from_dir <- function(dirpath, pattern, recursive = TRUE) {
  file_df <- data.frame(full_path = dir(dirpath,
    pattern = paste0(pattern, ".*tsv"),
    full.names = TRUE,
    recursive = recursive
  )) %>%
    tidyr::separate(full_path,
      into = c("path", "file"),
      sep = paste0("/", pattern),
      remove = FALSE
    ) %>%
    dplyr::mutate(file = paste0(pattern, file))
  return(file_df)
}

#' Filter pathways in a data frame
#'
#' @param pathway_data data frame with data for pathway analysis. Should include
#' a columns for enrichment/score, p-value and a pathway name
#' @param score_column name of score column, supplied as variable.
#' @param pvalue_column name of p-value column, supplied as variable.
#' @param rank_by which column to rank the data on, supplied as variable.
#' @param descending logical, whether ranking should be done from highest to
#' lowest value. default: TRUE
#' @param score_threshold double, return data subset that pass the score
#' threshold. Will compare the threshold against absolute values from score
#' column. default: 0
#' @param pvalue_threshold double, return data subset that pass the p-value
#' threshold.
#'
#' @return filtered pathway_data data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
filter_pathways <- function(pathway_data,
                            score_column,
                            pvalue_column,
                            rank_by,
                            descending = TRUE,
                            score_threshold = 0,
                            pvalue_threshold = 0.05) {
  orientation <- rank_how(descending)

  pathway_data %>%
    dplyr::filter(abs({{ score_column }}) > score_threshold) %>%
    dplyr::filter({{ pvalue_column }} < pvalue_threshold) %>%
    dplyr::arrange(orientation(abs({{ rank_by }}))) %>%
    dplyr::mutate(pathway_rank = dplyr::row_number())
}


#' Batch read GSEA files into one data frame
#'
#' Combines get_gsea_files_from_dir, read_gsea and filter_pathways functions and
#' sets the column and threshold defaults.
#'
#' @param dir character string, directory with GSEA results to process. Will be
#' processed recursively. default: . (current directory)
#' @param pattern character string, regex pattern for filenames to include.
#' default: gsea_report_for
#' @param score_column name of score column, supplied as variable. default: nes
#' @param score_threshold double, return data subset that pass the score
#' threshold. Will compare the threshold against absolute values from score
#' column. default: 0
#' @param pvalue_column name of p-value column, supplied as variable. default:
#' fdr_qval
#' @param pvalue_threshold double, return data subset that pass the p-value
#' threshold.
#' @param rank_by which column to rank the data on, supplied as variable.
#' default: nes
#'
#' @return data frame with parsed GSEA data
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
batch_read_filter_gsea <- function(dir = ".",
                                   pattern = "gsea_report_for",
                                   score_column = nes,
                                   score_threshold = 0,
                                   pvalue_column = fdr_qval,
                                   pvalue_threshold = 0.05,
                                   rank_by = nes,
                                   recursive = TRUE) {
  # read the data frame of files
  files_df <- get_gsea_files_from_dir(dir, pattern, recursive = recursive)

  # process files
  purrr::map_dfr(unique(files_df$path), function(x) {
    filtered_df <- files_df %>%
      dplyr::filter(path == x)
    purrr::map_dfr(filtered_df$file, function(y) {
      read_gsea(paste0(x, "/", y))
    }) %>%
      filter_pathways(
        score_column = score_column,
        pvalue_column = pvalue_column,
        rank_by = rank_by,
        score_threshold = score_threshold,
        pvalue_threshold = pvalue_threshold
      )
  })
}

#' Read analysis file from ipa_reports_snakemake
#'
#' @param filename character string, text file with IPA pathways to process
#' @param rank_by which column to rank the data on, supplied as variable.
#' default: zscore
#' @param descending logical, whether ranking should be done from highest to
#' lowest value. default: TRUE
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
read_ipa <- function(filename, rank_by = zscore, descending = TRUE) {
  orientation <- rank_how(descending)

  read.table(filename, header = TRUE, sep = "\t") %>%
    dplyr::arrange(orientation(abs({{ rank_by }}))) %>%
    dplyr::mutate(pathway_rank = dplyr::row_number())
}

#' Batch read IPA files into one data frame
#'
#' @param pathway_files_basepath character string, directory path where the
#' result files are located.
#' @param pattern character string, regex pattern for filenames to include.
#' @param rank_by which column to rank the data on, supplied as variable.
#' default: zscore
#' @param descending logical, whether ranking should be done from highest to
#' lowest value. default: TRUE
#'
#' @return data frame with IPA pathways
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
collate_ipa_pathways <- function(pathway_files_basepath,
                                 pattern,
                                 rank_by = zscore,
                                 descending = TRUE,
                                 recursive = TRUE) {
  pathway_files <- list.files(
    pathway_files_basepath,
    pattern = pattern,
    full.names = FALSE,
    recursive = recursive
  )

  ipa_pathways <- purrr::map_dfr(pathway_files, function(x) {
    # TODO add df type to ipa snakemake pipeline
    read_ipa(x, rank_by = rank_by, descending = descending) %>%
      dplyr::mutate(file = x)
  })
  return(ipa_pathways)
}

#' Batch read files from  clusterprofiler_reports_snakemake into one data frame
#'
#' @param pathway_files_basepath character string, directory path where the
#' result files are located.
#' @param pattern character string, regex pattern for filenames to include.
#' @param padj_threshold double, p-value significance threshold fo filtering
#' pathways
#'
#' @importFrom magrittr %>%
#'
#' @return data frame
#' @export
#'
#' @examples
collate_cp_pathways <- function(pathway_files_basepath,
                                pattern,
                                pval_threshold = 0.05,
                                recursive = TRUE) {
  pathway_files <- list.files(pathway_files_basepath,
    pattern = pattern,
    full.names = FALSE,
    recursive = recursive
  )

  cp_unified_colnames <- c(
    "ID", "Description", "GeneRatio/NES", "pvalue",
    "p.adjust", "SYMBOL", "ENTREZID", "log2FoldChange"
  )

  diffexp_pathways <- purrr::map_dfr(pathway_files, function(x) {
    readr::read_delim(paste0(pathway_files_basepath, x),
      delim = "\t", escape_double = FALSE,
      trim_ws = TRUE,
      show_col_types = FALSE
    ) %>%
      setNames(cp_unified_colnames) %>%
      dplyr::mutate(ID = as.character(ID)) %>%
      # TODO might not reflect actual columns, check with cp reports package
      dplyr::select(-SYMBOL, -ENTREZID, -log2FoldChange) %>%
      unique()
  }) %>%
    dplyr::filter(p.adjust < pval_threshold)
  return(diffexp_pathways)
}
