% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_pathways_ipa}
\alias{read_pathways_ipa}
\title{Read analysis file from ipa_reports_snakemake}
\usage{
read_pathways_ipa(filename, rank_by = zscore, descending = TRUE)
}
\arguments{
\item{filename}{character string, text file with IPA pathways to process}

\item{rank_by}{which column to rank the data on, supplied as variable.
default: zscore}

\item{descending}{logical, whether ranking should be done from highest to
lowest value. default: TRUE}
}
\value{
data frame
}
\description{
Read analysis file from ipa_reports_snakemake
}
