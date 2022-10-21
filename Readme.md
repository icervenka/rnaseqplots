# RNASeq visualization 

Collection of scripts for creating publication ready graph for bulk RNASeq and pathway analysis. Contains several types of graph with included default parameters to simplify graph creation and export. Currently this project exists as a collection scripts and a run file which will contain your code to generate graphs. In the future, this will be transformed into R package instead.

## How to use
For now the easiest way is to clone the project, change the entries in `paths.json` file to point to the relevant data and desired output directory, specify analysis and export parameters in the `config.json` file, place the code for the graphs in the `run.R` file and execute. Supplied `run.R` file contains examples for all graph types as well as their export.

## Input
Loacation of input files is specified in `paths.json` file either as a relative or an absolute path. Scripts accept both comma- and tab-delimited text files and xlsx files in some instances. See example `data` for more details.

### Expression data
Must contain some gene ID column (Ensemble, Entrez or Symbol) and one column of normalized gene expression for each sample. Column names with samples have to match the metadata.

### Differential expression data
Should contain gene ID column, some measure of expression level for MA plots (eg. baseMean from DESeq), a column for log2(fold-change) and columns with raw and adjusted pvalues for volcano plots, heatmaps etc.

### Metadata
Must contain columns for sample names and group

### DESeq2 dds objects
Exported dds object with `saveRDS` after running `DESeq()` function.

### Cuffdiff analysis
Path to folder where `cufdiff` results are stored.

### GSEA pathway analysis
Results of GSEA pathway analysis, script only makes use of summary results present in `gsea_report_for_*.tsv` files. Files can be organized in subdirectories, the script will scan them recursively and will add the directory path as a column in the final data frame.

### Clusterprofiler pathway analysis
This scripts processes `csv` output files from [clusterProfiler report generator](https://github.com/icervenka/clusterprofiler_reports_snakemake).

### Dire (dcode.org) transcription factor analysis
Processing of output from [dire](https://dire.dcode.org/) transcription factor analysis. Accepts also xlsx files. Must contain `Occurence` and `Importance` columns. 

### Qiagen IPA analysis
Processes ouput from [IPA report generator](https://github.com/icervenka/ipa_reports_snakemake). To learn more about IPA please visit [Qiagen website](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/).

### Gene lists
Json file containing custom gene list to be used in plots throughout the analysis. Gene names/IDs are case sensitive and must match the appropriate column in your expression data or differential expression data.

### Pathway lists
Json file containg lists of patways to plot in pathway barplots or scatterplots. Pathways can be either specified by their full name (case sensitive), or by rank they were assigned during analysis.

## Output

## Parameters

## Supported Graphs

### PCA plots

### MA plot

### Volcano plots

### Heatmaps

### Venn diagrams

### Dire transcription factor plots

### Pathway analysis plots

## Required packages

List of required packages can be found in `scripts/load_packages.R` file. 

## Issues

