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
library(cummeRbund)
library(DESeq2)
library(limma)
library(matrixStats)
library(broom)
library(stringr)
library(purrr)
library(dplyr)


#library(psych)
#library(clusterProfiler)
#library(ggfortify)

# upgrade of addFeatures from cummeRbund package, was using deprecated functions
.addFeatures<-function(object,features,level="genes",...){
  if(!is.data.frame(features)){
    stop("features must be a data.frame")
  }
  colnames(features)[1]<-slot(object,level)@idField
  colnames(features)<-make.db.names(object@DB,colnames(features),unique=T)
  dbWriteTable(object@DB,slot(object,level)@tables$featureTable,features,row.names=F,overwrite=T)
  indexQuery<-paste("CREATE INDEX ",slot(object,level)@idField," ON ", 
                    slot(object,level)@tables$featureTable," (",slot(object,level)@idField,")",sep="")
  res<-dbExecute(object@DB,indexQuery)
}
setMethod("addFeatures",signature(object="CuffSet"),.addFeatures)

# user defined reading functions -----------------------------------------------
read_cuffdiff_diff = function(dir) {
  cuff = readCufflinks(dir)
  annot = read.delim(paste0(dir,"/gene_exp.diff"), sep = "\t",header=T,na.string="-") %>% 
    dplyr::select(gene_id, gene)
  addFeatures(cuff,annot,level="genes")
  
  # munging of differential expression data
  diff = diffData(genes(cuff)) %>% as_tibble
  diff = diff %>% dplyr::filter(value_1 > 0 & value_2 > 0) %>% dplyr::filter(status == "OK")
  diff = diff %>% dplyr::mutate(comparison = paste0(sample_1, "_", sample_2))
  diff = merge(diff, annot %>% unique)
  return(diff)
}

merge_cuffdiff_fc = function(x, y, annot, by = "gene") {
  lfcs = merge(x %>% dplyr::select(gene, gene_id, log2_fold_change, q_value), 
               y %>% dplyr::select(gene, gene_id, log2_fold_change, q_value), 
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
  if(grepl(".xlsx", filename)) {
    data = read_excel(filename, sheet = sheet_name) %>% 
      select(-`#`)
  } else {
    data = read_delim(filename) %>%
      select(-`#`)
  }
  
  return(data)
}

# user defined plotting functions ----------------------------------------------
plot_pca_deseq = function(dds, condition = "condition") {
  vsd = DESeq2::vst(dds, blind = F)
  
  pcaData = DESeq2::plotPCA(vsd, intgroup=c(condition), returnData=TRUE)
  percentVar = round(100 * attr(pcaData, "percentVar"))
  segments = pcaData %>% 
    dplyr::group_by(!!as.symbol(condition)) %>% 
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pcaData = merge(pcaData, segments, by = condition)
  
  no_colors = pcaData[,condition] %>% unique() %>% length()
  
  p = pcaData %>%
    ggplot2::ggplot(aes(PC1, PC2, color=!!as.symbol(condition), shape=!!as.symbol(condition))) +
    geom_point(size=3) +
    geom_segment(aes(x = PC1, y = PC2, xend = xend, yend = yend), size = 0.5, linetype = "dashed") + 
    geom_point(data = segments, aes(x = xend, y = yend), size = 2) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    theme_bw() +
    scale_color_manual(values=viridis::viridis(no_colors+1)) +
    theme(legend.position = "top")
  
  return(p)
}

plot_pca_common = function(data, metadata, color_by = "group") {
  pca = prcomp(t(data), scale = T)
  pcaData <- as.data.frame(pca$x) %>% rownames_to_column("id")
  pcaData = merge(pcaData, metadata, by.x = "id", by.y = "sample")
  segments = pcaData %>% group_by(!!as.symbol(color_by)) %>% summarise(xend = mean(PC1), yend = mean(PC2))
  pcaData = merge(pcaData, segments, by = color_by)
  
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0(colnames(pca), "(", paste(as.character(percentage), "%", ")", sep=""))
  
  no_colors = pcaData[,color_by] %>% unique() %>% length()
  
  p = pcaData %>%
    ggplot(aes(PC1, PC2, color=!!as.symbol(color_by), shape=!!as.symbol(color_by))) +
    geom_point(size=3) +
    geom_segment(aes(x = PC1, y = PC2, xend = xend, yend = yend), size = 0.75) + 
    geom_point(data = segments, aes(x = xend, y = yend), size = 2) + 
    xlab(paste0("PC1: ",percentage[1]," variance")) +
    ylab(paste0("PC2: ",percentage[2]," variance")) + 
    theme_bw() +
    scale_color_manual(values=viridis(no_colors+1)[-length(viridis(no_colors+1))]) +
    theme(legend.position = "top")
  
  return(p)
}

plot_corr = function(gene_expr, param_values, param_values) {
  p = cbind.data.frame(expr = gene_expr, param = param_values, color_by = param_values) %>%
    ggplot(aes(x = expr, y = param, color = color_by)) +
    geom_point()
  p
}

plot_lfc_scatter = function(lfc_data) {
  p = lfc_data %>%
    dplyr::filter(p_chi < 0.1) %>%
    ggplot(aes(x = log2_fold_change.x, y = log2_fold_change.y, text = gene)) + 
    geom_point(aes(size = -log10(p_chi), color = err_sq)) + 
    scale_colour_continuous_sequential("viridis", rev = F) + 
    theme_bw() + 
    labs(size="-log10(p-value)", colour="squared\nerror")
  return(p)
}

plot_volcano_cuffdiff = function(data, sample1, sample2, customization, savename = NULL) {
  p = data %>% dplyr::filter(sample_1 == sample1 & sample_2 == sample2) %>% 
    ggplot(aes(x = log2_fold_change, y = -log10(p_value))) + 
    geom_point(aes(color = significant)) + 
    customization
  if(!is.null(savename)) {
    ggsave(savename,
           p,
           units = "mm",
           dpi = 600,
           width = 60, 
           height = 65)
  }
}

plot_venn2 = function(x, y, names, savename) {
  venn.diagram(
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

plot_dire = function(df) {
  p = df %>%
    ggplot(aes(x = Occurrence, y = Importance)) +
    geom_point() +
    theme_bw()
  return(p)
}

plot_dire_labeled = function(df, occurrence_threshold = 0.05, 
                             importance_threshold = 0.05, 
                             dot_size = 3.5) {
  p = df %>%
    plot_dire() +
    geom_label_repel(data = . %>% 
                       filter(Occurrence > occurrence_threshold | Importance > importance_threshold), 
                     aes(x = Occurrence, y = Importance, label = `Transcription Factor`),
                     size = dot_size)
  return(p)
}

# heatmap all
plot_heatmap = function(expression_data, metadata, palette = "RdBu") {
  expression_data %>%
    pheatmap(scale = "row",
             color = colorRampPalette(rev(RColorBrewer::brewer.pal(
               n = 7, name = palette)))(100),
             cluster_cols = T,
             show_rownames = F,
             cellwidth = 20,
             border_color = "white",
             treeheight_row = 15,
             treeheight_col = 20,
             annotation_col = metadata %>% 
               tibble::column_to_rownames("sample"))
  return(p)
}

plot_sample_heatmap = function(..., num_genes = 5000) {
  dots = list(...)
  expression_data = dots[1]
  metadata = dots[2]
  palette = dots[3]
  
  expr = expression_data %>%
    select(matches(metadata$sample)) %>% as.matrix
  
  p = plot_heatmap(expr[order(rowVars(expr)),][1:num_genes,],
                   metadata,
                   palette)
  return(p)
}

plot_diffexp_heatmap = function(..., geneid_colname, padj_colname) {
  dots = list(...)
  expression_data = dots[1]
  metadata = dots[2]
  palette = dots[3]
  
  p = expression_data %>%
    filter(!!as.symbol(geneid_colname) %in% (expression_data %>% 
                           filter (!!as.symbol(padj_colname) < 0.05) %>% 
                           pull(!!as.symbol(geneid_colname)))) %>%
    select(matches(metadata$sample)) %>%
    plot_heatmap(metadata, palette)
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
plot_volcano = function(fc_data, 
                        fc_colname = "log2FoldChange", 
                        pval_colname = "pvalue", 
                        padj_colname = "padj") {
  fc_data %>%
    mutate(padj = ifelse(is.na(!!as.symbol(padj_colname)), 1, !!as.symbol(padj_colname))) %>%
    mutate(log2FoldChange = !!as.symbol(fc_colname)) %>%
    mutate(significant = as.factor(case_when(padj < 0.05 & log2FoldChange > 0 ~ 1,
                                             padj < 0.05 & log2FoldChange < 0 ~ 2,
                                             TRUE ~ 0)))
  
  max_fc = max(fc_data$log2FoldChange)
  min_fc = min(fc_data$log2FoldChange)
  max_padj = max(-log10(padj))
  
  num_genes_up = fc_data %>% 
    filter(padj < 0.05 & log2FoldChange > 0) %>% 
    nrow()
  num_genes_down = fc_data %>% 
    filter(padj < 0.05 & log2FoldChange < 0) %>% 
    nrow()
  
  fc_data %>%
    ggplot(aes(x = log2FoldChange, 
               y = -log10(!!as.symbol(pval_colname)), 
               color = significant)) +
    geom_point(size = 2, alpha = 0.4) +
    geom_text(aes(x = max_fc, 
                  y = max_padj, 
                  label = paste0("\U1F815 ",num_genes_up)),
              color = "steelblue") + 
    geom_text(aes(x = min_fc, 
                  y = max_padj, 
                  label = paste0("\U1F817 ",num_genes_down)), 
              color = "darkred") +            
    scale_color_discrete(type = c("gray60", "darkred", "steelblue")) +
    theme_bw() + 
    theme(legend.position = "none")
}

plot_volcano_labeled = function(fc_data, gene_labels, symbol_colname = "SYMBOL", ...) {
  dots = list(...)
  fc_colname = dots[1]
  pval_colname = dots[2]
  padj_colname = dots[3]
  
  plot_volcano(fc_data, fc_colname, pval_colname, padj_colname) + 
    geom_label_repel(data = . %>% 
                       filter(!!as.symbol(symbol_colname) %in% gene_labels),
                     aes(label = !!as.symbol(symbol_colname)),
                     color = "black",
                     max.overlaps = 15)
}

# # heatmaps selected gene lists
# plot_gene_heatmap = function(df, metadata, gene_labels) {
#   df %>%
#     filter(SYMBOL %in% gene_labels) %>%
#     select(SYMBOL, matches(metadata$sample)) %>%
#     column_to_rownames("SYMBOL") %>%
#     pheatmap(scale = "row",
#              color = colorRampPalette(rev(RColorBrewer::brewer.pal(
#                n = 7, name = "RdBu")))(100),
#              cluster_cols = F,
#              show_rownames = F,
#              show_colnames = F,
#              cellwidth = 13,
#              cellheight = 13,
#              border_color = "white",
#              treeheight_row = 15,
#              annotation_col = metadata %>% column_to_rownames("sample"))
# }
# 
# plot_heatmap_fc = function(df, heatmap) {
#   labels = heatmap$tree_row$labels[heatmap$tree_row$order]
#   
#   df = df %>%
#     filter(SYMBOL %in% p$tree_row$labels[p$tree_row$order]) %>%
#     mutate(SYMBOL = factor(SYMBOL, levels = rev(p$tree_row$labels[p$tree_row$order])))
#   
#   df %>%
#     ggplot(aes(x = "gene", y = SYMBOL, fill = log2FoldChange, size = -log10(padj))) +
#     geom_point(shape = 21, color = "black") +
#     scale_fill_distiller(palette = "RdBu", limits = c(-1,1)*max(abs(df$log2FoldChange))) + 
#     theme_bw() + 
#     scale_y_discrete(position = "right") + 
#     xlab("") +
#     ylab("") + 
#     theme(axis.text.x = element_blank()) 
# }

# all_pathways = map_dfr(pathway_files, function(x) {
#   cp_unified_colnames = c("ID",	"Description",	"GeneRatio/NES",	"pvalue",	
#                           "p.adjust",	"SYMBOL",	"ENTREZID",	"log2FoldChange")
#   
#   read_delim(paste0(pathway_files_basepath, x),
#              delim = "\t", escape_double = FALSE,
#              trim_ws = TRUE)  %>% 
#     setNames(cp_unified_colnames) %>%
#     # TODO needs to be changed
#     mutate(source = gsub("_contrast.txt", "", x),
#            ID = as.character(ID)) %>%
#     # might not reflect actual columns
#     select(-SYMBOL, -ENTREZID, -log2FoldChange) %>%
#     unique()
# }) %>%
#   filter(p.adjust < 0.05)

collate_pathways = function(pathway_files_basepath, pattern) {
  pathway_files = list.files(pathway_files_basepath, 
                             pattern = pattern,
                             full.names = F)
  
  cp_unified_colnames = c("ID",	"Description",	"GeneRatio/NES",	"pvalue",	
                          "p.adjust",	"SYMBOL",	"ENTREZID",	"log2FoldChange")
  
  diffexp_pathways = map_dfr(pathway_files, function(x) {
    read_delim(paste0(pathway_files_basepath, x),
               delim = "\t", escape_double = FALSE,
               trim_ws = TRUE)  %>% 
      setNames(cp_unified_colnames) %>%
      # TODO needs to be changed
      mutate(source = gsub("_contrast.txt", "", x),
             ID = as.character(ID)) %>%
      # might not reflect actual columns
      select(-SYMBOL, -ENTREZID, -log2FoldChange) %>%
      unique()
  }) %>%
    filter(p.adjust < 0.05)
  
  return(diffexp_pathways)
}


plot_pathways_meta = function(df, top_pathways = 30) {
  pathways_summary = all_pathways %>%
    group_by(source) %>%
    summarise(n_pathways = n(), mean_enrichment = mean(`GeneRatio/NES`)) 
  
  pathways_summary %>%
    ggplot() +
    geom_histogram(aes(x = n_pathways), bins = 100, fill = "steelblue") + 
    theme_bw()
  
  # top pathway contributors
  pathways_summary %>%
    top_n(30, n_pathways) %>%
    arrange(n_pathways) %>%
    mutate(source = factor(source, levels = source)) %>%
    ggplot(aes(x = n_pathways, y = source)) +
    geom_bar(stat="identity", fill = "steelblue") + 
    theme_bw()
  
  # bottom pathway contributors
  pathways_summary %>%
    top_n(30, -n_pathways) %>%
    arrange(n_pathways) %>%
    mutate(source = factor(source, levels = source)) %>%
    ggplot(aes(x = n_pathways, y = source)) +
    geom_bar(stat="identity", fill = "steelblue") + 
    theme_bw()
}


plot_pathway_bargraph = function(data, pathway_source, top_n = 20, truncate_desc = 80) {
  data %>%
    filter(grepl(pathway_source, source)) %>%
    arrange(-abs(`GeneRatio/NES`)) %>%
    slice_head(n = top_n) %>%
    arrange(`GeneRatio/NES`) %>%
    mutate(Description = stringr::str_trunc(Description, truncate_desc)) %>%
    mutate(Description = factor(Description, levels = (.) %>% pull(Description))) %>%
    ggplot(aes(x = `GeneRatio/NES`, 
               y = Description, 
               color = -log10(p.adjust), 
               fill = -log10(p.adjust))) +
    geom_bar(stat="identity") + 
    ylab("") +
    theme_bw()
}


gene_lists = list(
  myh3 = c("Myh3"),
  tnn = c("Tnnc1", "Tnnc2", "Tnni1", "Tnni2", "Tnni3", "Tnnt1", "Tnnt2", "Tnnt3", "Tpm1", "Tpm2", "Tpm3", "Tpm4"),
  mito_strict = c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"),
  mito_assoc = c("Polrmt", "Tfam", "Tfb1m", "Tfb2m", "Mterf1a", "Mterf1b", "Mterf2", "Mterf3", "Mterf4", "Trmt10c", 
                 "Hsd17b10", "Prorp", "Elac2", "Fastk", "Fastkd1", "Fastkd2", "Fastkd4", "Lrpprc", "Slirp"),
  metabolism = c("Ppargc1a", "Mb", "Myog", "Mstn", "ND5", "Cyc1", "Sdha", "Cox5a", "Atp5a1")
)

pathway_lists = list(
  wp = c("WP1248", "WP295", "WP662", "WP434", "WP85", "WP6", "WP2185", "WP523", "WP113")
)

sample_expression <- read_delim("~/OneDrive/_Ruaslab/01_projects/01_hectd1/19_rnaseq/diffexp/deseq_default/degfiles/sample_expression.csv",
                                delim = "\t", escape_double = FALSE,
                                trim_ws = TRUE)

ko_vs_wt_all <- read_delim("~/OneDrive/_Ruaslab/01_projects/01_hectd1/19_rnaseq/diffexp/deseq_default/degfiles/ko_vs_wt_all.txt",
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE)

pathways <- read_delim("D:/igor/labmeeting_20211117/output/csv/ko_wt/Pathways_Wiki_GSEA_all_contrast.txt",
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)

pathways2 <- read_delim("D:/igor/labmeeting_20211117/output/csv/ko_wt/Pathways_Wiki_ORA_all_contrast.txt",
                        delim = "\t", escape_double = FALSE,
                        trim_ws = TRUE)
