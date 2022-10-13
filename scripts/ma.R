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
