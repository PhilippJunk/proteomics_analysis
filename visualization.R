library(tidyverse)
library(patchwork)
library(broom)
library(UpSetR)

# ideas:
# Volcano plots
# heatmap?

###############################################################################
vis_pca <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  df_wide <- p_df %>%
    pivot_wider(names_from = id, values_from = LFQ, values_fill = 0)
  
  pca <- df_wide %>%
    select(where(is.numeric)) %>%
    select(where(~ sd(.x) > 0)) %>%
    prcomp(scale=T, center=T)
  
  pc1_varexpl <- round(summary(pca)$importance[2,1] * 100, 2)
  pc2_varexpl <- round(summary(pca)$importance[2,2] * 100, 2)
  
  pca %>%
    broom::augment(df_wide) %>%
    inner_join(attr(p_df, 'annotation'), by= 'label') %>%
    ggplot(aes(x=.fittedPC1, y=.fittedPC2, color=group)) +
      geom_point() +
      labs(x = str_c('PC1 (', pc1_varexpl, '%)'),
           y  = str_c('PC2 (', pc2_varexpl, '%)')) +
      theme_bw() +
      NULL
}

###############################################################################
# TODO
# vis_umap <- function(p_df) {
#   # check inputs
#   p_df <- validate_proteomics_data(p_df)
# }

###############################################################################
vis_qc_scatter <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  labels <- attr(p_df, 'annotation') %>%
    pull(label) %>%
    unique %>%
    sort
  
  expand_grid(l1 = labels, l2 = labels) %>%
    pmap(function(l1, l2) {
      df1 <- p_df %>% 
        filter(label == l1) %>%
        rename(LFQ_1 = LFQ)
      df2 <- p_df %>%
        filter(label == l2) %>%
        rename(LFQ_2 = LFQ)
      inner_join(df1, df2, by='id') %>%
        ggplot(aes(x = LFQ_2, y = LFQ_1)) +
          geom_point() +
          geom_abline(slope=1, intercept=0) +
          labs(x = l2, y = l1)
    }) %>%
    wrap_plots(ncol = length(labels), nrow = length(labels), byrow=T)
}

###############################################################################
vis_qc_histo <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  ggplot(p_df, aes(x=LFQ, fill=label)) +
    geom_histogram() +
    facet_wrap(label ~ .) +
    theme(legend.position = 'none')
}

###############################################################################
vis_qc_count <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  p_df %>%
    count(label, name='count') %>% 
    ggplot(aes(y = label, x = count)) +
      geom_bar(stat='identity', color='black') +
      scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
}

###############################################################################
vis_upset <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  annotation <- attr(p_df, 'annotation')
  
  nsets <- annotation$group %>% unique %>% length
  nintersects = 2^nsets
  
  p_df %>%
    inner_join(annotation, 'label') %>%
    select(group, id) %>%
    distinct %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = group, values_from = present, values_fill = 0) %>% 
    as.data.frame %>% 
    upset(sets = annotation %>% pull(group) %>% unique, 
          order.by = "freq",
          nsets = nsets,
          nintersects = nintersects)
}

vis_upset_ch <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  annotation <- attr(p_df, 'annotation')
  
  combination_matrix <- inner_join(
    p_df, annotation, by = 'label') %>%
    group_by(group) %>%
    nest %>%
    mutate(data = set_names(data, group)) %>%
    pull(data) %>%
    map(function(df) {
      df %>%
        pull(id) %>%
        unique
    }) %>%
    ComplexHeatmap::make_comb_mat()
  
  ComplexHeatmap::UpSet(
    combination_matrix, 
    comb_order = order(-ComplexHeatmap::comb_size(combination_matrix)))
}

