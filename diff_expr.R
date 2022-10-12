library(tidyverse)
library(limma)

# TODO maybe include GSEA analysis in here?

###############################################################################
#' The available methods for statistical analysis.
#' 
#' @return A character vector with the available methods for statistical
#' analysis
stat_analysis_methods <- function() {
  c('limma', 'ttest')
}

#' Differential analysis of protemic data set with different methods.
#' 
#' @param p_df The data to run differential statistical analysis on. 
#'   Of class `proteomics_data`.
#' @param method The method used for differential analysis. See
#'   `stat_analysis_methods()` for available options.
#' @param contrasts The contrasts between condition to test in differential
#'   analysis. Expected to be a data frame with two columns, a and b, where 
#'   each row corresponds to a contrast to test.
#' @param adjust_method The method used for correcting for multiple hypothesis
#'  testing. See `p.adjust.methods` for available options.
#' @return A data frame with calculated log fold changes and adjusted pvalues
#'   for each specified contrast.
diff_expr <- function(
    p_df,
    method,
    contrasts,
    adjust_method = 'BH') {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  # extract from input data
  df_wide <- p_df %>%
    pivot_wider(names_from = 'label', values_from = 'LFQ', values_fill = NA)
  mx <- df_wide %>%
    select(where(is.numeric)) %>%
    as.matrix
  rownames(mx) <- df_wide %>% 
    select(!where(is.numeric)) %>% 
    pull
  annotation <- attr(p_df, 'annotation')
  # get groups from data and annotation
  cols <- mx %>% colnames
  groups <- annotation$group[order(match(annotation$label, cols))]
  
  if (method == 'limma') {
    # transform contrasts into limma format
    groups <- factor(make.names(groups))
    contrasts <- contrasts %>%
      transmute(str_glue('{make.names(a)} - {make.names(b)}')) %>%
      pull
    # run analysis
    out <- .stat_analysis_limma(mx, groups, contrasts)
  } else if (method == 'ttest') {
    NULL
  }
  
  # adjust pvalues and return output
  out %>%
    mutate(pval_adj = p.adjust(pval, adjust_method))
}


#' Differential analysis using LIMMA
#' 
#' @param mx An expression matrix.
#' @param groups A vector indicating which columns belong to which group
#' @param contrasts A vector of contrasts in limma format.
#' @return A data frame with raw pvalues.
.stat_analysis_limma <- function(
    mx,
    groups,
    contrasts) {
  # create design matrix and fix colnames
  design_mx <- stats::model.matrix(~ 0 + groups)
  colnames(design_mx) <- colnames(design_mx) %>%
    str_replace(deparse(substitute(groups)), '')
  # create contrast matrix
  contrast_mx <- limma::makeContrasts(contrasts = contrasts, levels = design_mx)
  # run limma
  fit <- limma::lmFit(mx, design_mx)
  contr <- limma::contrasts.fit(fit, contrast_mx)
  bayes <- limma::eBayes(contr)
  # parse output
  colnames(contrast_mx) %>%
    map(function(contr) {
      limma::topTable(bayes, coef=contr, number=Inf) %>%
        tibble::rownames_to_column('id') %>%
        mutate(contrast = contr)
    }) %>% 
    bind_rows %>%
    rename(logfc = logFC,
           pval = P.Value) %>%
    select(contrast, id, logfc, pval)
}


#' Differential analysis using T-tests.
#' 
#' @param mx An expression matrix.
#' @param groups A vector indicating which columns belong to which group
#' @param contrasts The contrasts between condition to test in differential
#'   analysis. Expected to be a data frame with two columns, a and b, where 
#'   each row corresponds to a contrast to test.
#' @return A data frame with raw pvalues.
.stat_analysis_ttest <- function(
    mx, 
    groups, 
    contrasts) {
  # for each contrast
  contrasts %>%
    pmap(function(a, b) {
      # extract groups from matrix
      mx1 <- mx[,which(groups == a)]
      mx2 <- mx[,which(groups == b)]
      # run tests
      p_values <- 1:nrow(mx1) %>%
        map(function(n) {
          # check if enough observations
          if (sum(!is.na(mx1[n,])) < 2 & sum(!is.na(mx2[n,])) < 2) {
            return(NA_real_)
          }
          tryCatch({
            stats::t.test(mx1[n,], mx2[n,]) %>%
              glance %>%
              pull(p.value)
          }, error = function(e) { return(NA_real_) })
        }) %>%
        flatten_dbl
      
      # get logfcs
      logfc <- unname(rowMeans(mx1, na.rm = T) - rowMeans(mx2, na.rm = T))
      # assemble output
      tibble(contrast = str_glue('{make.names(a)} - {make.names(b)}'),
             id = rownames(mx),
             logfc = logfc,
             pval = p_values)
    }) %>%
    bind_rows
}

###############################################################################
# helper functions

# construct contrasts against control group
construct_contrasts_control <- function(p_df, cntrl_group) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  all_groups <- unique(attr(p_df, 'annotation')$group)
  
  if (!cntrl_group %in% all_groups) {
    stop('Control group not found in protemic data set.',
         call. = FALSE)
  }
  
  tibble(a = all_groups[!all_groups == cntrl_group],
         b = cntrl_group)
}

###############################################################################
# helper function

# construct info whether comparisons are based on value-value, value-imput, 
# or imput-imput
diff_type <- function(p_df, contrasts) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  contrasts %>%
    pmap(function(a,b) {
      p_df %>%
        filter_data_group(c(a, b)) %>%
        inner_join(attr(p_df, 'annotation'), by = 'label') %>%
        select(id, group) %>%
        distinct %>%
        group_by(id) %>%
        summarise(diff_type = case_when(all(c(a, b) %in% group) ~ 'value_value',
                                        a %in% group ~ 'value_imput',
                                        b %in% group ~ 'imput_value',
                                        T ~ NA_character_),
        # summarise(diff_type = case_when(n() == 1 ~ 'value_imput',
        #                                 n() == 2 ~ 'value_value',
        #                                 T ~ NA_character_),
                  contrast = str_glue('{make.names(a)} - {make.names(b)}'),
                  .groups = 'drop')
      
        
    }) %>%
    bind_rows %>%
    pivot_wider(names_from = contrast, values_from = diff_type, 
                values_fill = 'imput_imput') %>%
    pivot_longer(-id, names_to = 'contrast', values_to = 'diff_type') %>%
    select(contrast, id, diff_type)
}
