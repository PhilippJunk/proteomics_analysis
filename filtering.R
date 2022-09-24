library(tidyverse)

###############################################################################
#' The available filtering methods.
#'
#' @return A character vector of available filtering methods.
filtering_methods <- function() {
  c('none', 'whole_dataset', 'any_group', 'each_group')
}


#' Filtering expression data set with different methods
#'
#' @param p_df The data to run filtering on. Of class `proteomics_data`.
#' @param method The method used for filtering. See `filtering_methods()` for
#'   available option.
#' @param threshold The threshold how many samples need to be present in the
#'   subset the different methods look at to retain a data point.
#' @return A filtered object of class `proteomics_data`.
filtering <- function(
    p_df,
    method,
    threshold) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  method <- match.arg(method, choices=filtering_methods())
  if (method %in% c('whole_dataset', 'any_group', 'each_group')) {
    if (!all(threshold >=0, threshold <= 1, is.numeric(threshold), length(threshold) == 1)) {
      stop('Threshold must be a numeric value between zero and one.',
           call. = FALSE)
    }
  }
  
  # extract data
  data <- as_tibble(p_df)
  annotation <- attr(p_df, 'annotation')
  
  # perform computation
  if (method == 'none') {
    data_filtered <- data
  } else if (method == 'whole_dataset') {
    data_filtered <- .filtering_whole_dataset(data, threshold)
  } else if (method == 'any_group') {
    data_filtered <- .filtering_any_group(data, threshold, annotation)
  } else if (method == 'each_group') {
    data_filtered <- .filtering_each_group(data, threshold, annotation)
  }
  
  # check if all groups are still present between data and annotation
  if (!all(annotation$label %in% unique(data_filtered$label))) {
    missing_labels <- annotation$label[!annotation$label %in% unique(data_filtered$label)]
    warning('Filtering removed the following labels completely from the data:',
            str_glue(' {missing_labels}'))
    annotation <- annotation %>%
      filter(!label %in% missing_labels)
  }
  
  # reform to object
  proteomics_data(
    data_filtered, annotation, 
    has_tech_repl = attr(p_df, 'has_tech_repl'),
    is_log2 = TRUE)
}


#' Filtering data set by presence of values across whole data set
#'
#' @param data The data frame to run filtering on.
#' @param threshold The threshold how many samples need to be present in the
#'   subset the different methods look at to retain a data point.
#' @return A filtered data frame.
.filtering_whole_dataset <- function(data, threshold) {
  n_samples <- data %>% pull(label) %>% unique %>% length
  data %>%
    group_by(id) %>%
    mutate(present_samples = n()) %>%
    ungroup %>%
    filter(present_samples/n_samples >= threshold)
}


#' Filtering data set by presence of values in at least one group in the data
#' set
#'
#' @param data The data frame to run filtering on.
#' @param threshold The threshold how many samples need to be present in the
#'   subset the different methods look at to retain a data point.
#' @param annotation A data frame that contains information about which samples
#'   in the input data frame belong to which condition.
#' @return A filtered data frame.
.filtering_any_group <- function(data, threshold, annotation) {
  data %>%
    # join with annotation
    inner_join(annotation %>%
                 group_by(group) %>%
                 mutate(n_group = n()) %>%
                 ungroup,
               by='label') %>%
    # calculate ratio of present sample for each peptide for each group
    group_by(id, group) %>%
    mutate(present_ratio = n()/n_group) %>%
    ungroup %>%
    # get best ratio for each peptide
    group_by(id) %>%
    mutate(best_present_ratio = max(present_ratio)) %>%
    ungroup %>%
    # filter by best ratio
    filter(best_present_ratio >= threshold)
}


#' Filtering data set by presence of values for each group in the data
#'
#' @param data The data frame to run filtering on.
#' @param threshold The threshold how many samples need to be present in the
#'   subset the different methods look at to retain a data point.
#' @param annotation A data frame that contains information about which samples
#'   in the input data frame belong to which condition.
#' @return A filtered data frame.
.filtering_each_group <- function(data, threshold, annotation) {
  df <- data %>%
    # join with annotation
    inner_join(annotation %>%
                 group_by(group) %>%
                 mutate(n_group = n()) %>%
                 ungroup,
               by='label') %>%
    # calculate ratio of present sample for each peptide for each group
    group_by(id, group) %>%
    mutate(present_ratio = n()/n_group) %>%
    ungroup %>%
    # filter by ratio
    filter(present_ratio >= threshold)
}

