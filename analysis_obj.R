#' The OOP framework for storing data in this analysis pipeline.
#' 
#' Contains constructors for the `proteomics_data` class.

library(tidyverse)
library(forcats)

###############################################################################
# class definitions, as suggested as best practice in "Advanced R"

# low level constructor
new_proteomics_data <- function(df, annotation, has_tech_repl) {
  stopifnot(is.data.frame(df))
  stopifnot(is.data.frame(annotation))
  stopifnot(is.logical(has_tech_repl))
  stopifnot(is.character(df$id), is.character(df$label), 
            is.character(annotation$label), is.character(annotation$group))
  stopifnot(is.numeric(df$LFQ), is.numeric(annotation$biol_repl))
  if (has_tech_repl) {
    stopifnot(is.numeric(annotation$tech_repl))
  }

  proteomics_data <- structure(
    df,
    annotation = annotation,
    class = c('proteomics_data', 'tbl_df', 'tbl', 'data.frame'),
    has_tech_repl = has_tech_repl)

  return(proteomics_data)
}

# validator
validate_proteomics_data <- function(proteomics_data) {
  df <- as_tibble(proteomics_data)
  annotation <- attr(proteomics_data, 'annotation')
  has_tech_repl <- attr(proteomics_data, 'has_tech_repl')

  if (!inherits(proteomics_data, 'proteomics_data')) {
    stop('Please provide an object of class `proteomics_data` as input. ', 
         'For more information, see `?proteomics_data`.',
         call. = FALSE)
  }
  
  if (!is.data.frame(annotation)) {
    stop(
      '`annotation` must be a `data.frame`.',
      call. = FALSE
    )
  }

  if (!is.logical(has_tech_repl)) {
    stop(
      '`has_tech_repl` must be a logical vector.',
      call. = FALSE
    )
  }
  
  if (!all(colnames(df) == c('id', 'label', 'LFQ'))) {
    stop(
      'Data must consists of the columns `id`, `label` and `LFQ`.',
      call. = FALSE
    )
  }
  
  if (has_tech_repl) {
    if (!all(colnames(annotation) == c('label', 'group', 'biol_repl', 'tech_repl'))) {
      stop('Annotation is expected to have the columns `label`, `group`, `biol_repl` and `tech_repl`.', 
           call. = FALSE)
    }
  } else {
    if (!all(colnames(annotation) == c('label', 'group', 'biol_repl'))) {
      stop('Annotation is expected to have the columns `label`, `group` and `biol_repl`.', 
           call. = FALSE)
    }
  }
  
  if (any(is.na(df$LFQ))) {
    stop('LFQ intensities must be non-missing.',
         call. = FALSE)
    
  }
      
  if (any(df$LFQ < 0)) {
    stop('LFQ intensities must be non-missing and greater than zero.',
         call. = FALSE)
  }
  
  if (!all(unique(df$label) %in% annotation$label)) {
    stop('Not all `label`s in data are present in `annotation`.', 
         call. = FALSE)
  }
  
  if (!all(annotation$label %in% unique(df$label))) {
    stop('Not all `label`s in `annotation` are present in data', 
         call. = FALSE)
  }
  
  if (any(duplicated(annotation$label))) {
    stop('`label`s in `annotation` must be unique.', 
         call. = FALSE)
  }
  
  return(proteomics_data)
}

# user-friendly helper
proteomics_data <- function(
    df, annotation, has_tech_repl, is_log2,
    df_id = 'id', df_label = 'label', df_LFQ = 'LFQ', 
    annotation_label = 'label', annotation_group = 'group', 
    annotation_biol_repl = 'biol_repl', annotation_tech_repl ='tech_repl') {
  
  df <- df %>% 
    select({df_id}, {df_label}, {df_LFQ}) %>%
    rename(id = {df_id},
           label = {df_label},
           LFQ = {df_LFQ}) %>%
    select(id, label, LFQ) %>%
    mutate(id = make.names(id),
           label = make.names(label))
  if (!is_log2) {
    stopifnot(all(df$LFQ > 2))
    df <- df %>%
      mutate(LFQ = log2(LFQ))
  }
  if (has_tech_repl) {
    annotation <- annotation %>%
      select({annotation_label}, {annotation_group}, 
             {annotation_biol_repl}, {annotation_tech_repl}) %>%
      rename(label = {annotation_label},
             group = {annotation_group},
             biol_repl = {annotation_biol_repl},
             tech_repl = {annotation_tech_repl}) %>%
      select(label, group, biol_repl, tech_repl)
  } else {
    annotation <- annotation %>%
      select({annotation_label}, {annotation_group}, 
             {annotation_biol_repl}) %>%
      rename(label = {annotation_label},
             group = {annotation_group},
             biol_repl = {annotation_biol_repl}) %>%
      select(label, group, biol_repl)
  }
  annotation %>% 
    mutate(group = make.names(group),
           label = make.names(label))
  
  validate_proteomics_data(new_proteomics_data(
    df, annotation,has_tech_repl))
}

###############################################################################
# helper functions

# wrapper around inner join with annotation attr
join_annotation <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  p_df %>%
    inner_join(attr(p_df, 'annotation'), by = 'label')
}

# safely collapsing technical replicates
collapse_tech_repl <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  if (!attr(p_df, 'has_tech_repl')) {
    stop('The `proteomics_data` object must have technical replicates in order to collapse them.',
         call. = FALSE)
  }
  
  annotation <- attr(p_df, 'annotation') %>%
    mutate(new_label = str_c(group, biol_repl, sep='_'))
  
  p_df <- p_df %>%
    # transform back from log2 for collapsing
    mutate(LFQ = 2^LFQ) %>%
    inner_join(annotation, by = 'label') %>%
    group_by(id, new_label) %>%
    summarise(LFQ = median(LFQ), .groups='drop')
  
  annotation <- annotation %>%
    select(new_label, group, biol_repl) %>%
    distinct
    
  proteomics_data(p_df, annotation, has_tech_repl = FALSE, is_log2 = FALSE,
                  df_label = 'new_label', annotation_label = 'new_label')
} 

# filter data by group and id
# df_filter expected to have group and id column!
filter_data_group_id <- function(p_df, df_filter) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  df_filter <- df_filter %>%
    select(id, group)

  # filter data
  p_df <- p_df %>%
    inner_join(attr(p_df, 'annotation'), by='label') %>%
    inner_join(df_filter, by = c('id', 'group'))
  # filter annotation
  annotation <- attr(p_df, 'annotation') %>%
    inner_join(df_filter %>% select(group) %>% distinct, by = 'group') 
  
  proteomics_data(
    p_df, annotation, 
    has_tech_repl = attr(p_df, 'has_tech_repl'), is_log2 = TRUE)
}

# filter data by groups
filter_data_group <- function(p_df, groups) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)

  # filter data
  p_df <- p_df %>%
    inner_join(attr(p_df, 'annotation'), by='label') %>%
    filter(group %in% groups)
  # filter annotation
  annotation <- attr(p_df, 'annotation') %>%
    filter(group %in% groups)
  
  proteomics_data(
    p_df, annotation, 
    has_tech_repl = attr(p_df, 'has_tech_repl'), is_log2 = TRUE)
}

# remove samples from data
remove_samples <- function(p_df, samples) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  # filter data
  p_df <- p_df %>%
    filter(!label %in% samples)
  # filter annotation
  annotation <- attr(p_df, 'annotation') %>%
    filter(!label %in% samples)
  
  proteomics_data(
    p_df, annotation,
    has_tech_repl = attr(p_df, 'has_tech_repl'), is_log2 = TRUE)
}

deconstruct_groups <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  annotation <- attr(p_df, 'annotation')

  sets <- p_df %>%
    inner_join(annotation, by='label') %>%
    mutate(group = factor(group) %>% fct_infreq %>% fct_rev) %>%
    select(id, group) %>%
    distinct %>%
    arrange(group) %>%
    group_by(id) %>%
    summarise(set = str_c(group, collapse = '__'),
              .groups = 'drop')
  
  sets %>%
    count(set) %>% 
    arrange(-n) %>%
    pull(set) %>%
    set_names %>%
    map(function(s) {
      sets %>%
        filter(set == s) %>%
        pull(id)
    })
}
  

# TODO function that remove label or group from df and annotation
# TODO? function that collapses biological replicates