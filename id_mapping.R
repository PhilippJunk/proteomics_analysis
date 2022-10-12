library(tidyverse)

#' Collapse UNIPROT IDs
#' 
#' This function takes a vector of UNIPROT IDs as produced by Maxquant
#' (in the column "Majority Protein IDs"), collapses isoforms. If after that,
#' multiple IDs are found, there are different modes to select only one. 
#' 
#' @param uniprot_ids A vector of UNIPROT IDs
#' @param multiple_mode Determines how to deal with entries that can not be
#'   broken down into a single entry. Possible values are "none", "manual",
#'   and "first" (default is "first").
#' @param manual_file In multiple_mode "manual", the file used for manual
#'   replacements. Required data frame needs to be saved in in a tab separated
#'   file with the columns "input" and "output".
#' @param separator_in The separator between multiple UNIPROT IDs and isoforms
#'   in an input (default ";").
#' @param separator_out In multiple_mode "none", the separator between
#'   multiple UNIPROT IDs when writing output (default ";").
#' @return Vector of collapes UNIPROT IDs
collapse_uniprot_ids <- function(
  uniprot_ids, 
  multiple_mode = 'first',
  manual_file, 
  separator_in = ';',
  separator_out = ';') {
  
  # check inputs
  if (!multiple_mode %in% c('none', 'manual', 'first')) {
    stop('Given mode not valid: ', multiple_mode)
  }
  
  if (multiple_mode == 'manual') {
    # import data set
    manual <- read.table(manual_file, header=T, sep='\t')
  }
  
  # only iterate over unique elements
  uniprot_ids_unique <- unique(uniprot_ids)
  out_unique <- vector('character', length(uniprot_ids_unique))
  
  for (i in 1:length(uniprot_ids_unique)) {
    id <- uniprot_ids_unique[i]
    out_id <- c()
    # split by input separator
    for (split_id in unlist(str_split(id, separator_in))) {
      # if starts with CON__, skip
      if (str_detect(split_id, 'CON__')) {
        next()
      }
      #if starts with REV__, skip
      else if (str_detect(split_id, 'REV__')) {
        next()
      }
      # if dash found, extract everything in front of it
      else if (str_detect(split_id, '-\\d+')) {
        out_id <- c(out_id, str_extract(split_id, '^[^-]+(?=-\\d+)'))
      } 
      # if none of those, extract everything
      else {
        out_id <- c(out_id, split_id)
      }
    }
    # take unique ids
    out_id <- unique(out_id)
    
    # if more than one UNIPROT id was found, proceed
    # depending on mode
    if (multiple_mode == 'none') {
      out_id <- str_c(out_id, collapse = separator_out)
    } else if (multiple_mode == 'first') {
      out_id <- out_id[1]
    } else if (multiple_mode == 'manual') {
      out_id_concat <- str_c(out_id, collapse = separator_out)
      if (length(out_id) == 1) {
        # pass
        invisible()
      } else if (out_id_concat %in% manual$input) {
        out_id <- manual %>% 
          filter(input == out_id_concat) %>%
          pull(output)
      } else {
        warning(out_id, ' not found in manual replacement. ',
                'Selecting first element.')
        out_id <- out_id[1]
      }
    }
    out_unique[i] <- out_id
  }
  
  # apply unique mapping to input
  left_join(data.frame(input = uniprot_ids),
            data.frame(input = uniprot_ids_unique,
                       output = out_unique), by = 'input') %>%
    pull(output)
}


#' Map UNIPROT IDs to HGNC gene names
#' 
#' Given a vector of UNIPROT IDs , extract corresponding HGNC gene 
#' names where applicable.
#' 
#' @param uniprot_ids A vector of UNIPROT IDs.
#' @param hgnc_database The database file of the HGNC genenames. Expected to
#'   be tab-separated, and contain the columns "hgnc_id", "symbol" and 
#'   "uniprot_ids".
#' @param multiple_mode Determines how cases are treated where one UNIPROT ID
#'   corresponds to multiple HGNC gene names. Possible values are "first" and
#'   "manual" and "none" (default "first").
#' @param manual_file In multiple_mode "manual", the file used for manual
#'   replacements. Required data frame needs to be saved in in a tab separated
#'   file with the columns "input" and "output".
#' @param separator_in The separator used to split multiple UNIPROT ID entries
#'   per input. Default is ";".
#' @param separator_out The separator used for separating multiple HGNC entries
#'   associated with a single UNIPROT ID, if multiple_mode is "none". The 
#'   default is ";".
#' @return Vector of HGNC gene names
map_uniprot_hgnc <- function(
  uniprot_ids,
  hgnc_database,
  multiple_mode = 'first',
  manual_file,
  separator_in = ';',
  separator_out = ';') {
  
  # check inputs
  if (!multiple_mode %in% c('none', 'manual', 'first')) {
    stop('Given mode not valid: ', multiple_mode)
  }
  
  if (multiple_mode == 'manual') {
    # import data set
    manual <- read.table(manual_file, header=T, sep='\t')
  }
  
  # For element in uniprot_id, find assignment in HGNC
  # There can be multiple HGNC IDs associated with one UNIPROT ID
  # There can also be multiple UNIPROT IDs associated with one HGNC ID
  
  # load hgnc database
  hgnc_db <- read_tsv(hgnc_database,
                      show_col_types = F) %>%
    select(hgnc_id, symbol, uniprot_ids) %>%
    separate_rows(uniprot_ids) %>%
    filter(!is.na(uniprot_ids))
  # TODO still produces parsing problems... annoying but well
  
  # make input unique
  uniprot_ids_unique <- unique(uniprot_ids)
  out_unique <- vector('character', length(uniprot_ids_unique))
  
  # perform mapping
  for (i in 1:length(uniprot_ids_unique)) {
    id <- uniprot_ids_unique[i]
    hits <- c()
    # split by semi-colon
    for (split_id in unlist(str_split(id, separator_in))) {
      hits <- c(hits, hgnc_db %>%
                  filter(uniprot_ids == split_id) %>%
                  pull(symbol))
    } 
    # take unique ids
    hits <- unique(hits)
    
    # if more than one UNIPROT id was found, proceed
    # depending on mode
    if (multiple_mode == 'none') {
      hits <- str_c(hits, collapse = separator_out)
    } else if (multiple_mode == 'first') {
      hits <- hits[1]
    } else if (multiple_mode == 'manual') {
      hits_concat <- str_c(hits, collapse = separator_out)
      if (length(hits) == 1) {
        # pass
        invisible()
      } else if (hits_concat %in% manual$input) {
        hits <- manual %>% 
          filter(input == hits_concat) %>%
          pull(output)
      } else {
        warning(out_id, ' not found in manual replacement. ',
                'Selecting first element.')
        hits <- out_id[1]
      }
    }
    out_unique[i] <- hits
  }
  
  # apply unique mapping to input
  left_join(data.frame(input = uniprot_ids),
            data.frame(input = uniprot_ids_unique,
                       output = out_unique), by = 'input') %>%
    pull(output)
}


#' Map HGNC gene names to SYSGO gene names
#' 
#' Given a vector of single HGNC IDs, extract corresponding SYSGO gene 
#' names where applicable.
#' 
#' @param hgnc_ids A vector of UNIPROT IDs.
#' @param sysgo_database The database file of the HGNC genenames. Expected to
#'   be tab-separated, and contain the columns "hgnc_id", "symbol" and 
#'   "uniprot_ids".
#' @return Vector of HGNC gene names
map_hgnc_sysgo <- function(
  hgnc_ids,
  sysgo_database) {
  # load SysGO mapping
  sysgo_mapping <- xlsx_cells(sysgo_database) %>% 
    filter(sheet == 'new symbols and aliases') %>%
    filter(row >= 2) %>%
    filter(data_type %in% c('character', 'numeric')) %>%
    mutate(character = case_when(
      data_type == 'character' ~ character,
      data_type == 'numeric' ~ as.character(numeric))) %>%
    select(character, row, col) %>%
    rename(name = character) %>%
    as.data.frame
  
  # make input unique
  hgnc_ids_unique <- unique(hgnc_ids)
  out_unique <- vector(mode='character', length=length(hgnc_ids_unique))
  
  # perform mapping
  for (i in 1:length(hgnc_ids_unique)) {
    id <- hgnc_ids_unique[i]
    match <- sysgo_mapping %>%
      filter(name == id)
    
    # check match
    # in order to map a match back to a unique name, we have to
    # extract the row of the match (match_row) and then extract
    # the first column for that row
    if (nrow(match) == 0) { # no match
      warning('Found no match for HGNC gene name ',
              id, ' in SysGO mapping.', sep='')
      out_unique[i] <- ''
    } else if (nrow(match) > 1) { # multiple matches
      # check if one of them has a hit in column 1
      if (nrow(match %>% filter(col == 1)) == 1) {
        match_row <- match %>%
          filter(col == 1) %>%
          pull(row)
        out_unique[i] <- sysgo_mapping %>%
          filter(row == match_row & col == 1) %>%
          pull(name)
      } else { # no hit in column 1, choose first entry
        match_row <- match %>% pull(row) %>% head(1)
        out_unique[i] <- sysgo_mapping %>%
          filter(row == match_row & col == 1) %>%
          pull(name)
        warning('Found multiple matches for HGNC gene name ', 
                id, ', used first one: ', out_unique[i], 
                sep='')
      }
    } else { # one match
      match_row <- match %>% pull(row)
      out_unique[i] <- sysgo_mapping %>%
        filter(row == match_row & col == 1) %>%
        pull(name)
    }
  }
  # apply unique mapping to input
  left_join(data.frame(input = hgnc_ids),
            data.frame(input = hgnc_ids_unique,
                       output = out_unique), by = 'input') %>%
    pull(output)
  
}

