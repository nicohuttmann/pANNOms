#' Title
#'
#' @param accession 
#' @param data_UniProt (optional) predownloaded UniProt data with 
#' download_UniProt_data() or download_UniProt_data_1o()
#' @param taxon_ids 
#' @param max.query maximum
#' @param keep.empty Keep Ids without data?
#' @param separate.multiple.sites split disconnected features into multiple rows
#' @param split.position add columns 'from' and 'to' 
#' @param silent Suppress messages?
#'
#' @returns
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
get_InterPro_data_from_UniProt <- function(accession, 
                                           data_UniProt = NULL, 
                                           taxon_ids = NULL, 
                                           max.query = 100, 
                                           keep.empty = F, 
                                           separate.multiple.sites = T, 
                                           split.position = T, 
                                           export_as_UniProt = F, 
                                           silent = T) {
  
  # Only query unique accessions 
  accession_query <- unique(accession)
  
  # Get UniProt annotations 
  if (is.null(data_UniProt)) {
    if (!is.null(taxon_ids)) {
      data_UniProt <- accession_query %>% 
        str_extract("^.+?(?=;|$)") %>% 
        get_UniProt_data_o(fields = c("uniparc_id", 
                                       "sequence"), 
                            taxon_ids = taxon_ids)
    } else {
      data_UniProt <- accession_query %>% 
        str_extract("^.+?(?=;|$)") %>% 
        get_UniProt_data(fields = c("uniparc_id", 
                                    "sequence"), 
                         max.query = max.query)
    }
  }
  
  # Download UniParc data individually 
  data_UniParc <- data_UniProt %>% 
    dplyr::filter(!is.na(UniParc)) %>% 
    dplyr::pull(UniParc) %>% 
    setNames(., .) %>% 
    purrr::map(
      \(x) download_UniParc_sequence_features(
        uniparc_id = x, 
        keep.empty = keep.empty, 
        separate.multiple.sites = separate.multiple.sites, 
        split.position = split.position, 
        return.sequence = T, 
        silent = silent), 
      .progress = "Donwloading InterPro data from UniParc.")
  
  # Check if sequences match 
  data_sequence_check <- data_UniProt %>% 
    left_join(data_UniParc %>% 
                map(\(x) tibble::tibble(sequence_UniParc = x[["sequence"]])) %>% 
                bind_rows(.id = "UniParc"), 
              by = "UniParc") %>% 
    # Remove UniParc entries without sequence and without filters 
    dplyr::filter(!is.na(sequence_UniParc)) %>% 
    dplyr::mutate(sequence_match = Sequence == sequence_UniParc)
  
  # In theory, this should not happen :)
  if (any(!data_sequence_check$sequence_match)) {
    warning("Some sequences do not match between UniProt and UniParc.")
  }
  
  # Combine UniProt Ids with UniParc features 
  data_output <- data_UniProt %>% 
    left_join(data_UniParc %>% 
                map("data_features") %>% 
                bind_rows(.id = "UniParc"), 
              by = "UniParc") %>% 
    dplyr::select(-c(Sequence)) %>% 
    mutate(InterPro = paste0('INTERPRO ', 
                             stringr::str_replace(location, '-', '..'), 
                             '; /name="', name, 
                             '"; /interpro_id="', InterProId, 
                             '"; /evidence="', 
                             UniParc, '|', 
                             database, ':', database, ':', databaseId, '"'), 
           .after = 'Entry')
  
  # If export as UniProt format 
  if (export_as_UniProt) {
    
    if (!separate.multiple.sites) 
      stop("separate.multiple.sites must be set TRUE.")
    
    data_output <- data_output %>% 
      dplyr::summarise(InterPro = paste(InterPro, collapse = "; "), 
                       .by = "Entry")
    
  }
  
  
  return(data_output)
  
}


#' Download sequence features data for UniParc sequences 
#'
#' @param uniparc_id UniParc Id (one at a time)
#' @param keep.empty Keep Ids without data?
#' @param separate.multiple.sites split disconnected features into multiple rows
#' @param split.position add columns 'from' and 'to' 
#' @param return.sequence keep the protein sequence column in the output 
#' @param silent Suppress messages? 
#'
#' @returns
#' @export
#'
#' @examples
#' \dontrun{
#' download_UniParc_sequence_features("UPI000004C26F") 
#' }
download_UniParc_sequence_features <- function(uniparc_id, 
                                               keep.empty = F, 
                                               separate.multiple.sites = T, 
                                               split.position = T, 
                                               return.sequence = F, 
                                               silent = T) {
  
  if (!hasArg(uniparc_id)) 
    stop("Please provide a UniParc Id.")
  
  if (length(uniparc_id) != 1)
    stop("Please provide 1 Id at a time.")
  
  
  # Following code from https://www.uniprot.org/api-documentation/uniparc 
  require(httr2)
  
  base_url <- paste0("https://rest.uniprot.org/uniparc/", 
                     uniparc_id, "/light")
  params <- list(
    #fields = "upi,organism,length"
  )
  
  req <- httr2::request(base_url)
  req %>% 
    httr2::req_headers(
      accept = "text/plain;format=tsv")
  
  req %>% 
    httr2::req_url_query(!!!params)
  
  resp <- httr2::req_perform(req)
  
  if (httr2::resp_status(resp) != 200) {
    stop(sprintf("Error %d: %s", 
                 httr2::resp_status(resp), 
                 httr2::resp_body_string(resp)))
  }
  
  # Extract json data 
  data_json <- httr2::resp_body_json(resp)
  
  
  # Extract features 
  if (!"sequenceFeatures" %in% names(data_json)) {
    if (!silent) {
      warning(paste0("No features found. Returning empty data frame for ", 
                     uniparc_id, 
                     "."))
    }
    return(tibble::tibble())
  }
  
  data_extracted <- data_json$sequenceFeatures %>% 
    tibble::enframe() %>% 
    tidyr::unnest_wider(value) 
  
  if (!"interproGroup" %in% names(data_extracted)) 
    data_extracted <- data_extracted %>% 
    mutate(interproGroup = list(NULL))
  
  if (!"locations" %in% names(data_extracted)) 
    data_extracted <- data_extracted %>% 
    mutate(locations = list(NULL))
  
  data_extracted <- data_extracted %>% 
    dplyr::mutate(location = 
                    purrr::map_chr(locations, 
                                   \(x) purrr::map(x, \(y) paste(y[c("start", "end")], 
                                                                 collapse = "-")) %>% 
                                     unlist() %>% paste(collapse = ";"))) %>% 
    dplyr::mutate(name = purrr::map_chr(interproGroup, \(x) ifelse(length(x) > 0,
                                                                   x[["name"]],
                                                                   ""))) %>% 
    dplyr::mutate(InterProId = purrr::map_chr(interproGroup, \(x) ifelse(length(x) > 0,
                                                                         x[["id"]],
                                                                         ""))) %>% 
    dplyr::select(name, 
                  location, 
                  InterProId, 
                  databaseId, 
                  database) %>% 
    dplyr::arrange(as.numeric(str_extract(location, "\\d+")))
  
  # Remove entries 
  if (!keep.empty) {
    data_extracted <- data_extracted %>% 
      dplyr::filter(name != "")
  }
  
  # Separate entries 
  if (separate.multiple.sites) {
    data_extracted <- data_extracted %>% 
      separate_longer_delim(location, ";")
    # Split position 
    if (split.position) {
      data_extracted <- data_extracted %>% 
        dplyr::mutate(from = as.numeric(str_extract(location, "\\d+")), 
                      to = as.numeric(str_extract(location, "(?<=-)\\d+")), 
                      .before = "location") %>% 
        dplyr::arrange(from, to)
    }
  }
  
  # Add sequence 
  if (return.sequence) 
    data_extracted <- list(data_features = data_extracted, 
                           sequence = data_json$sequence$value)
  
  
  return(data_extracted)
  
}

