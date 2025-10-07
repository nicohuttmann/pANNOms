

library(tidyverse)
library(pOmics3)
library(pANNOms)



fasta <- Biostrings::readAAStringSet("dev/uniprotkb_ecoli_AND_model_organism_8333_2024_06_06.fasta")

proteins <- names(fasta) %>% 
  str_extract("(?<=^..\\|).+?(?=\\|)")

taxons <- names(fasta) %>% 
  str_extract("(?<=OX=).+?(?= )") %>% 
  unique() %>% 
  as.numeric()


data_UniProt <- get_UniProt_data_o(accession = proteins, 
                                   fields = c("uniparc_id", 
                                               "sequence"), 
                                    taxon_id = taxons)

data_UniParc <- get_InterPro_data_from_UniProt(data_UniProt$Entry, 
                                               taxon_ids = taxons)

saveRDS(data_UniParc, "dev/InterPro.Rds")

library(tidyverse)
library(pOmics3)

uniprot_id <- list_EC50$Alanine$data_peptides$Protein.Group %>% 
  unique()


UniProt_fields()


uniparc_id <- get_UniProt_data_1o(str_extract(uniprot_id, 
                                              ".+?(?=;|$)"), 
                                  fields = "uniparc_id", 
                                  taxon_id = 83333)$UniParc

data_up <- get_UniProt_data("P0A6Y8", fields = "sequence")

# Check sequence
data$sequence$value == data_up$Sequence

accession <- a$data_peptides$Protein.Group %>% 
  unique() %>% 
  head()



get_InterPro_data_from_UniProt("P00760", split.position = T)


a <- readRDS(choose.files())
b <- readRDS(choose.files())




t <- get_InterPro_data_from_UniProt(a$data_peptides$Protein.Group %>% 
                                      unique() %>% 
                                      head(), split.position = T)




t <- get_InterPro_data_from_UniProt(list_EC50$Alanine$data_peptides$Protein.Group %>% 
                                      unique() %>% 
                                      head(), 
                                    export_as_UniProt = T)

test <- list_EC50$Alanine$data_peptides %>% 
  #filter(row_number() < 1000) %>% 
  left_join(get_InterPro_data_from_UniProt(.$Protein.Group, 
                                           export_as_UniProt = T, 
                                           taxon_id = 83333), 
            by = c("Protein.Group" = "Entry"))

test <- get_InterPro_data_from_UniProt(list_EC50$Alanine$data_peptides$Protein.Group[1:1000], 
                                       export_as_UniProt = T)


list_EC50$Alanine$data_peptides %>% 
  pull(Protein.Group) %>% 
  unique() %>% 
  get_InterPro_data_from_UniProt(export_as_UniProt = T)




extract_UniProt_seqinfo(b$data_protein_annotation, "Binding site")


#' Title
#'
#' @param accession 
#' @param data_UniProt 
#' @param taxon_id 
#' @param max.query 
#' @param keep.empty 
#' @param separate.multiple.sites 
#' @param split.position 
#' @param silent 
#'
#' @returns
#' @export
#'
#' @examples
get_InterPro_data_from_UniProt <- function(accession, 
                                           data_UniProt = NULL, 
                                           taxon_id = NULL, 
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
    if (!is.null(taxon_id)) {
      data_UniProt <- accession_query %>% 
        str_extract("^.+?(?=;|$)") %>% 
        get_UniProt_data_1o(fields = c("uniparc_id", 
                                       "sequence"), 
                            taxon_id = taxon_id)
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



tictoc::tic()
r <- map(cc(uniparc_id[1:20]), download_UniParc_sequence_features)
tictoc::toc()

download_UniParc_sequence_features("UPI000004C26F")


#' Title
#'
#' @param uniparc_id 
#' @param keep.empty 
#' @param separate.multiple.sites 
#' @param split.position 
#' @param return.sequence 
#' @param silent 
#'
#' @returns
#' @export
#'
#' @examples
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

