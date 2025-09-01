

library(tidyverse)

# Following code from https://www.uniprot.org/api-documentation/uniparc 
require(httr2)

base_url <- "https://ted.cathdb.info/api/v1/uniprot/summary/P0A6Y8?skip=0&limit=100"
base_url <- "https://www.ebi.ac.uk/interpro/api/protein/reviewed/P99999?extra_features=extra_features"
base_url <- "https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/reviewed/P0A6Y8/"


  #paste0("https://rest.uniprot.org/uniparc/", 
   #                uniparc_id, "/light")
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

data_extracted <- data_json$data %>% 
  tibble::enframe() %>% 
  tidyr::unnest_wider(value) 

data_extracted <- data_json$metadata %>% 
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

