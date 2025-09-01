



# From https://degronopedia.com/degronopedia/degron_motifs
db_degrons <- openxlsx::read.xlsx("https://degronopedia.com/degronopedia/download/data/DEGRONOPEDIA_degron_dataset.xlsx", 
                                  sheet = "Degrons")


data <- get_UniProt_data_1o(fields = "sequence")

data_deg <- data %>% 
  mutate(deg1 = str_detect(Sequence, 
                           db_degrons %>% 
                             filter(Degron == "RxxG") %>% 
                             pull(Degron_regex)), 
         n_deg1 = str_count(Sequence, 
                            db_degrons %>% 
                              filter(Degron == "RxxG") %>% 
                              pull(Degron_regex)), 
         motif_deg1 = str_extract(Sequence, 
                                  db_degrons %>% 
                                    filter(Degron == "RxxG") %>% 
                                    pull(Degron_regex)))




data_domains <- get_UniProt_data_1o(fields = "ft_domain")

