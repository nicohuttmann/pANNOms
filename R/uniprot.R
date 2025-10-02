#' Return list of vectors of UniProt data fields 
#'
#' @returns
#' @export
#'
#' @examples
#' UniProt_fields()
UniProt_fields <- function() {
  
  UniProt_fields <- list(
    `Names & Taxonomy` = c("Entry Name" = "id",
                           "Gene Names" = "gene_names",
                           "Gene Names (ordered locus)" = "gene_oln",
                           "Gene Names (ORF)" = "gene_orf",
                           "Gene Names (primary)" = "gene_primary",
                           "Gene Names (synonym)" = "gene_synonym",
                           "Organism" = "organism_name",
                           "Organism (ID)" = "organism_id",
                           "Protein names" = "protein_name",
                           "Proteomes" = "xref_proteomes",
                           "Taxonomic lineage" = "lineage",
                           "Taxonomic lineage (Ids)" = "lineage_ids",
                           "Virus hosts" = "virus_hosts"), 
    Sequences = c("Alternative products (isoforms)" = "cc_alternative_products",
                  "Alternative sequence" = "ft_var_seq",
                  "Erroneous gene model prediction" = "error_gmodel_pred",
                  "Fragment" = "fragment",
                  "Gene encoded by" = "organelle",
                  "Length" = "length",
                  "Mass" = "mass",
                  "Mass spectrometry" = "cc_mass_spectrometry",
                  "Natural variant" = "ft_variant",
                  "Non-adjacent residues" = "ft_non_cons",
                  "Non-standard residue" = "ft_non_std",
                  "Non-terminal residue" = "ft_non_ter",
                  "Polymorphism" = "cc_polymorphism",
                  "RNA Editing" = "cc_rna_editing",
                  "Sequence" = "sequence",
                  "Sequence caution" = "cc_sequence_caution",
                  "Sequence conflict" = "ft_conflict",
                  "Sequence uncertainty" = "ft_unsure",
                  "Sequence version" = "sequence_version"), 
    Function = c("Absorption" = "absorption",
                 "Active site" = "ft_act_site",
                 "Binding site" = "ft_binding",
                 "Catalytic activity" = "cc_catalytic_activity",
                 "Cofactor" = "cc_cofactor",
                 "DNA binding" = "ft_dna_bind",
                 "EC number" = "ec",
                 "Activity regulation" = "cc_activity_regulation",
                 "Function [CC]" = "cc_function",
                 "Kinetics" = "kinetics",
                 "Pathway" = "cc_pathway",
                 "pH dependence" = "ph_dependence",
                 "Redox potential" = "redox_potential",
                 "Rhea ID" = "rhea",
                 "Site" = "ft_site",
                 "Temperature dependence" = "temp_dependence"), 
    Miscellaneous = c("Annotation" = "annotation_score",
                      "Caution" = "cc_caution",
                      "Keywords" = "keyword",
                      "Keyword ID" = "keywordid",
                      "Miscellaneous [CC]" = "cc_miscellaneous",
                      "Protein existence" = "protein_existence",
                      "Reviewed" = "reviewed",
                      "Tools" = "tools",
                      "UniParc" = "uniparc_id",
                      "Comments" = "comment_count",
                      "Features" = "feature_count"), 
    Interaction = c("Interacts with" = "cc_interaction",
                    "Subunit structure" = "cc_subunit"), 
    Expression = c("Developmental stage" = "cc_developmental_stage",
                   "Induction" = "cc_induction",
                   "Tissue specificity" = "cc_tissue_specificity"), 
    `Gene Ontology (GO)` = c("Gene Ontology (biological process)" = "go_p",
                             "Gene Ontology (cellular component)" = "go_c",
                             "Gene Ontology (GO)" = "go",
                             "Gene Ontology (molecular function)" = "go_f",
                             "Gene Ontology IDs" = "go_id"), 
    `Pathology & Biotech` = c("Allergenic Properties" = "cc_allergen",
                              "Biotechnological use" = "cc_biotechnology",
                              "Disruption phenotype" = "cc_disruption_phenotype",
                              "Involvement in disease" = "cc_disease",
                              "Mutagenesis" = "ft_mutagen",
                              "Pharmaceutical use" = "cc_pharmaceutical",
                              "Toxic dose" = "cc_toxic_dose"), 
    `Subcellular location` = c("Intramembrane" = "ft_intramem",
                               "Subcellular location [CC]" = "cc_subcellular_location",
                               "Topological domain" = "ft_topo_dom",
                               "Transmembrane" = "ft_transmem"), 
    `PTM / Processing` = c("Chain" = "ft_chain",
                           "Cross-link" = "ft_crosslnk",
                           "Disulfide bond" = "ft_disulfid",
                           "Glycosylation" = "ft_carbohyd",
                           "Initiator methionine" = "ft_init_met",
                           "Lipidation" = "ft_lipid",
                           "Modified residue" = "ft_mod_res",
                           "Peptide" = "ft_peptide",
                           "Post-translational modification" = "cc_ptm",
                           "Propeptide" = "ft_propep",
                           "Signal peptide" = "ft_signal",
                           "Transit peptide" = "ft_transit"), 
    Structure = c("3D" = "structure_3d",
                  "Beta strand" = "ft_strand",
                  "Helix" = "ft_helix",
                  "Turn" = "ft_turn"), 
    Publications = c("PubMed ID" = "lit_pubmed_id",
                     "DOI ID" = "lit_doi_id"), 
    `Data of` = c("Date of creation" = "date_created",
                  "Date of last modification" = "date_modified",
                  "Date of last sequence modification" = "date_sequence_modified",
                  "Entry version" = "version"), 
    `Family & Domains` = c("Coiled coil" = "ft_coiled",
                           "Compositional bias" = "ft_compbias",
                           "Domain [CC]" = "cc_domain",
                           "Domain [FT]" = "ft_domain",
                           "Motif" = "ft_motif",
                           "Protein families" = "protein_families",
                           "Region" = "ft_region",
                           "Repeat" = "ft_repeat",
                           "Sequence similarities" = "cc_similarity",
                           "Zinc finger" = "ft_zn_fing"), 
    # External Resources 
    xref_Sequence = c("CCDS" = "xref_ccds",
                      "EMBL" = "xref_embl",
                      "GeneRIF" = "xref_generif",
                      "PIR" = "xref_pir",
                      "RefSeq" = "xref_refseq"), 
    `xref_3D structure` = c("AlphaFoldDB" = "xref_alphafolddb",
                            "BMRB" = "xref_bmrb",
                            "EMDB" = "xref_emdb",
                            "PCDDB" = "xref_pcddb",
                            "PDB" = "xref_pdb_full",
                            "PDBsum" = "xref_pdbsum",
                            "SASBDB" = "xref_sasbdb",
                            "SMR" = "xref_smr"), 
    `xref_Protein-protein interaction` = c("BioGRID" = "xref_biogrid_full",
                                           "CORUM" = "xref_corum",
                                           "ComplexPortal" = "xref_complexportal_full",
                                           "DIP" = "xref_dip",
                                           "ELM" = "xref_elm",
                                           "IntAct" = "xref_intact_full",
                                           "MINT" = "xref_mint",
                                           "STRING" = "xref_string"), 
    xref_Chemisty = NA_character_, 
    `xref_Protein family/group` = NA_character_, 
    xref_PTM = NA_character_, 
    `xref_Genetic variation` = NA_character_, 
    `xref_2D gel` = NA_character_, 
    xref_Proteomic = NA_character_, 
    `xref_Protocols and materials` = NA_character_, 
    `xref_Genome annotation` = c("Ensembl" = "xref_ensembl",
                                 "EnsemblBacteria" = "xref_ensemblbacteria",
                                 "EnsemblFungi" = "xref_ensemblfungi",
                                 "EnsemblMetazoa" = "xref_ensemblmetazoa",
                                 "EnsemblPlants" = "xref_ensemblplants",
                                 "EnsemblProtists" = "xref_ensemblprotists",
                                 "GeneID" = "xref_geneid",
                                 "Gramene" = "xref_gramene",
                                 "KEGG" = "xref_kegg",
                                 "MANE-Select" = "xref_mane-select",
                                 "PATRIC" = "xref_patric",
                                 "UCSC" = "xref_ucsc",
                                 "VectorBase" = "xref_vectorbase",
                                 "WBParaSite" = "xref_wbparasite"), 
    `xref_Organism-specific` = NA_character_, 
    xref_Phylogenomic = NA_character_, 
    `xref_Enzyme and pathway` = NA_character_, 
    xref_Miscellaneous = NA_character_, 
    `xref_Gene expression` = NA_character_, 
    `xref_Family and domain` = c("AntiFam" = "xref_antifam_full",
                                 "CDD" = "xref_cdd_full",
                                 "DisProt" = "xref_disprot",
                                 "FunFam" = "xref_funfam_full",
                                 "Gene3D" = "xref_gene3d_full",
                                 "HAMAP" = "xref_hamap_full",
                                 "IDEAL" = "xref_ideal",
                                 "InterPro" = "xref_interpro_full",
                                 "NCBIfam" = "xref_ncbifam_full",
                                 "PANTHER" = "xref_panther_full",
                                 "PIRSF" = "xref_pirsf_full",
                                 "PRINTS" = "xref_prints_full",
                                 "PROSITE" = "xref_prosite_full",
                                 "Pfam" = "xref_pfam_full",
                                 "SFLD" = "xref_sfld_full",
                                 "SMART" = "xref_smart_full",
                                 "SUPFAM" = "xref_supfam_full"))
  
  return(UniProt_fields)
  
}


#' Download UniProt data for given protein accessions and data fields 
#' (see available fields with UniProt_fields())
#'
#' @param accession vector of UniProt accessions 
#' @param fields UniProt data fields to query
#' @param max.query maximum number of accessions to query at once; if the 
#' the number exceeds max.query, the query is split up in multiple parts 
#'
#' @returns
#' @export
#'
#' @examples
get_UniProt_data <- function(accession, 
                             fields = c("accession", 
                                        "gene_names", 
                                        "organism_name"), 
                             max.query = 100) {
  
  # Add accession as field
  if (!"accession" %in% fields) fields <- c("accession", fields)
  
  # Check accession for ;
  if (any(stringr::str_detect(accession, ";"))) 
    stop("There are accessions contanining a semicolon (;);, please remove or correct protein groups with multiple Ids.")
  
  # Only query unique accessions 
  accession_query <- unique(accession)
  
  
  # Formulate query/ies and download data 
  if (length(accession_query) > max.query) {
    
    l <- length(accession_query)
    from <- seq(1, l, max.query)
    to <- c(seq(1, l, max.query)[-1] - 1, l)
    
    data_download <- purrr::map2(from, to, 
                                 \(from, to) get_UniProt_data(accession_query[from:to], 
                                                              fields = fields, 
                                                              max.query = max.query)) %>% 
      bind_rows()
    
  } else {
    
    query_url <- paste0("https://rest.uniprot.org/uniprotkb/stream?", 
                        "format=tsv", 
                        "&fields=", 
                        paste(fields, collapse = "%2C"), 
                        "&query=", 
                        paste(
                          paste0("accession%3A", accession_query), 
                          collapse = "+OR+"))
    
    data_download <- vroom::vroom(query_url, 
                                  delim = "\t", 
                                  col_types = readr::cols())
    
  }
  
  # Merge given accessions and downloaded data 
  data_output <- dplyr::left_join(tibble::tibble(Entry = accession), 
                                  data_download, 
                                  by = "Entry")
  
  # Return tibble with accessions as Entry and data columns 
  return(data_output)
  
}


#' Download UniProt data for given protein accessions, multiple taxonomy 
#' identifiers and data fields (faster for 1000s of proteins; see available fields with 
#' UniProt_fields())
#'
#' @param accession vector of UniProt accessions 
#' @param fields UniProt data fields to query
#' @param taxon_id taxonomy Id
#'
#' @returns
#' @export
#'
#' @examples
get_UniProt_data_1o <- function(accession, 
                                fields = c("accession", 
                                           "gene_names", 
                                           "organism_name"), 
                                taxon_id = c(human = 9606, 
                                             mouse = 10900, 
                                             E.coliK12 = 83333)) {
  
  # Check organism identifier
  if (length(taxon_id) > 1) 
    warning("More than one taxon_id provided, only using the first one.")
  
  # Add accession as field
  if (!"accession" %in% fields) fields <- c("accession", fields)
  
  
  # Formulate query and download data 
  query_url <- paste0("https://rest.uniprot.org/uniprotkb/stream?", 
                      "format=tsv", 
                      "&fields=", 
                      paste(fields, collapse = "%2C"), 
                      "&query=%28taxonomy_id%3A", 
                      taxon_id[1], "%29")
  
  data_download <- vroom::vroom(query_url, 
                                delim = "\t", 
                                col_types = readr::cols())
  
  
  # Merge given accessions and downloaded data 
  if (hasArg(accession)) {
    data_output <- dplyr::left_join(tibble::tibble(Entry = accession), 
                                    data_download, 
                                    by = "Entry")
  } else {
    data_output <- data_download
  }
  
  
  # Return tibble with accessions as Entry and data columns 
  return(data_output)
  
}


#' Download UniProt data for given protein accessions, taxonomy identifiers and 
#' data fields (faster for 1000s of proteins; see available fields with 
#' UniProt_fields())
#'
#' @param accession vector of UniProt accessions 
#' @param fields UniProt data fields to query
#' @param taxon_ids taxonomy Ids
#'
#' @returns
#' @export
#'
#' @examples
get_UniProt_data_o <- function(accession, 
                               fields = c("accession", 
                                          "gene_names", 
                                          "organism_name"), 
                               taxon_ids = c(human = 9606, 
                                            mouse = 10900, 
                                            E.coliK12 = 83333)) {
  
  # Check organism identifier
  # if (length(taxon_id) > 1) 
  #   warning("More than one taxon_id provided, only using the first one.")
  
  data_download <- purrr::map(
    taxon_ids, 
    \(x) {
      
      # Formulate query and download data 
      query_url <- paste0("https://rest.uniprot.org/uniprotkb/stream?", 
                          "format=tsv", 
                          "&fields=", 
                          paste(fields, collapse = "%2C"), 
                          "&query=%28model_organism%3A", 
                          x, "%29")
      
      vroom::vroom(query_url, 
                   delim = "\t", 
                   col_types = readr::cols())}) %>% 
    dplyr::bind_rows()
  
  
  # Merge given accessions and downloaded data 
  if (hasArg(accession)) {
    data_output <- dplyr::left_join(tibble::tibble(Entry = accession), 
                                    data_download, 
                                    by = "Entry")
  } else {
    data_output <- data_download
  }
  
  
  # Return tibble with accessions as Entry and data columns 
  return(data_output)
  
}
