#### Script to format the raw data for the analysis
# load dependencies
source("./scripts/load_deps.r")


#### INPUT
## read in the count matrices and fix naming order
raw_mat <-
  read_tsv(
    paste0(raw_data_folder, "/rsem_matrices/rsem.merged.gene_counts_deconvo.tsv")
  ) %>%
  # set NA's as 0's by default and enforce integers (for compact data format)
  mutate(across(where(is.numeric), ~ as.integer(replace_na(.x, 0)))) %>%
  # remove the corresponding transcript ID's and the Undetermined sample
  select(-"transcript_id(s)") %>%
  # fix the sample names
  rename_with(
    ~ .x %>%
      str_replace("HU10", "HUT10") %>%
      str_remove("(?<=\\HUT)[0]+(?=\\d)")
  ) %>%
  # fix missing sample names from the added samples
  rename_with(~ str_c(.x, "_UF"), matches("HUT9|HUT25|HUT29|HUT32|HUT27"))


## read in the metadata and consolidate samplenames
pheno <-
  read_tsv(paste0(raw_data_folder, "/mens_phase.tsv")) %>%
  # apparently some samples in this table are not present in the expression data
  mutate(
    # add a group column based on eithe UF or biopsy
    group = ifelse(str_detect(samplename, "UF"), "UF", "biopsy"),
    # rename groups to short group standards in this analysis
    cyclephase = case_when(
      cyclephase == "Proliferative" ~ "pro",
      cyclephase == "LH2_3" ~ "pre",
      cyclephase == "LH7_8" ~ "rec",
      cyclephase == "LH11_13" ~ "post",
      T ~ cyclephase
    )
  ) %>%
  { .[match(colnames(raw_mat)[-1], .$samplename), ] }


#### FILTER
removals <- c(
  "HUT10_UF", "HUT1_UF", "HUT1_biopsy", "HUT35_UF",
  "HUT42_UF", "HUT71_biopsy_2", "HUT71_biopsy_3"
)

# switch and filter some samples
raw_mat <-
  raw_mat %>%
  rename(
    ## this takes too much time to figure out how to switch programmatically
    ## just switch manually
    HUT23_biopsy = HUT23_UF,
    HUT23_UF = HUT23_biopsy
  ) %>%
  select(-removals[removals %in% colnames(raw_mat)])

pheno <-
  pheno %>%
  # but remove an samples actually not in the raw_mat (supposedly 8 of them)
  # and also quarantee that the samplenames are in the same order
  # as raw_mat colnames (needed for DESeq2)
  { .[match(colnames(raw_mat)[-1], .$samplename), ] }

## join count matrix with pheno for visualization or downstream analysis
expr_mat <- raw_mat %>% transform_count_mat(pheno)


#### OUTPUT
## Emit formatted data files
data_subfolder <- paste0(data_folder, "/filtered/")
raw_mat %>% write_feather(paste0(data_subfolder, "counts_raw.feather"))
expr_mat %>% write_feather(paste0(data_subfolder, "counts_format.feather"))
pheno %>% write_tsv(paste0(data_subfolder, "phenotype.tsv"))

## also write out raw file with switched and formatted annotation
raw_mat %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  rename(gene_id = external_gene_name) %>%
  write_feather(paste0(data_subfolder, "annot_raw.feather"))
