### Script to format the raw data for the analysis
# load dependencies
source("./scripts/load_deps.r")


#### INPUT
## read in the count matrices and fix naming order
raw_mat <-
  read_tsv(
    paste0(raw_data_folder, "/rsem.merged.gene_counts_deconvo.tsv")
  ) %>%
  # set NA's as 0's by default and enforce integers (for compact data format)
  mutate(across(where(is.numeric), ~ as.integer(replace_na(.x, 0)))) %>%
  # remove the corresponding transcript ID's and the Undetermined sample
  select(-"transcript_id(s)")

## read in the TPM matrices and fix naming order
tpm_mat <-
  read_tsv(
    paste0(raw_data_folder, "/rsem.merged.gene_tpm_deconvo.tsv")
  ) %>%
  # set NA's as 0's by default and enforce integers (for compact data format)
  mutate(across(where(is.numeric), ~ as.integer(replace_na(.x, 0)))) %>%
  # remove the corresponding transcript ID's and the Undetermined sample
  select(-"transcript_id(s)")

## read in the metadata and consolidate samplenames
pheno <-
  read_tsv(paste0(raw_data_folder, "/E-MTAB-15505.sdrf.txt")) %>%
  # rewire samples to fit "samplename", "cyclepase" (pre-post), "group" (uf-biopsy)
  transmute(
    "samplename" = .$"Source Name",
    "cyclephase" = .$"Characteristics[cycle phase]",
    "group" = .$"Characteristics[sampling site]",
  ) %>%
  # rename them back to previous annotation, so that following analysis steps rename again
  mutate(
    # add a group column based on eithe UF or tissue
    group = ifelse(str_detect(samplename, "UF"), "UF", "biopsy"),
    # rename groups to short group standards in this analysis
    cyclephase = case_when(
      cyclephase == "P" ~ "pro",
      cyclephase == "ES" ~ "pre",
      cyclephase == "MS" ~ "rec",
      cyclephase == "LS" ~ "post",
      T ~ cyclephase
    )
  ) %>%
  { .[match(colnames(raw_mat)[-1], .$samplename), ] }


#### FILTER
# S01, S22, S19, S04 were excluded due to poor clustering
# S22, S19, S04 were excluded due to low biotype proportion
# S07_UF was swapped with S07_tissue based on clustering, but removed in the end
removals <- c(
  "S22_UF", "S22_tissue", "S22_tissue_1", "S22_tissue_2", "S22_tissue_3",
  "S01_UF", "S01_tissue", "S07_tissue", "S07_UF",
  "S19_UF", "S19_tissue", "S04_UF", "S04_tissue"
)

# switch and filter some samples
raw_mat_unfilt <- raw_mat
raw_mat <-
  raw_mat %>%
  select(-removals[removals %in% colnames(raw_mat)])

# switch and filter some samples
tpm_mat_unfilt <- tpm_mat
tpm_mat <-
  tpm_mat %>%
  select(-removals[removals %in% colnames(tpm_mat)])


pheno_unfilt <- pheno
pheno <-
  pheno %>%
  # but remove an samples actually not in the raw_mat (supposedly 8 of them)
  # and also quarantee that the samplenames are in the same order
  # as raw_mat colnames (needed for DESeq2)
  { .[match(colnames(raw_mat)[-1], .$samplename), ] }

## join count matrix with pheno for visualization or downstream analysis
expr_mat <- raw_mat %>% transform_count_mat(pheno)
expr_tpm <- tpm_mat %>% transform_count_mat(pheno)


#### OUTPUT
## Emit formatted data files
data_subfolder <- paste0(data_folder, "/filtered/")
raw_mat %>% write_feather(paste0(data_subfolder, "counts_raw.feather"))
tpm_mat %>% write_feather(paste0(data_subfolder, "tpm_raw.feather"))
raw_mat_unfilt %>% write_feather(paste0(data_subfolder, "counts_raw_unfilt.feather"))
tpm_mat_unfilt %>% write_feather(paste0(data_subfolder, "tpm_raw_unfilt.feather"))
expr_mat %>% write_feather(paste0(data_subfolder, "counts_format.feather"))
expr_tpm %>% write_feather(paste0(data_subfolder, "tpm_format.feather"))
pheno %>% write_tsv(paste0(data_subfolder, "phenotype.tsv"))
pheno_unfilt %>% write_tsv(paste0(data_subfolder, "phenotype_unfilt.tsv"))

## also write out raw file with switched and formatted annotation
raw_mat %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  rename(gene_id = external_gene_name) %>%
  write_feather(paste0(data_subfolder, "annot_raw.feather"))

raw_mat_unfilt %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  rename(gene_id = external_gene_name) %>%
  write_feather(paste0(data_subfolder, "annot_raw_unfilt.feather"))

tpm_mat %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  rename(gene_id = external_gene_name) %>%
  write_feather(paste0(data_subfolder, "annot_tpm.feather"))

tpm_mat_unfilt %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  rename(gene_id = external_gene_name) %>%
  write_feather(paste0(data_subfolder, "annot_tpm_unfilt.feather"))
