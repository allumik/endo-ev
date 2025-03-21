#### Settings and functions for the environment

# set this manually if it gets overwritten during load
select <- dplyr::select

## convert ENSG datasets to common gene names
geneid_converter <- function(dat, biomart_id) {
  dat %>%
    left_join(
      biomart_id %>% select(c("external_gene_name", "ensembl_gene_id")),
      by = "ensembl_gene_id"
    ) %>%
    select(-ensembl_gene_id) %>%
    filter(!is.na("external_gene_name")) %>%
    ## Then we summarize gene counts per ensembl ID
    pivot_longer(
      where(is_bare_numeric) | where(is.logical),
      names_to = "samplename",
      values_to = "raw_counts"
    ) %>%
    group_by(external_gene_name, samplename) %>%
    summarize(raw_counts = sum(raw_counts)) %>%
    ungroup %>%
    pivot_wider(
      names_from = "samplename",
      values_from = "raw_counts"
    )
}

transform_count_mat <- function(
  count_mat, pheno, id_columns = c("samplename", "cyclephase", "group")
) {
  count_mat %>%
    ## select only genes with more than 0 counts overall
    filter(rowSums(across(where(is.numeric))) > 0) %>%
    pivot_longer(
      where(is.numeric),
      names_to = "samplename",
      values_to = "counts"
    ) %>%
    ## add the pheno data by inner join (ie only samples in pheno file)
    inner_join(pheno, by = "samplename") %>%
    ## and turn it into a expr table where samples x genes
    pivot_wider(
      id_cols = id_columns,
      names_from = "gene_id",
      values_from = "counts"
    )
}


join_model_pca <- function(pca_data, score_table) {
  #' This function takes the PCA data from the `beready_model` method
  #' and joins it with the results table for visualization.
  full_join(pca_data, score_table, by = c("sample" = "Sample Name")) %>%
    rename(label = `Predicted Group`) %>%
    mutate(
      label =
        case_when(
          group == "pre" ~ "Pre-receptive",
          group == "rec" ~ "Receptive",
          group == "post" ~ "Post-receptive",
          T ~ label),
      group = ifelse(group == "Input", "Input", "Reference")
    )
}

#### Apply TPM normalization
## 1. Gene counts divide with gene length.
## 2. Sum all the ratios from divisions.
## 3. Divide each transcript ratio with ratio sum and then multiply by 1000000.
## Micheal Love form: https://support.bioconductor.org/p/91218/
tpm_norm <- function(counts, gene_len) {
  x <- as.matrix(counts) / gene_len
  return(t(t(x) * 1e6 / colSums(x)))
}

cpm_norm <- function(counts) {
  return(t(t(counts) * 1e6 / colSums(counts)))
}

## Get gene lengths from the GTF reference
## adapted from: https://www.biostars.org/p/91218/#91233
gene_len_file <- paste0(dirname(data_folder), "/ref_genomes/gencode.v27.annotation.gtf.gz")
if (!file.exists(gene_len_file))
  download.file(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz",
    gene_len_file
  )

gene_len <-
  # Load in the GTF GRCh38 GenCode 27 reference
  rtracklayer::import.gff(
    gene_len_file,
    format = "gtf",
    feature.type = "exon"
  ) %>%
  { S4Vectors::split(., .@elementMetadata$gene_name) } %>%
  IRanges::reduce() %>%
  unlist(use.names = T) %>%
  ## get the ranges element from GRanges
  .@ranges %$%
  tibble("gene_id" = .@NAMES, "width" = .@width) %>%
  distinct %>%
  ## summarize the transcript level sizes
  group_by(gene_id) %>%
  summarize(len = sum(width))


## Download biomaRt annotations for all the genes
# and save them for quick access later on
annot_file <- paste0(raw_data_folder, "/annot_table.csv")
if (!file.exists(annot_file))
  biomaRt::useMart(
    "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "https://grch37.ensembl.org"
  ) %>% {
    biomaRt::getBM(
      mart = .,
      attributes = c("ensembl_gene_id", "hgnc_symbol", "description")
    )
  } %>%
    rename(external_gene_name = "hgnc_symbol") %>%
    write_csv(annot_file)

annot <- read_csv(annot_file)