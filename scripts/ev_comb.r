#### R script for reading in Vigano/NOTED/HUT datasets
## and harmonising the dataset for further downstream analysis

## 1. Load datasets and unify metadata
## 2. Join them and apply batch normalisation with sva-combatch
## 3. Look at the clustering of the joined dataset before and after sva

source("./scripts/load_deps.r")

# load the filtered set of samples and previous functions from backend
hut_mat <-
  read_feather(paste0(data_folder, "/filtered/counts_raw.feather")) %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  # set NA's as 0's by default and enforce integers (for compact data format)
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  filter(!is.na(external_gene_name))

hut_pheno <- read_tsv(paste0(data_folder, "/filtered/phenotype.tsv"))


# Also load in the Vigano EV dataset
ev_comb <-
  read_excel(
    paste0(raw_data_folder, "/EV_RNAseq_Vigano.xlsx"),
    range = "D3:AB54784"
  ) %>%
  # inner_join(ev_comb, by = "GeneID") %>%
  rename(external_gene_name = GeneID)

# combine the phenotype
ev_comb_pheno <-
  read_excel(
    paste0(raw_data_folder, "/EV_RNAseq_Vigano.xlsx"),
    range = "E3:AB3"
  ) %>% {
    tibble(
      samplename = colnames(.),
      cyclephase = c(rep("pre", 12), rep("rec", 12)),
      group = "UF"
    )
  } %>%
  left_join(read_csv(paste0(raw_data_folder, "/vigano_batches.csv")))


#### Join the datasets
# join the dataset itself and don't filter out 
comb_mat <-
  list("HUT" = hut_mat, "Vigano" = ev_comb) %>%
  reduce(inner_join, by = "external_gene_name") %>%
  select(!ends_with("biopsy_1")) # !ends_with("biopsy") & 

# harmonise the pheno file and add batch parameter
# grouping says the sample type (EV/Biopsy)
# cyclephase the grouping of the cycle timing (pre/rec/post)
comb_pheno <-
  list("HUT" = hut_pheno, "Vigano" = ev_comb_pheno) %>%
  bind_rows(.id = "dataset") %>%
  mutate(dataset = ifelse(dataset == "Vigano", batch, dataset)) %>%
  select(c("samplename", "cyclephase", "dataset", "group")) %>%
  # align the sample order with the matrix
  filter(samplename %in% colnames(comb_mat)) %>%
  arrange(match(samplename, colnames(comb_mat)[-1]))

comb_pheno_uf <- comb_pheno %>% filter(group == "UF") %>% select(-grou)


# apply batch normalisation on the dataframes
comb_batch_all <-
  comb_mat %>%
  column_to_rownames("external_gene_name") %>%
  as.matrix %>%
  sva::ComBat_seq(
    batch = comb_pheno$dataset,
    group = setNames(comb_pheno$cyclephase, comb_pheno$group)
  ) %>%
  as_tibble(rownames = "gene_id")

# apply batch normalisation on the dataframes - only UF
comb_batch_uf <-
  comb_mat %>%
  select(!ends_with("biopsy")) %>% # select only UF samples
  column_to_rownames("external_gene_name") %>%
  as.matrix %>%
  sva::ComBat_seq(
    batch = comb_pheno_uf$dataset,
    group = comb_pheno_uf$cyclephase
  ) %>%
  as_tibble(rownames = "gene_id")


#### write the expresiion matrix and the phenotype data out
write_feather(comb_batch_all, paste0(data_folder, "/combined/comb_all_batch.feather"))
write_feather(comb_batch_uf, paste0(data_folder, "/combined/comb_uf_batch.feather"))
write_feather(comb_mat, paste0(data_folder, "/combined/comb_all_raw.feather"))
write_tsv(comb_pheno_uf, paste0(data_folder, "/combined/comb_uf_pheno.tsv"))
write_tsv(comb_pheno, paste0(data_folder, "/combined/comb_all_pheno.tsv"))