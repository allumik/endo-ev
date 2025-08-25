#### Differential analysis runner for the ev-biomarker studies
## Depends on back.R and back_diff.R

source("./scripts/load_deps.r")
source("./scripts/de_functions.r")

## And the data too
data_subfolder <- paste0(data_folder, "/filtered/")
raw_mat <- read_feather(paste0(data_subfolder, "counts_raw.feather"))
pheno <- 
  read_tsv(paste0(data_subfolder, "phenotype.tsv")) %>% 
  # rename the group to tissue for compatilibity 
  mutate(group = ifelse(str_detect(samplename, "UF"), "UF", "tissue"))
annot <- read_csv(paste0(raw_data_folder, "/annot_table.csv"))

## Diff analysis comparisons within `group`
comps <- list(
  "pro_vs_pre" = c("pro", "pre"),
  "pre_vs_rec" = c("pre", "rec"),
  "rec_vs_post" = c("rec", "post")
)

pval <- .05
diffterm <- diff_terms(pheno)

## Filter out counts with edgeR's cpm method
filter_ids <-
  cpm_filter_grouped(raw_mat, setNames(pheno$samplename, pheno$group))

named_groups <-
  setNames(
    as.character(unique(diffterm$group)),
    as.character(unique(diffterm$group))
  )

## generate vst normalised counts
expr_vst <-
  dds_maker(raw_mat, diffterm, ~ group) %>%
  DESeq2::vst(blind = T) %>%
  SummarizedExperiment::assay() %>%
  t %>%
  as_tibble(rownames = "sample")

# for UF and biopsy, go over
deseq_res_phase <- lapply(
  named_groups,
  function(group_name) {
    dds_maker(
      raw_mat %>% select(c("gene_id", contains(group_name))),
      diffterm %>% filter(group == group_name),
      ~cyclephase,
      filter_ids = filter_ids
    ) %>%
      deseq_results(comps, annot, groupterm = "cyclephase") %>%
      bind_rows(.id = "comparison")
  }
) %>%
  bind_rows(.id = "methods")

voom_res_phase <- lapply(
  named_groups,
  function(group_name) {
    edger_maker(
      raw_mat %>% select(c("gene_id", contains(group_name))),
      diffterm %>% filter(group == group_name),
      "cyclephase",
      filter_ids
    ) %>%
      voom_maker(comps) %>%
      voom_results(comps, annot) %>%
      bind_rows(.id = "comparison")
  }
) %>%
  bind_rows(.id = "methods")


# also run the differential analysis for between the UF and Biopsy
# but same cyclephase group
named_groups <-
  setNames(
    as.character(unique(diffterm$cyclephase)),
    as.character(unique(diffterm$cyclephase))
  )

samples_in_group <- lapply(
  named_groups,
  function(cyclephase_name)
    diffterm %>% filter(cyclephase == cyclephase_name) %>% rownames
)

comps <- list(
  "UF_vs_biopsy" = c("UF", "tissue")
)

deseq_res_group <- lapply(
  named_groups,
  function(group_name) {
    dds_maker(
      raw_mat %>% select(c("gene_id", samples_in_group[[group_name]])),
      diffterm %>% filter(cyclephase == group_name),
      ~ group,
      filter_ids = filter_ids
    ) %>%
      deseq_results(comps, annot, groupterm = "group") %>%
      .[[1]] # unwrap the first element
  }
) %>%
  bind_rows(.id = "comparison")

voom_res_group <- lapply(
  named_groups,
  function(group_name) {
    edger_maker(
      raw_mat %>% select(c("gene_id", samples_in_group[[group_name]])),
      diffterm %>% filter(cyclephase == group_name),
      "group",
      filter_ids
    ) %>%
      voom_maker(comps) %>%
      voom_results(comps, annot) %>%
      .[[1]] # unwrap the first element
  }
) %>%
  bind_rows(.id = "comparison")

### OUTPUT
# write out the DE result
de_folder <- paste0(data_subfolder, "/de/")
if (!dir.exists(de_folder)) dir.create(de_folder)
deseq_res_phase %>% 
  mutate(methods = ifelse(methods == "tissue", "biopsy", methods)) %>%
  write_feather(paste0(de_folder, "deseq_phases.feather"))
voom_res_phase %>% 
  mutate(methods = ifelse(methods == "tissue", "biopsy", methods)) %>%
  write_feather(paste0(de_folder, "voom_phases.feather"))
deseq_res_group %>% write_feather(paste0(de_folder, "deseq_methods.feather"))
voom_res_group %>% write_feather(paste0(de_folder, "voom_methods.feather"))
