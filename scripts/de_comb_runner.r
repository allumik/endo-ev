#### Differential analysis runner for the ev-biomarker studies, combined dataset

source("./scripts/load_deps.r")
source("./scripts/de_functions.r")

## And the data too
data_subfolder <- paste0(data_folder, "/combined/")
comb_batch <- read_feather(paste0(data_subfolder, "comb_uf.feather"))
comb_pheno <- read_tsv(paste0(data_subfolder, "comb_uf_pheno.tsv"))
annot <- read_csv(paste0(raw_data_folder, "/annot_table.csv"))

## Diff analysis comparisons within `group`
comps <- list(
  "pro_vs_pre" = c("pro", "pre"),
  "pre_vs_rec" = c("pre", "rec"),
  "rec_vs_post" = c("rec", "post")
)

pval <- .05

diffterm <- diff_terms(comb_pheno)

## Filter out counts with edgeR's cpm method
filter_ids <-
  cpm_filter_grouped(
    comb_batch,
    setNames(comb_pheno$samplename, comb_pheno$group)
  )

named_groups <-
  setNames(
    as.character(unique(diffterm$group)),
    as.character(unique(diffterm$group))
  )


## generate vst normalised counts
expr_vst <-
  dds_maker(comb_batch, diffterm, ~ group) %>%
  DESeq2::vst(blind = T) %>%
  SummarizedExperiment::assay() %>%
  t %>%
  as_tibble(rownames = "sample")

# for UF and biopsy separate, go over the cyclephases
deseq_res_phase <- lapply(
  named_groups,
  function(group_name) {
    print(group_name)
    dds_maker(
      comb_batch %>% select(c("gene_id", diffterm %>% filter(group == group_name) %>% rownames)),
      diffterm %>% filter(group == group_name),
      ~cyclephase,
      filter_ids = filter_ids
    ) %>%
      deseq_results(
        comps,
        annot,
        groupterm = "cyclephase",
        annotation_join_on = "external_gene_name",
        pval = pval
      ) %>%
      bind_rows(.id = "comparison")
  }
) %>%
  bind_rows(.id = "methods")

voom_res_phase <- lapply(
  named_groups,
  function(group_name) {
    print(group_name)
    voom_top <-
      edger_maker(
        comb_batch %>% select(c("gene_id", diffterm %>% filter(group == group_name) %>% rownames)),
        diffterm %>% filter(group == group_name),
        "cyclephase",
        filter_ids
      ) %>%
      voom_maker(comps) %>%
      voom_results(
        comps,
        annot,
        annotation_join_on = "external_gene_name",
        pval = pval
      ) %>%
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
  "UF_vs_biopsy" = c("UF", "biopsy")
)

deseq_res_group <- lapply(
  named_groups,
  function(group_name) {
    print(group_name)
    dds_maker(
      comb_batch %>% select(c("gene_id", samples_in_group[[group_name]])),
      diffterm %>% filter(cyclephase == group_name),
      ~ group,
      filter_ids = filter_ids
    ) %>%
      deseq_results(
        comps,
        annot,
        groupterm = "group",
        annotation_join_on = "external_gene_name",
        pval = pval
      ) %>%
      .[[1]] # unwrap the first element
  }
) %>%
  bind_rows(.id = "comparison")

voom_res_group <- lapply(
  named_groups,
  function(group_name) {
    print(group_name)
    edger_maker(
      comb_batch %>% select(c("gene_id", samples_in_group[[group_name]])),
      diffterm %>% filter(cyclephase == group_name),
      "group",
      filter_ids
    ) %>%
      voom_maker(comps) %>%
      voom_results(
        comps,
        annot,
        annotation_join_on = "external_gene_name",
        pval = pval
      ) %>%
      .[[1]] # unwrap the first element
  }
) %>%
  bind_rows(.id = "comparison")

### OUTPUT
# write out the DE result
de_folder <- paste0(data_subfolder, "/de")
if (!dir.exists(de_folder)) dir.create(de_folder)
deseq_res_phase %>% write_feather(paste0(de_folder, "deseq_phases.feather"))
voom_res_phase %>% write_feather(paste0(de_folder, "voom_phases.feather"))
deseq_res_group %>% write_feather(paste0(de_folder, "deseq_methods.feather"))
voom_res_group %>% write_feather(paste0(de_folder, "voom_methods.feather"))
