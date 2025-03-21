#### Differential analysis runner for the ev-biomarker studies
## Depends on back.R and back_diff.R

source("./scripts/load_deps.r")
source("./scripts/de_functions.r")

## And the data too
data_subfolder <- paste0(data_folder, "/filtered/")
raw_mat <- read_feather(paste0(data_subfolder, "counts_raw.feather"))
pheno <- read_tsv(paste0(data_subfolder, "phenotype.tsv"))
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
    print(group_name)
    deseq_top <-
      dds_maker(
        raw_mat %>% select(c("gene_id", contains(group_name))),
        diffterm %>% filter(group == group_name),
        ~cyclephase,
        filter_ids = filter_ids
      ) %>%
      deseq_results(comps, annot, groupterm = "cyclephase")
    deseq_top_genes <- extract_top(deseq_top)
    return(c(deseq_top, deseq_top_genes))
  }
)

voom_res_phase <- lapply(
  named_groups,
  function(group_name) {
    print(group_name)
    voom_top <-
      edger_maker(
        raw_mat %>% select(c("gene_id", contains(group_name))),
        diffterm %>% filter(group == group_name),
        "cyclephase",
        filter_ids
      ) %>%
      voom_maker(comps) %>%
      voom_results(comps, annot)
    voom_top_genes <- extract_top(voom_top)
    return(c(voom_top, voom_top_genes))
  }
)


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
    deseq_top <-
      dds_maker(
        raw_mat %>% select(c("gene_id", samples_in_group[[group_name]])),
        diffterm %>% filter(cyclephase == group_name),
        ~ group,
        filter_ids = filter_ids
      ) %>%
      deseq_results(comps, annot, groupterm = "group")
    deseq_top_genes <- extract_top(deseq_top)
    return(c(deseq_top, deseq_top_genes))
  }
)

voom_res_group <- lapply(
  named_groups,
  function(group_name) {
    print(group_name)
    voom_top <-
      edger_maker(
        raw_mat %>% select(c("gene_id", samples_in_group[[group_name]])),
        diffterm %>% filter(cyclephase == group_name),
        "group",
        filter_ids
      ) %>%
      voom_maker(comps) %>%
      voom_results(comps, annot)
    voom_top_genes <- extract_top(voom_top)
    return(c(voom_top, voom_top_genes))
  }
)