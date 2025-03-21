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
  cpm_filter_grouped(comb_batch, setNames(comb_pheno$samplename, comb_pheno$group))

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
    deseq_top <-
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
      )
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
        comb_batch %>% select(c("gene_id", diffterm %>% filter(group == group_name) %>% rownames)),
        diffterm %>% filter(group == group_name),
        "cyclephase",
        filter_ids
      ) %>%
      voom_maker(comps) %>%
      voom_results(comps, annot, annotation_join_on = "external_gene_name", pval = pval)
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
      )
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
      )
    voom_top_genes <- extract_top(voom_top)
    return(c(voom_top, voom_top_genes))
  }
)


#### Formatting some data for the raport
## transform some matrices to HUGO
expr_mat_test <-
  raw_mat_test %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  # set NA's as 0's by default and enforce integers (for compact data format)
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  filter(!is.na(external_gene_name)) %>%
  rename(gene_id = external_gene_name) %>%
  transform_count_mat(pheno_test)

expr_tpm_test <-
  raw_tpm_test %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  # set NA's as 0's by default and enforce integers (for compact data format)
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  filter(!is.na(external_gene_name)) %>%
  rename(gene_id = external_gene_name) %>%
  transform_count_mat(pheno_test)


## re/precalculate some large nasty bits
global_expr_pca <-
  comb_batch_tpm %>%
  transform_count_mat(comb_pheno) %>%
  bind_rows(expr_tpm_test[colnames(expr_tpm_test) %in% colnames(.)]) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  recipe() %>%
  update_role(samplename, new_role = "id") %>%
  update_role(c("group", "cyclephase"), new_role = "other") %>%
  ## eliminate some genes with zero expression
  step_zv(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 3) %>%
  prep(strings_as_factors = F) %>%
  bake(new_data =
    comb_batch_tpm %>% 
    transform_count_mat(comb_pheno) %>%
    bind_rows(expr_tpm_test[colnames(expr_tpm_test) %in% colnames(.)]) %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  )

global_expr_umap <-
  comb_batch_tpm %>%
  transform_count_mat(comb_pheno) %>%
  bind_rows(expr_tpm_test[colnames(expr_tpm_test) %in% colnames(.)]) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  recipe() %>%
  update_role(samplename, new_role = "id") %>%
  update_role(c("group", "cyclephase"), new_role = "other") %>%
  ## eliminate some genes with zero expression
  step_zv(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
  step_umap(all_numeric(), num_comp = 2) %>%
  prep(strings_as_factors = F) %>%
  bake(new_data =
    comb_batch_tpm %>%
    transform_count_mat(comb_pheno) %>%
    bind_rows(expr_tpm_test[colnames(expr_tpm_test) %in% colnames(.)]) %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  )



uf_out_join <-
  unique(c(
    deseq_res_phase$UF$pro_vs_pre$locus,
    deseq_res_phase$UF$pre_vs_rec$locus,
    deseq_res_phase$UF$rec_vs_post$locus,
    voom_res_phase$UF$pro_vs_pre$locus,
    voom_res_phase$UF$pre_vs_rec$locus,
    voom_res_phase$UF$rec_vs_post$locus
  )) %>% {
    tibble(
      locus = .,
      deseq_pro_vs_pre = . %in% deseq_res_phase$UF$pro_vs_pre$locus,
      deseq_pre_vs_rec = . %in% deseq_res_phase$UF$pre_vs_rec$locus,
      deseq_rec_vs_post = . %in% deseq_res_phase$UF$rec_vs_post$locus,
      voom_pro_vs_pre = . %in% voom_res_phase$UF$pro_vs_pre$locus,
      voom_pre_vs_rec = . %in% voom_res_phase$UF$pre_vs_rec$locus,
      voom_rec_vs_post = . %in% voom_res_phase$UF$rec_vs_post$locus
    )
  }


bio_out_join <-
  unique(c(
    deseq_res_phase$biopsy$pro_vs_pre$locus,
    deseq_res_phase$biopsy$pre_vs_rec$locus,
    deseq_res_phase$biopsy$rec_vs_post$locus,
    voom_res_phase$biopsy$pro_vs_pre$locus,
    voom_res_phase$biopsy$pre_vs_rec$locus,
    voom_res_phase$biopsy$rec_vs_post$locus
  )) %>% {
    tibble(
      locus = .,
      deseq_pro_vs_pre = . %in% deseq_res_phase$biopsy$pro_vs_pre$locus,
      deseq_pre_vs_rec = . %in% deseq_res_phase$biopsy$pre_vs_rec$locus,
      deseq_rec_vs_post = . %in% deseq_res_phase$biopsy$rec_vs_post$locus,
      voom_pro_vs_pre = . %in% voom_res_phase$biopsy$pro_vs_pre$locus,
      voom_pre_vs_rec = . %in% voom_res_phase$biopsy$pre_vs_rec$locus,
      voom_rec_vs_post = . %in% voom_res_phase$biopsy$rec_vs_post$locus
    )
  }


out_join <-
  unique(c(
    deseq_res_group$pro$UF_vs_biopsy$locus,
    deseq_res_group$pre$UF_vs_biopsy$locus,
    deseq_res_group$rec$UF_vs_biopsy$locus,
    deseq_res_group$post$UF_vs_biopsy$locus,
    voom_res_group$pro$UF_vs_biopsy$locus,
    voom_res_group$pre$UF_vs_biopsy$locus,
    voom_res_group$rec$UF_vs_biopsy$locus,
    voom_res_group$post$UF_vs_biopsy$locus
  )) %>% {
    tibble(
      locus = .,
      deseq_pro = . %in% deseq_res_group$pro$UF_vs_biopsy$locus,
      deseq_pre = . %in% deseq_res_group$pre$UF_vs_biopsy$locus,
      deseq_rec = . %in% deseq_res_group$rec$UF_vs_biopsy$locus,
      deseq_post = . %in% deseq_res_group$post$UF_vs_biopsy$locus,
      voom_pro = . %in% voom_res_group$pre$UF_vs_biopsy$locus,
      voom_pre = . %in% voom_res_group$pre$UF_vs_biopsy$locus,
      voom_rec = . %in% voom_res_group$rec$UF_vs_biopsy$locus,
      voom_post = . %in% voom_res_group$post$UF_vs_biopsy$locus
    )
  }


### OUTPUT
# write out the DE results