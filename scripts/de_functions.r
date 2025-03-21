#### Backend scripts for the cERtain project
##

#### Differential analysis steps

## top genes boolean table generator
extract_top <- function(deg_list) {
  unlist(unname(sapply(deg_list, function(x) x$locus))) %>% {
    mutate(
      as_tibble(sapply(deg_list, function(x) . %in% x$locus)),
      locus = .
    )
  } %>%
    distinct(locus, .keep_all = T)
}

## terms for the standrard differential analysis
diff_terms <- function(phenotype_data) {
  phenotype_data %>%
    mutate(
      group = factor(group),
      cyclephase = factor(cyclephase)
    ) %>%
    column_to_rownames("samplename")
}

## do cpm filtering on raw counts with edgeR `cpm` function
## returns only the
cpm_filter <- function(counts_data, cpm_lim = 5, mean_or_sum = "mean") {
  counts_data %>%
    column_to_rownames("gene_id") %>%
    edgeR::DGEList() %>%
    edgeR::cpm() %>%
    ## Output thresholded boolean vector
    {
      if(mean_or_sum == "mean") rowMeans(.) > cpm_lim
      else if(mean_or_sum == "sum") rowSums(.) > cpm_lim
      else(print("Check your `mean_or_sum` argument"))
    }
}

## Do CPM filtering using grouped datasets with phenotype
## Outputs a tibble with CPM > cpm_limit boolean columns for each group
cpm_filter_per_group <-
  function(counts_data, named_samples, cpm_lim = 5, group_lim = 3) {
    tmp <- foreach::foreach(g = unique(names(named_samples)), .combine = cbind, .inorder = T) %do% {
      counts_data %>%
        select(c("gene_id", as.character(named_samples[names(named_samples) == g]))) %>%
        cpm_filter(cpm_lim)
    }
    colnames(tmp) <- unique(names(named_samples))
    return(as_tibble(tmp, rownames = "gene_id"))
  }

## Conviniencs function for filtering out genes with at least > cpm_lim counts in any group
## Outputs vector
cpm_filter_grouped <-
  function(counts_data, named_samples, cpm_lim = 5, group_lim = 3) {
    cpm_filter_per_group(counts_data, named_samples, cpm_lim, group_lim) %>%
      column_to_rownames("gene_id") %>%
      rowSums(.) > 1
  }



## DESEq2
## Prepare DESeq2 object
dds_maker <- function(counts_data, phenotype_data, criteria, ...) {
  DESeq2::DESeqDataSetFromMatrix(
    countData =
      counts_data %>% column_to_rownames("gene_id"),
    colData = phenotype_data,
    design = criteria
  ) %>%
    {
      ## filter genes with lower total amount of cpm_lim
      if (hasArg(filter_ids)) .[filter_ids, ] else (.)
    } %>%
    # parallel does not work with 16G of RAM
    DESeq2::DESeq(parallel = F)
}

# get DESeq2 object results in the form of DEG table with other stuff
deseq_results <- function(
  dds_object, 
  comparison_list, 
  annotation, 
  pval = 0.05, 
  groupterm = "group",
  annotation_join_on = "ensembl_gene_id"
) {
  lapply(
    comparison_list,
    function(x) {
      DESeq2::results(dds_object, contrast = c(groupterm, x)) %>%
        as_tibble(rownames = annotation_join_on) %>%
        filter(padj < pval) %>%
        select(-c("lfcSE", "stat")) %>%
        left_join(
          annotation,
          by = annotation_join_on
        ) %>%
        rename("locus" = all_of(annotation_join_on))
    }
  )
}


## EdgeR
edger_maker <- function(counts_data, phenotype_data, criteria, filter_ids) {
  edge_mat <-
    edgeR::DGEList(
      counts = counts_data %>% column_to_rownames("gene_id"),
      samples = phenotype_data,
      group = phenotype_data[, criteria]
    )

  groups <- phenotype_data[, criteria]
  ## calculate normalization factors and estiamte dispersion
  # edge_mat[edgeR::filterByExpr(edge_mat), keep.lib.sizes = F] %>%
  edge_mat[filter_ids, keep.lib.sizes = T] %>%
    edgeR::calcNormFactors() %>%
    edgeR::estimateDisp(model.matrix(~ 0 + groups), robust = T)
}


## Diff analysis edger
edger_results <- function(
  edge_obj,
  comparison_list,
  annotation,
  pval = 0.05,
  annotation_join_on = "ensembl_gene_id"
) {
  lapply(
    comparison_list,
    function(x)
      edgeR::exactTest(edge_obj, pair = x) %>%
        ## select by FDR
        edgeR::topTags(n = Inf, p.value = pval) %>%
        {
          # handle empty table clause too
          if (nrow(.) == 0) {
            tibble(!!annotation_join_on := character()) %>%
              left_join(annotation, by = annotation_join_on) %>%
              rename(locus = all_of(annotation_join_on))
          } else {
            .$table %>%
              as_tibble(rownames = annotation_join_on) %>%
              left_join(annotation, by = annotation_join_on) %>%
              rename(locus = all_of(annotation_join_on))
          }
        }
  )
}



#### limma-voom
## apply voom on reused edgeR matrix
voom_maker <- function(edge_obj, comparison_list) {
  voom_fit <-
    limma::voom(edge_obj, edge_obj$design) %>%
    limma::lmFit(edge_obj$design)

  ## Diff analysis voom
  ## specify contrasts
  voom_diff <-
    limma::makeContrasts(
      contrasts =
        sapply(unname(comparison_list), function(x) paste(x, collapse = "-")),
      levels =
        colnames(coef(voom_fit)) %>% str_remove("groups")
    )

  ## rename the colnames for the next step
  colnames(voom_fit$coefficients) <- rownames(voom_diff)

  ## estimate contrasts and perform bayesian smoothing
  voom_fit %>%
    limma::contrasts.fit(voom_diff) %>%
    limma::eBayes()
}

voom_results <- function(
  voom_obj,
  comparison_list,
  annotation,
  pval = 0.05,
  annotation_join_on = "ensembl_gene_id"
) {
  lapply(
    setNames(seq(length(comparison_list)), names(comparison_list)),
    function(x)
      limma::topTreat(
        voom_obj,
        coef = x,
        number = Inf,
        adjust.method = "BH",
        p.value = pval
      ) %>%
        as_tibble(rownames = annotation_join_on) %>%
        left_join(
          annotation,
          by = annotation_join_on
        ) %>%
        rename("locus" = all_of(annotation_join_on))
  )
}
