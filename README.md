# Non-invasive Transcriptomic Cell Profiling of the Human Endometrium with Generative Deep Learning


## Abstract

This repository contains scripts and Quarto documents for the analysis of uterine fluid extracellular vesicles (UF-EV) using generative deep models, specifically the BulkTrajBlend architecture from `omicverse` package.


## Repository Structure

Best efforts were made to try to organise things in the following way:

* `./preproc_scripts/` - Scripts to run nf-core/rnaseq preprocessing pipeline on the raw samples in a SLURM HPC.

* `./analysis/` - RMarkdown and Quarto documents to generate interactive analysis reports. The main analysis workflow used in the manuscript is described in the `ev_article.qmd`, while other documents in that folder were used for explorative data analysis and to run statistical testing.

* `./scripts/` - Scripts for processing raw read counts emitted by the nf-core/rnaseq pipeline.

  * `ev_raw.r` & `ev_comb.r` are used for preprocessing the read tables for our dataset and combined.

  * `de_runner.r` & `de_comb_runner.r` are used for running differential analysis (not included in the manuscript).

  * `preproc_sc.qmd` & `preproc_st.qmd` are document versions of the scripts used to run single cell atlas preprocessing, loading in the UF-EV datasets, training models and running inference for deconvolution and project to spatial transcriptomic datasets.

* R and Python scripts to take the read count matrices emitted by nf-core/rnaseq pipeline in `$RAW_DATA_FOLDER` (not included in this repository) and output phenotype files after formatting to `$DATA_FOLDER`.


## Setup

1. Using a `conda`-like environment manager, recreate the analysis environment with

```
conda env create -f env.yml -y
conda env create -f ev-proc.yml -y
conda env create -f renv.yml -y
```

2. Then, copy the `.env_template` to `.env` and populate the environment variables to suit your situation. 

3. Next, start running the preprocessing scripts found in `./preproc_scripts/` to quantify reads.

4. Follow the `./run_wf.sh` script.


## Citations

