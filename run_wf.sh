#!/bin/bash


## run the preprocessing with scripts in ./preproc_scripts/

## prepare the .env

## run processing and model training part - recommended to have a ~16GB VRAM GPU for this
conda run -n endo-ev-r RScript scripts/ev_raw.r
conda run -n endo-ev-r RScript scripts/ev_comb.r
conda run -n ev-proc quarto render scripts/preproc_sc.qmd
conda run -n ev-proc quarto render scripts/preproc_st.qmd

## the analysis steps - NB! run ev_deconvolution.qmd interactively first to generate accessory documents
conda run -n endo-ev quarto render analysis/ev_deconvolution.qmd --to html
conda run -n endo-ev quarto render analysis/ev2space_project.qmd --to html
conda run -n endo-ev quarto render analysis/ev2space_stats.qmd --to html
conda run -n endo-ev quarto render analysis/ev_article.qmd --to html
