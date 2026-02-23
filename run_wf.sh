#!/bin/bash


## run the preprocessing with scripts in ./preproc_scripts/

## prepare the .env

## run processing and model training part - recommended to have a ~16GB VRAM GPU for this
pixi run -e r-analysis RScript scripts/ev_raw.r
pixi run -e r-analysis RScript scripts/ev_comb.r
pixi run -e proc quarto render scripts/preproc_sc.qmd
pixi run -e proc quarto render scripts/preproc_st.qmd

## the analysis steps - NB! run ev_deconvolution.qmd interactively first to generate accessory documents
pixi run -n analysis quarto render analysis/ev_deconvolution.qmd --to html
pixi run -n analysis quarto render analysis/ev2space_project.qmd --to html
pixi run -n analysis quarto render analysis/ev2space_stats.qmd --to html
pixi run -n analysis quarto render analysis/ev_article.qmd --to html
