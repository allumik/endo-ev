# Endometrial receptivity profiling and deconvolution of Extracellular Vesicules

This study builds on previous work done by Vigano et al (TODO: add cite) and by Meltsov et al. (TODO: add cite) to move towards a non-invasive method for receptivity testing. Best efforts were made to try to organise things in the following way:

* Scripts to run nf-core/rnaseq preprocessing pipeline on the raw samples in the SLURM HPC.

* R scripts to take the read count matrices emitted by nf-core in `raw_data/` and output phenotype files after formatting to `data/`.

* RMarkdown files to generate interactive analysis reports in `analysis/`.

Environment manager of choice was `micromamba`, rebuild `renv-ev` environment with (TODO: the micromamba command).
