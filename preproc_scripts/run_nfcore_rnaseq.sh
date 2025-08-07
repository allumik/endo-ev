#!/bin/bash

#SBATCH --job-name=EVREM-mRNA
#SBATCH --mail-user=alvinmeltsov@gmail.com
#SBATCH --mail-type=ALL

#SBATCH --output=logs/%j.log
## these parameters are for the executor, not child processes
## for processes, check nextflow.config
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00

# you are gonna need to install micromamba and nf-core as environment per
# https://nf-co.re/docs/nf-core-tools/installation
# and then preload the rnaseq pipeline with
# nf-core download rnaseq -r 3.14.0 --outdir nf-core-rnaseq
micromamba activate nf-core

NXF_OPTS='-Xms1g -Xmx4g'
module load bioinf/java/17.0.8_7
module load singularity
## don't load nextflow module as its too old :(
## you have to byob
# module load nextflow

## nf-core/rnaseq v3.14 (https://nf-co.re/rnseq/3.14.0)
# --aligner star_rsem     <- alternatives [star_rsem|star_salmon|hisat2]
# --save_reference        <- Save the reference genomes used
# --save_align_intermeds  <- Save the intermediate BAM files
# --save_trimmed          <- I guess you figured it out already...
# -profile singularity    <- Use singularity containers
# nextflow run nf-core/rnaseq \
#    -r 3.14.0 \
nextflow run ./nfcore_rnaseq/3_14_0/ \
    --input samplesheet_mrna.csv \
    --genome GRCh37 \
    --aligner star_rsem \
    --email "alvinmeltsov@gmail.com" \
    --outdir "./results_mrna/" \
    -c nextflow.config \
    -profile singularity # \
    # -resume ## if you wish to resume, not start again
