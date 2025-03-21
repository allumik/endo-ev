#!/bin/bash

#SBATCH --job-name=EVREM-miRNA
#SBATCH --mail-user=alvinmeltsov@gmail.com
#SBATCH --mail-type=ALL

#SBATCH --output=logs/%j.log
#SBATCH --partition=main
## these parameters are for the executor, not child processes
## for processes, check nextflow.config
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00

NXF_OPTS='-Xms1g -Xmx4g'
module load java
module load singularity/3.8.5
## don't load nextflow module as its too old :(
## you have to byob
# module load nextflow

## nf-core/smrnaseq v2.2.3 (https://nf-co.re/smrnaseq/2.2.3)
# --aligner star_rsem     <- alternatives [star_rsem|star_salmon|hisat2]
# --protocol		  <- illumina/nextfleq/qiaseq etc
# --filter_contamination  <- filter out cdna/rrna/trna etc
# --save_reference        <- Save the reference genomes used
# --save_align_intermeds  <- Save the intermediate BAM files
# --save_trimmed          <- I guess you figured it out already...
# -profile singularity    <- Use singularity containers
nextflow run nf-core/smrnaseq \
    -r 2.2.3 \
    --input samplesheet_mirna.csv \
    --genome GRCh37 \
    --aligner star_rsem \
    --protocol "illumina" \
    --filter_contamination \
    --cdna "./ref_data/" /
    --save_reference \
    --save_trimmed \
    --save_align_intermeds \
    --email "alvinmeltsov@gmail.com" \
		--outdir "./results_mirna/" \
    -c "./nextflow.config" \
    -profile singularity
    #-resume ## if you wish to resume, not start
