# Some prelim lines to gather the fastq's from dragen working directories.
# (rsync -nv) for dry-run
# rsync -uv [folder/with/files/**/*.fastq.gz] .
# rm ./Undetermined*

## This file is part of nf-core/rnaseq pipeline: https://nf-co.re/rnaseq/
# This script will look at the folder with fastq files and autogenerate samplesheet

module load python

FASTQDIR="./fastq_mrna/"
FILE="fastq_dir_to_samplesheet.py"
SAMPLESHEET="samplesheet_mrna.csv"
 
## Download the script if file not present
[ ! -f "$FILE" ] && wget -L \
  https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py

## Run the script
chmod u+x $FILE
./$FILE $FASTQDIR $SAMPLESHEET --strandedness reverse

# and replace the "mRNA" in "name" column (1st) with "biopsy"
sed -i 's/mRNA/biopsy/1#' $SAMPLESHEET
sed -i 's/_S[[:digit:]]\+//1#' $SAMPLESHEET

# add some thing manually, such as "_1/2/3" for the control samples etc
