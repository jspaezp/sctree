#!/bin/bash -l

cd $PBS_O_WORKDIR
TOKEN="734551512:AAEpkwbPaR-T0Vqr6FL3TBKKvhXm3J6prEo"
CHAT_ID="455502653"

export R_LIBS=~/Rlibs:$R_LIBS
module load gcc hdf5

BASE_FILE=seurat_std_workflow.Rmd
DIR_10X={{DIR_10X}}
OUT_FILENAME={{OUT_FILENAME}}
BASE_OUT_ID=$(echo "${OUT_FILENAME}" | sed "s/\W//g")

mkdir -p build_dir
mkdir -p output_dir
TMP_INTERMEDIATES_DIR=$(mktemp -d -p build_dir)

wget "https://api.telegram.org/bot${TOKEN}/sendMessage" \
     --post-data "chat_id=${CHAT_ID}&text='JOB_STARTED ${OUT_FILENAME}'" \
     --quiet -O /dev/null

R --no-restore-data --no-save -e "rmarkdown::render( \
  '${BASE_FILE}', clean = FALSE, \
  intermediates_dir = '${TMP_INTERMEDIATES_DIR}', \
  output_file= '${BASE_OUT_ID}_.html', \
  output_dir = './output_dir/', \
  params = list( \
  out_filename = '${OUT_FILENAME}', \
  directory_10x = '${DIR_10X}'))"

wget "https://api.telegram.org/bot${TOKEN}/sendMessage" \
     --post-data "chat_id=${CHAT_ID}&text='JOB_FINISHED ${OUT_FILENAME}'" \
     --quiet -O /dev/null
