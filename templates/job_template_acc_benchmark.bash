#!/bin/bash -l

cd $PBS_O_WORKDIR
TOKEN="734551512:AAEpkwbPaR-T0Vqr6FL3TBKKvhXm3J6prEo"
CHAT_ID="455502653"

export R_LIBS=~/Rlibs:$R_LIBS
module load gcc hdf5

BASE_FILE=./benchmarking_ml.Rmd
SEURAT_RDS={{SEURAT_RDS}}
CLUSTER_ID={{CLUSTER_ID}}
DATASET_ID={{DATASET_ID}}
MARKERS_RDS={{MARKERS_RDS}}

mkdir -p build_dir
mkdir -p output_dir
TMP_INTERMEDIATES_DIR=$(mktemp -d -p build_dir)

wget "https://api.telegram.org/bot${TOKEN}/sendMessage" \
     --post-data "chat_id=${CHAT_ID}&text='JOB_STARTED {{SEURAT_RDS}} {{CLUSTER_ID}}'" \
     --quiet -O /dev/null

R --no-restore-data --no-save -e "rmarkdown::render( \
  '${BASE_FILE}', clean = FALSE, \
  intermediates_dir = '${TMP_INTERMEDIATES_DIR}', \
  output_file= '{{DATASET_ID}}_CLUSTER_{{CLUSTER_ID}}.html', \
  output_dir = './output_dir/', \
  params = list( \
  SEURAT_RDS = '{{SEURAT_RDS}}', \
  MARKERS_RDS = '{{MARKERS_RDS}}', \
  GLOBAL_REPEATS = 2, \
  NUM_CROSS_VAL = 5, \
  DATA_FRACTION = 1, \
  FILTER_MEMBRANE = TRUE, \
  NUMTHREADS = 20, \
  DATASET_ID = '{{DATASET_ID}}', \
  CLUSTER = '{{CLUSTER_ID}}'))"

wget "https://api.telegram.org/bot${TOKEN}/sendMessage" \
     --post-data "chat_id=${CHAT_ID}&text='JOB_FINISHED {{SEURAT_RDS}} {{CLUSTER_ID}}'" \
     --quiet -O /dev/null
