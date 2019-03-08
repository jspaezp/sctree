#!/bin/bash -l 

export R_LIBS=~/Rlibs:$R_LIBS

# Call it from the command line and 
# will use the first argument as the parent RDS

SEURAT_RDS=$1
MARKERS_RDS=$(echo ${SEURAT_RDS} | sed "s/seurat_/markers_/g")
DATASET_ID=$(echo ${SEURAT_RDS} | sed "s/seurat_//g" | sed "s/.RDS//g")

[[ -z $SEURAT_RDS ]] && echo "Please provide an argument (seurat object in disk) to the script" && exit

set -x 
set -e

VARIABLES=$(R --silent --slave -e "cat(unique(make.names(readRDS('${SEURAT_RDS}')@ident)))")

mkdir -p logs
mkdir -p built_jobs

for var in ${VARIABLES};
do
  echo "${var}"
  JOBNAME=$(echo "./built_jobs/_${DATASET_ID}${var}_BUILT_JOB_.bash")
  sed \
    -e "s/{{CLUSTER_ID}}/${var}/g" \
    -e "s/{{DATASET_ID}}/${DATASET_ID}/g" \
    -e "s/{{SEURAT_RDS}}/${SEURAT_RDS}/g" \
    -e "s/{{MARKERS_RDS}}/${MARKERS_RDS}/g" \
    ./templates/job_template_acc_benchmark.bash > "${JOBNAME}"
  echo "${JOBNAME}"
  qsub -e ./logs -o ./logs -l nodes=1:ppn=24 -l walltime=4:00:00 "${JOBNAME}"
done

