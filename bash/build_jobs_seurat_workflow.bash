#!/bin/bash -l 

# Call it from the command line and 
# will use the first argument as the parent RDS
# 

ERRORMSG="Please Provide 2 arguments, first the directory to read \
  the 10x dataset from and then the name that will be used for the outputs"

DIR_10X=$1
OUT_FILENAME=$2

[[ -z $DIR_10X ]] && echo "${ERRORMSG}" && exit
[[ -z $OUT_FILENAME ]] && echo "${ERRORMSG}" && exit

set -x 
set -e

mkdir -p logs
mkdir -p built_jobs

sed -e "s+{{OUT_FILENAME}}+${OUT_FILENAME}+g" \
  -e "s+{{DIR_10X}}+${DIR_10X}+g" \
  ./templates/job_template_seurat_workflow.bash > \
  "./built_jobs/BUILT_JOB_${OUT_FILENAME}.bash"

echo "Submitting BUILT_JOB_${OUT_FILENAME}.bash"
qsub -e ./logs -o ./logs "./built_jobs/BUILT_JOB_${OUT_FILENAME}.bash"

