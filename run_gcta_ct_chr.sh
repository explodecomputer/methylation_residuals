#!/bin/bash

#PBS -N methres
#PBS -o job_reports/methres-output
#PBS -e job_reports/methres-error
#PBS -t 1-2000
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash
#PBS -q veryshort

set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

script="${HOME}/repo/methylation_residuals/run_gcta_ct_chr.R"
split="1000"

R --no-save --args ${i} ${split} < ${script}
