#!/bin/bash

#PBS -N aggregate
#PBS -o job_reports/aggregate-output
#PBS -e job_reports/aggregate-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -q ghlc
#PBS -t 1-91

set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

Rscript ~/repo/methylation_residuals/aggregate_hsq.R ${i}
