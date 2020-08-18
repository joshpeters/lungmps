#! /bin/bash

#$ -cwd
#$ -q broad
#$ -P broad
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_vmem=24g
#$ -l h_rt=96:00:00
#$ -e liger_error.txt
#$ -o liger_out.txt
#$ -t 1-30
#$ -N output

cd <working_dir>
source /broad/software/scripts/useuse
reuse UGER
reuse -q R-3.5

Rscript <working_dir>/iNMF_script.R \
  $SGE_TASK_ID \
  output \
  20 \
  5
