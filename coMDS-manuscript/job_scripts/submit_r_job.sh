#!/bin/bash

#$ -N job_name
#$ -pe smp 20
#$ -o logs/$JOB_NAME.txt

module load R

cd ../
Rscript scripts/${1}.R "${@:2}"
