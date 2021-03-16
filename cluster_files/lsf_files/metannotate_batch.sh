#! /usr/bin/env bash

#BSUB -J metannotate_batch_run
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 60:00
#BSUB -o metannotate_batch_run_%J.out
#BSUB -e metannotate_batch_run_%J.err

log_dir="$(pwd)"
log_file="logs/metannotate.log.txt"
num_jobs=20

snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs $num_jobs --use-conda &> $log_dir/$log_file

echo "finished with exit code $? at: `date`"

