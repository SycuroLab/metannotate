#!/bin/bash

#SBATCH --partition=cpu2019
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=5G
#SBATCH --error=metannotate_sbatch_run.%J.err
#SBATCH --output=metannotate_sbatch_run.%J.out

log_dir="$(pwd)"
log_file="logs/metannotate-analysis.log.txt"
num_jobs=10

echo "started at: `date`"

snakemake --cluster-config cluster.json --cluster 'sbatch --partition={cluster.partition} --cpus-per-task={cluster.cpus-per-task} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --time={cluster.time} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}' --jobs $num_jobs --use-conda &> $log_dir/$log_file

echo "finished with exit code $? at: `date`"

