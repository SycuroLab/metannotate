# metannotate

Snakemake pipeline for functionally annotating microbial communities from metagenomic shotgun sequencing data using [HUMAnN3](http://huttenhower.sph.harvard.edu/humann3).

HUMAnN3 uses a tiered approach - it first maps reads to clade-specific marker genes to identify species present in samples, then maps reads to functionally annotated pangenomes of identified species, and finally aligns unclassified reads to a protein database using a translated search (DIAMOND). 

## Overview

Input: 

* Cleaned and filtered reads from shotgun metagenome sequencing. Can be paired or unpaired.

Output: 

This pipeline produces 7 files, which can be found in the `results` directory upon completion of all steps.

* `merged_genefamilies.tsv`, a table of gene family abundances.
* `merged_genefamilies_relab.tsv`, a table of gene family abundances, normalized using relative abundance.
* `merged_genefamilies_cpm.tsv`, a table of gene family abundances, normalized using copies per million.
* `merged_pathabundance.tsv`, a table of metabolic pathway abundances.
* `merged_pathabundance_relab.tsv`, a table of metabolic pathway abundances, normalized by relative abundance.
* `merged_pathabundance_cpm.tsv`, a table of metabolic pathway  abundances, normalized by copies per million.
* `merged_pathcoverage.tsv`, a table of the presence (1) and absence (0) of pathways in a community, independent of their quantitative abundance.

## Pipeline summary

<!---Insert DAG here if needed? -->

<!--- % ![Rulegraph](./metaphlan_files/rulegraph.png) -->

### Steps

0) If reads are paired, merge forward and reverse reads.

1) Run HUMAnN2 on samples. 

2) Merge output files.

3) Normalize output.

## Installation

To use this pipeline, navigate to your project directory and clone this repository into that directory using the following command:

```
git clone https://github.com/SycuroLab/metannotate.git metannotate
```

Note: you need to have **conda** and **snakemake** installed in order to run this. To install conda, see the instructions [here](https://github.com/ucvm/synergy/wiki). 

To install snakemake using conda, run the following line:

```
conda install -c bioconda -c conda-forge snakemake
```

See the snakemake installation [webpage](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for further details.

## Config file

All the parameters required to run this pipeline are specified in a config file, written in yaml. See/modify the provided example file with your custom parameters, called `config.yaml`. This is the only file that should be modified before running the pipeline. Make sure to follow the syntax in the example file in terms of when to use quotations around parameters.

## HUMAnN3 Databases

Two databases are needed to run this pipeline - a nucleotide database (ChocoPhlAn) and a protein database (UniRef).

The HUMAnN3 Databases were downloaded and installed previously for ease of use.

Location: 
`arc.ucalgary.ca`

Directory Paths:

The nuclotide database path:
nuc_db `/bulk/IMCshared_bulk/shared/dbs/humann_dbs/chocophlan`

The protein database path:
prot_db `/bulk/IMCshared_bulk/shared/dbs/humann_dbs/uniref`

If there are newer versions of the databases that you want to use for your project you can download the newer versions;

To download them, follow the instructions [here](https://huttenhower.sph.harvard.edu/humann/). 

Place the paths of each database in its corresponding parameter nuc_db for the nucleotide database and prot_db for the protein database in the `config.yaml` file.

## Data and list of files

Specify the full path to the directory that contains your data files in the config file. You also need to have a list of sample names which contains the names of the samples to run the pipeline on, one sample per line. You can run this pipeline on any number or subset of your samples. Sample names should include everything up to the suffix of the file names of the fastq files. Specify the path and name of your list in the config file.

## Description of parameters
| Parameter | Description | Example |
| -------------- | --------------- | ------------ |
| list_files | Full path and name of your sample list. | `"/export/home/hramay/projects/Arrieta/PC1000/antibiotic_babies/analysis/metqc/list_files.txt"` |
| path | Location of input files. | `"/export/home/hramay/projects/Arrieta/PC1000/antibiotic_babies/analysis/metqc/output/bmtagger/"` |
| metaphlan_results_path | Location of processed metaphlan files | `"/export/home/hramay/projects/Arrieta/PC1000/antibiotic_babies/analysis/metaphlan/output/metaphlan/"` |
| num_threads | Number of threads to use for humann3. | `8` |
| paired | Are the input reads paired? | `TRUE` |
| for | If paired, suffix of forward reads. | `"_filtered_1.fastq"` |
| rev | If paired, suffix of reverse reads. | `"_filtered_2.fastq"` |
| suff | If unpaired, suffix of reads. | `"_bmt_merged_ELC_trimmed_filtered.fastq"`
| nuc_db | Location of nucleotide database. | `"/bulk/IMCbinf_bulk/hramay/projects/databases/humann_dbs/chocophlan"` |
| prot_db | Location of protein database. | `"/bulk/IMCbinf_bulk/hramay/projects/databases/humann_dbs/uniref"` |

## Running the pipeline on ARC (SLURM cluster)

Test the pipeline by running `snakemake -np`. This command prints out the commands to be run without actually running them. 

To run the pipeline on the ARC compute cluster, enter the following command from the project directory:

```
sbatch < metannotate_sbatch.sh
```

The above command submits jobs to ARC, one for each sample and step of the metannotate pipeline.

Note: the file `cluster.json` contains the parameters for the SLURM job submission system that ARC uses. In most cases, this file should not be modified. Use the `cluster.json` file in the `cluster_files/slurm_files/` folder. 

The ARC Cluster Guide can be found here:
https://rcs.ucalgary.ca/index.php/ARC_Cluster_Guide

The General Guidelines and Policies can be found here:
https://rcs.ucalgary.ca/index.php/General_Cluster_Guidelines_and_Policies


## Running the pipeline on Synergy (LSF cluster)

Test the pipeline by running `snakemake -np`. This command prints out the commands to be run without actually running them. 

To run the pipeline on the Synergy compute cluster, enter the following command from the project directory:

```
snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 500 --use-conda
```
The above command submits jobs to Synergy, one for each sample and step of the metaphlan pipeline. Note: the file `cluster.json` in the `cluster_files/lsf_files/` folder contains the parameters for the LSF job submission system that Synergy uses. In most cases, this file should not be modified.

## Results and log files

Snakemake will create a directory for the results of the pipeline as well as a directory for log files. Log files of each step of the pipeline will be written to the `logs` directory. Intermediate output files can be found in the `output` directory. 




