# *************************************
# * Snakefile for metannotate pipeline *
# *************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:
        "results/humann2_genefamilies.tsv",
        "results/humann2_pathabundance.tsv",
        "results/humann2_pathcoverage.tsv"

rule humann2:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
        m = "data/merged/{sample}_merged.fastq",
        genefam = "output/{sample}_merged_genefamilies.tsv",
        pathcov = "output/{sample}_merged_pathcoverage.tsv",
        pathabun = "output/{sample}_merged_pathabundance.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            """
            cat {input.r1} {input.r2} > {output.m}
            humann2 --input {output.m} --output output --threads 16 --nucleotide-database /home/aschick/refs/humann2/chocophlan --protein-database /home/aschick/refs/humann2/uniref --metaphlan-options="--bowtie2db /home/aschick/miniconda3/envs/humann2/bin/databases"
            """

rule normalize:
    input:
        genefam = "output/{sample}_merged_genefamilies.tsv",
        pathabun = "output/{sample}_merged_pathabundance.tsv"
    output:
        genefam = "output/{sample}_genefamilies_norm.tsv",
        pathabun = "output/{sample}_pathabundance_norm.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            """
            humann2_renorm_table --input {input.genefam} --output {output.genefam} --units relab
            humann2_renorm_table --input {input.pathabun} --output {output.pathabun} --units relab
            """

rule merge:
    input:
        genefam = expand("output/{sample}_genefamilies_norm.tsv", sample=SAMPLES),
        pathcov = expand("output/{sample}_merged_pathcoverage.tsv", sample=SAMPLES),
        pathabun = expand("output/{sample}_pathabundance_norm.tsv", sample=SAMPLES)
    output:
        genefam = "results/humann2_genefamilies.tsv",
        pathcov = "results/humann2_pathcoverage.tsv",
        pathabun = "results/humann2_pathabundance.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            """
            humann2_join_tables --input output/ --output {output.genefam} --file_name genefamilies_norm
            humann2_join_tables --input output/ --output {output.pathcov} --file_name pathcoverage
            humann2_join_tables --input output/ --output {output.pathabun} --file_name pathabundance_norm
            """
