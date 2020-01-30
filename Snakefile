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
        "results/merged_genefamilies_norm.tsv",
        "results/merged_pathabundance_norm.tsv",
        "results/merged_pathcoverage.tsv"

rule humann2:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
        m = "data/merged/{sample}.fastq",
        genefam = "output/{sample}_genefamilies.tsv",
        pathcov = "output/{sample}_pathcoverage.tsv",
        pathabun = "output/{sample}_pathabundance.tsv"
    params:
        db1 = "--bowtie2db output/databases/{sample}_database",
        db2 = "output/database/{sample}_database"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            """
            cat {input.r1} {input.r2} > {output.m}
            humann2 --input {output.m} --threads 16 --output output --nucleotide-database /home/aschick/refs/humann2/chocophlan --protein-database /home/aschick/refs/humann2/uniref --metaphlan-options="{params.db1}"
            rm -rf {params.db2}
            """

rule merge:
    input:
        genefam = expand("output/{sample}_genefamilies.tsv", sample=SAMPLES),
        pathcov = expand("output/{sample}_pathcoverage.tsv", sample=SAMPLES),
        pathabun = expand("output/{sample}_pathabundance.tsv", sample=SAMPLES)
    output:
        genefam = "results/merged_genefamilies.tsv",
        pathcov = "results/merged_pathcoverage.tsv",
        pathabun = "results/merged_pathabundance.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            "humann2_join_tables --input output/ --output {output.genefam} --file_name genefamilies; humann2_join_tables --input output/ --output {output.pathcov} --file_name pathcoverage; humann2_join_tables --input output/ --output {output.pathabun} --file_name pathabundance"

rule normalize:
    input:
        genefam = "results/merged_genefamilies.tsv",
        pathabun = "results/merged_pathabundance.tsv"
    output:
        genefam = "results/merged_genefamilies_norm.tsv",
        pathabun = "results/merged_pathabundance_norm.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            "humann2_renorm_table --input {input.genefam} --output {output.genefam} --units relab; humann2_renorm_table --input {input.pathabun} --output {output.pathabun} --units relab"
