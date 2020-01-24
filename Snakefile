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
        "output/humann2_genefamilies.tsv",
        "output/humann2_pathabundance.tsv",
        "output/humann2_pathcoverage.tsv"

rule humann2:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
        m = "data/merged/output/{sample}_merged.fastq",
        genefam = "output/{sample}_genefamilies.tsv",
        pathcov = "output/{sample}_pathcoverage.tsv",
        pathabun = "output/{sample}_pathabundance.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            """
            cat {input.r1} {input.r2} > data/merged/{output.m}
            humann2 --input {output.m} --output output/
            """

rule normalize:
    input:
        genefam = "output/{sample}_genefamilies.tsv",
        pathabun = "output/{sample}_pathabundance.tsv"
    output:
        genefam = "output/{sample}_genefamilies_norm.tsv",
        pathabun = "output/{sample}_pathabundance_norm.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            """
            humann2_renorm_table --input {input.genefam} --output {output.genefam} --units relab
            humann2_renorm_table --input {input.pathabun} --output {output.pathabund} --units relab
            """

rule merge:
    input:
        genefam = expand("output/{sample}_genefamilies_norm.tsv", sample=SAMPLES),
        pathcov = expand("output/{sample}_pathcoverage.tsv", sample=SAMPLES),
        pathabun = expand("output/{sample}_pathabundance_norm.tsv", sample=SAMPLES)
    output:
        genefam = "results/humann2_genefamilies.tsv",
        pathcov = "results/humann2_pathcoverage.tsv",
        pathabun = "results/humann2_pathabundance.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            """
            humann2_join_tables --input output/ --output {output.genefam} --filename genefamilies_norm
            humann2_join_tables --input output/ --output {output.pathcov} --filename pathcoverage
            humann2_join_tables --input output/ --output {output.pathabun} --filename pathabundance_norm
            """
