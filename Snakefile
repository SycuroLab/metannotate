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
        "results/merged_genefamilies_relab.tsv",
        "results/merged_pathabundance_relab.tsv",
        "results/merged_pathcoverage.tsv",
        "results/merged_genefamilies_cpm.tsv",
        "results/merged_pathabundance_cpm.tsv"

rule merge_reads:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
        "data/merged/{sample}.fastq"
    shell:
            "cat {input.r1} {input.r2} > {output.m}"

rule humann2:
    input:
        "data/merged/{sample}.fastq" if config["paired"] else config["path"]+"{sample}"+config["suff"]
    output:
        genefam = "output/{sample}_genefamilies.tsv",
        pathcov = "output/{sample}_pathcoverage.tsv",
        pathabun = "output/{sample}_pathabundance.tsv"
    params:
        db1 = "--bowtie2db output/databases/{sample}_database",
        db2 = "output/database/{sample}_database",
        s = "{sample}"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            """
            humann2 --input {input} --threads 16 --output output --output-basename {params.s} --nucleotide-database {config[nuc_db]} --protein-database {config[prot_db]} --metaphlan-options="{params.db1}"
            rm -rf {params.db2}
            """

rule merge_output:
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
        gfrelab = "results/merged_genefamilies_relab.tsv",
        parelab = "results/merged_pathabundance_relab.tsv",
        gfcpm = "results/merged_genefamilies_cpm.tsv",
        pacpm = "results/merged_pathabundance_cpm.tsv"
    conda: "utils/envs/humann2_env.yaml"
    shell:
            "humann2_renorm_table --input {input.genefam} --output {output.gfrelab} --units relab; humann2_renorm_table --input {input.pathabun} --output {output.parelab} --units relab; humann2_renorm_table --input {input.genefam} --output {output.gfcpm}; humann2_renorm_table --input {input.pathabun} --output {output.pacpm}"
