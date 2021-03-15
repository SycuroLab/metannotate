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
        "output/merged_genefamilies_relab.tsv",
        "output/merged_pathabundance_relab.tsv",
        "output/merged_pathcoverage.tsv",
        "output/merged_genefamilies_cpm.tsv",
        "output/merged_pathabundance_cpm.tsv"

rule merge_reads:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
        "data/merged/{sample}.fastq"
    shell:
            "cat {input.r1} {input.r2} > {output}"

rule humann3:
    input:
        "data/merged/{sample}.fastq" if config["paired"] else config["path"]+"{sample}"+config["suff"]
    output:
        genefam = "output/{sample}_genefamilies.tsv",
        pathcov = "output/{sample}_pathcoverage.tsv",
        pathabun = "output/{sample}_pathabundance.tsv"
    params:
        db1 = "--bowtie2db output/databases/{sample}_database",
        db2 = "output/database/{sample}_database",
        metaphlan= config["metaphlan_results_path"],
        s = "{sample}"
    conda: "utils/envs/humann3_env.yaml"
    shell:
            """
            humann3 --input {input} --threads 2 --output output --output-basename {params.s} --nucleotide-database {config[nuc_db]} --protein-database {config[prot_db]} --taxonomic-profile {params.metaphlan}
# This has been taken out becaue we can use the metaphlan outputs from the metaphlan run and don't need to re-do them.
#--metaphlan-options="{params.db1}"
            rm -rf {params.db2}
            """

 

rule merge_output:
    input:
        genefam = expand("output/{sample}_genefamilies.tsv", sample=SAMPLES),
        pathcov = expand("output/{sample}_pathcoverage.tsv", sample=SAMPLES),
        pathabun = expand("output/{sample}_pathabundance.tsv", sample=SAMPLES)
    output:
        genefam = "output/merged_genefamilies.tsv",
        pathcov = "output/merged_pathcoverage.tsv",
        pathabun = "output/merged_pathabundance.tsv"
    conda: "utils/envs/humann3_env.yaml"
    shell:
            "humann_join_tables --input output/ --output {output.genefam} --file_name genefamilies; humann_join_tables --input output/ --output {output.pathcov} --file_name pathcoverage; humann_join_tables --input output/ --output {output.pathabun} --file_name pathabundance"

rule normalize:
    input:
        genefam = "output/merged_genefamilies.tsv",
        pathabun = "output/merged_pathabundance.tsv"
    output:
        gfrelab = "output/merged_genefamilies_relab.tsv",
        parelab = "output/merged_pathabundance_relab.tsv",
        gfcpm = "output/merged_genefamilies_cpm.tsv",
        pacpm = "output/merged_pathabundance_cpm.tsv"
    conda: "utils/envs/humann3_env.yaml"
    shell:
            "humann_renorm_table --input {input.genefam} --output {output.gfrelab} --units relab; humann_renorm_table --input {input.pathabun} --output {output.parelab} --units relab; humann_renorm_table --input {input.genefam} --output {output.gfcpm}; humann_renorm_table --input {input.pathabun} --output {output.pacpm}"
