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
        config["output_dir"] + "/final_results_relab/merged_genefamilies_uniref_renamed_relab_stratified.tsv",
        config["output_dir"] + "/final_results_relab/merged_genefamilies_ko_renamed_relab_stratified.tsv",
        config["output_dir"] + "/final_results_relab/merged_genefamilies_rxn_renamed_relab_stratified.tsv",
        config["output_dir"] + "/final_results_relab/merged_genefamilies_eggnog_renamed_relab_stratified.tsv",
        config["output_dir"] + "/final_results_relab/merged_pathabundance_relab_stratified.tsv", 
        config["output_dir"] + "/final_results_relab/merged_genefamilies_uniref_renamed_relab_unstratified.tsv",
        config["output_dir"] + "/final_results_relab/merged_genefamilies_ko_renamed_relab_unstratified.tsv",
        config["output_dir"] + "/final_results_relab/merged_genefamilies_rxn_renamed_relab_unstratified.tsv",
        config["output_dir"] + "/final_results_relab/merged_genefamilies_eggnog_renamed_relab_unstratified.tsv",
        config["output_dir"] + "/final_results_relab/merged_pathabundance_relab_unstratified.tsv", 
        config["output_dir"] + "/final_results_raw/merged_genefamilies_uniref_renamed_stratified.tsv",
        config["output_dir"] + "/final_results_raw/merged_genefamilies_ko_renamed_stratified.tsv",
        config["output_dir"] + "/final_results_raw/merged_genefamilies_rxn_renamed_stratified.tsv",
        config["output_dir"] + "/final_results_raw/merged_genefamilies_eggnog_renamed_stratified.tsv",
        config["output_dir"] + "/final_results_raw/merged_pathabundance_stratified.tsv", 
        config["output_dir"] + "/final_results_raw/merged_genefamilies_uniref_renamed_unstratified.tsv",
        config["output_dir"] + "/final_results_raw/merged_genefamilies_ko_renamed_unstratified.tsv",
        config["output_dir"] + "/final_results_raw/merged_genefamilies_rxn_renamed_unstratified.tsv",
        config["output_dir"] + "/final_results_raw/merged_genefamilies_eggnog_renamed_unstratified.tsv",
        config["output_dir"] + "/final_results_raw/merged_pathabundance_unstratified.tsv", 
        config["output_dir"] + "/final_results_cpm/merged_genefamilies_uniref_renamed_cpm_stratified.tsv",
        config["output_dir"] + "/final_results_cpm/merged_genefamilies_ko_renamed_cpm_stratified.tsv",
        config["output_dir"] + "/final_results_cpm/merged_genefamilies_rxn_renamed_cpm_stratified.tsv",
        config["output_dir"] + "/final_results_cpm/merged_genefamilies_eggnog_renamed_cpm_stratified.tsv",
        config["output_dir"] + "/final_results_cpm/merged_pathabundance_cpm_stratified.tsv", 
        config["output_dir"] + "/final_results_cpm/merged_genefamilies_uniref_renamed_cpm_unstratified.tsv",
        config["output_dir"] + "/final_results_cpm/merged_genefamilies_ko_renamed_cpm_unstratified.tsv",
        config["output_dir"] + "/final_results_cpm/merged_genefamilies_rxn_renamed_cpm_unstratified.tsv",
        config["output_dir"] + "/final_results_cpm/merged_genefamilies_eggnog_renamed_cpm_unstratified.tsv",
        config["output_dir"] + "/final_results_cpm/merged_pathabundance_cpm_unstratified.tsv" 

# Using this rule if you are merging two lanes otherwise comment it out and use the other merge_reads rule.
rule merge_reads:
    input:
        r11 = config["path"]+"{sample}_L001"+config["for"],
        r12 = config["path"]+"{sample}_L001"+config["rev"],
        r21 = config["path"]+"{sample}_L002"+config["for"],
        r22 = config["path"]+"{sample}_L002"+config["rev"]
    output:
        r1 = config["output_dir"]+"/merged_data/{sample}"+config["for"],
        r2 = config["output_dir"]+"/merged_data/{sample}"+config["rev"],
        merged_fastq = config["output_dir"]+"/merged_data/{sample}.fastq"
    shell:
        "cat {input.r11} {input.r21} > {output.r1}; "
        "cat {input.r12} {input.r22} > {output.r2}; "
        "cat {output.r1} {output.r2} > {output.merged_fastq}; "

#rule merge_reads:
#    input:
#        r1 = config["path"]+"{sample}"+config["for"],
#        r2 = config["path"]+"{sample}"+config["rev"]
#    output:
#        merged_fastq = config["output_dir"]+"/merged_data/{sample}.fastq"
#    shell:
#            "cat {input.r1} {input.r2} > {output.merged_fastq}"

rule humann3:
    input:
        config["output_dir"]+"/merged_data/{sample}.fastq" if config["paired"] else config["path"]+"{sample}"+config["suff"]
    output:
        genefam = config["output_dir"] + "/raw/{sample}_genefamilies.tsv",
        pathcov = config["output_dir"] + "/raw/{sample}_pathcoverage.tsv",
        pathabun = config["output_dir"] + "/raw/{sample}_pathabundance.tsv"
    params:
        num_threads = config["num_threads"],
        db1 = "--bowtie2db " + config["output_dir"] + "/raw/{sample}_database",
        db2 = config["output_dir"] + "/databases/{sample}_database",
        metaphlan= config["metaphlan_results_path"],
        s = "{sample}",
        output=config["output_dir"] + "/raw"

    conda: "humann3"
    shell:
            """
            humann3 --input {input} --threads {params.num_threads} --output {params.output} --output-basename {params.s} --nucleotide-database {config[nuc_db]} --protein-database {config[prot_db]} --taxonomic-profile {params.metaphlan}
# This has been taken out becaue we can use the metaphlan outputs from the metaphlan run and don't need to re-do them.
#--metaphlan-options="{params.db1}"
            rm -rf {params.db2}
            """

rule regroup:
     input: genefam=config["output_dir"] + "/raw/{sample}_genefamilies.tsv"
     output:genefam_eggnog=config["output_dir"] + "/regrouped/{sample}_genefamilies_eggnog.tsv",
            genefam_ko=config["output_dir"] + "/regrouped/{sample}_genefamilies_ko.tsv",
            genefam_rxn=config["output_dir"] + "/regrouped/{sample}_genefamilies_rxn.tsv"
     params: eggnog=config["eggnog_uniref90"],
             ko=config["ko_uniref90"]
     conda: "humann3"
     shell:
       """ humann_regroup_table --input {input.genefam}     --output {output.genefam_ko}  -c {params.ko} ; 
           humann_regroup_table --input {input.genefam}     --output {output.genefam_eggnog}  -c {params.eggnog};
           humann_regroup_table --input {input.genefam}     --output {output.genefam_rxn}  -g uniref90_rxn;"""


rule rename:
      input:genefam=config["output_dir"] + "/raw/{sample}_genefamilies.tsv",
            genefam_eggnog=config["output_dir"] + "/regrouped/{sample}_genefamilies_eggnog.tsv",
            genefam_ko=config["output_dir"] + "/regrouped/{sample}_genefamilies_ko.tsv",
            genefam_rxn=config["output_dir"] + "/regrouped/{sample}_genefamilies_rxn.tsv"

      output: 
            genefam_name=config["output_dir"] + "/renamed/uniref/{sample}_genefamilies_renamed.tsv",
            genefam_ko_name=config["output_dir"] + "/renamed/ko/{sample}_genefamilies_ko_renamed.tsv",
            genefam_rxn_name=config["output_dir"] + "/renamed/rxn/{sample}_genefamilies_rxn_renamed.tsv",
            genefam_eggnog_name=config["output_dir"] + "/renamed/eggnog/{sample}_genefamilies_eggnog_renamed.tsv"
      conda: "humann3"
      params: uniref=config["uniref90_name"], 
              eggnog=config["eggnog_name"],
              ko=config["ko_name"]
      shell:
       """ 
	humann_rename_table -i {input.genefam} -c {params.uniref} -o {output.genefam_name};
	humann_rename_table -i {input.genefam_ko} -c {params.ko} -o {output.genefam_ko_name};
	humann_rename_table -i {input.genefam_eggnog} -c {params.eggnog} -o {output.genefam_eggnog_name};
	humann_rename_table -i {input.genefam_rxn} -n metacyc-rxn -o {output.genefam_rxn_name} """
     


rule merge_output:
    input:
            genefam_name=expand(config["output_dir"] + "/renamed/uniref/{sample}_genefamilies_renamed.tsv",sample=SAMPLES),
            genefam_eggnog_name=expand(config["output_dir"] + "/renamed/eggnog/{sample}_genefamilies_eggnog_renamed.tsv",sample=SAMPLES),
            genefam_ko_name=expand(config["output_dir"] + "/renamed/ko/{sample}_genefamilies_ko_renamed.tsv",sample=SAMPLES),
            genefam_rxn_name=expand(config["output_dir"] + "/renamed/rxn/{sample}_genefamilies_rxn_renamed.tsv",sample=SAMPLES),
            pathabun = expand(config["output_dir"] + "/raw/{sample}_pathabundance.tsv", sample=SAMPLES),
            pathcov = expand(config["output_dir"] + "/raw/{sample}_pathcoverage.tsv", sample=SAMPLES)
    params:
        renamed_uniref_dir=config["output_dir"] + "/renamed/uniref/",
        renamed_ko_dir=config["output_dir"] + "/renamed/ko/",
        renamed_rxn_dir=config["output_dir"] + "/renamed/rxn/",
        renamed_eggnog_dir=config["output_dir"] + "/renamed/eggnog/",
        renamed_raw_dir=config["output_dir"] + "/raw/",
    output:
        genefam_uniref = config["output_dir"] + "/merged/merged_genefamilies_uniref_renamed.tsv",
        genefam_eggnog = config["output_dir"] + "/merged/merged_genefamilies_eggnog_renamed.tsv",
        genefam_ko = config["output_dir"] + "/merged/merged_genefamilies_ko_renamed.tsv",
        genefam_rxn = config["output_dir"] + "/merged/merged_genefamilies_rxn_renamed.tsv",       
        pathabun = config["output_dir"] + "/merged/merged_pathabundance.tsv",
        pathcov = config["output_dir"] + "/merged/merged_pathcoverage.tsv"
    conda: "humann3"
    shell:
         """
            humann_join_tables --input {params.renamed_uniref_dir} --output {output.genefam_uniref} --file_name genefamilies; 
            humann_join_tables --input {params.renamed_ko_dir} --output {output.genefam_ko} --file_name genefamilies; 
            humann_join_tables --input {params.renamed_rxn_dir} --output {output.genefam_rxn} --file_name genefamilies; 
            humann_join_tables --input {params.renamed_eggnog_dir} --output {output.genefam_eggnog} --file_name genefamilies; 
            humann_join_tables --input {params.renamed_raw_dir} --output {output.pathcov} --file_name pathcoverage; 
            humann_join_tables --input {params.renamed_raw_dir}  --output {output.pathabun} --file_name _pathabundance; """


rule relab:
    input:
        genefam_uniref = config["output_dir"] + "/merged/merged_genefamilies_uniref_renamed.tsv",
        genefam_ko = config["output_dir"] + "/merged/merged_genefamilies_ko_renamed.tsv",
        genefam_rxn = config["output_dir"] + "/merged/merged_genefamilies_rxn_renamed.tsv",
        genefam_eggnog = config["output_dir"] + "/merged/merged_genefamilies_eggnog_renamed.tsv",
        pathabun = config["output_dir"] + "/merged/merged_pathabundance.tsv"
    output:
        genefam_uniref_relab = config["output_dir"] + "/merged/merged_genefamilies_uniref_renamed_relab.tsv",
        genefam_ko_relab = config["output_dir"] + "/merged/merged_genefamilies_ko_renamed_relab.tsv",
        genefam_rxn_relab = config["output_dir"] + "/merged/merged_genefamilies_rxn_renamed_relab.tsv",
        genefam_eggnog_relab = config["output_dir"] + "/merged/merged_genefamilies_eggnog_renamed_relab.tsv",
        pathabun_relab = config["output_dir"] + "/merged/merged_pathabundance_relab.tsv" 
    conda: "humann3"
    shell:
         """
         humann_renorm_table --input {input.genefam_uniref} --output {output.genefam_uniref_relab} -s n --units relab;
         humann_renorm_table --input {input.genefam_ko} --output {output.genefam_ko_relab} -s n --units relab;
         humann_renorm_table --input {input.genefam_rxn} --output {output.genefam_rxn_relab} -s n --units relab;
         humann_renorm_table --input {input.genefam_eggnog} --output {output.genefam_eggnog_relab} -s n --units relab;
         humann_renorm_table --input {input.pathabun} --output {output.pathabun_relab} -s n --units relab; """


rule cpm:
    input:
        genefam_uniref = config["output_dir"] + "/merged/merged_genefamilies_uniref_renamed.tsv",
        genefam_ko = config["output_dir"] + "/merged/merged_genefamilies_ko_renamed.tsv",
        genefam_rxn = config["output_dir"] + "/merged/merged_genefamilies_rxn_renamed.tsv",
        genefam_eggnog = config["output_dir"] + "/merged/merged_genefamilies_eggnog_renamed.tsv",
        pathabun = config["output_dir"] + "/merged/merged_pathabundance.tsv"
    output:
        genefam_uniref_cpm = config["output_dir"] + "/merged/merged_genefamilies_uniref_renamed_cpm.tsv",
        genefam_ko_cpm = config["output_dir"] + "/merged/merged_genefamilies_ko_renamed_cpm.tsv",
        genefam_rxn_cpm = config["output_dir"] + "/merged/merged_genefamilies_rxn_renamed_cpm.tsv",
        genefam_eggnog_cpm = config["output_dir"] + "/merged/merged_genefamilies_eggnog_renamed_cpm.tsv",
        pathabun_cpm = config["output_dir"] + "/merged/merged_pathabundance_cpm.tsv" 
    conda: "humann3"
    shell:
         """
         humann_renorm_table --input {input.genefam_uniref} --output {output.genefam_uniref_cpm} -s n --units cpm;
         humann_renorm_table --input {input.genefam_ko} --output {output.genefam_ko_cpm} -s n --units cpm;
         humann_renorm_table --input {input.genefam_rxn} --output {output.genefam_rxn_cpm} -s n --units cpm;
         humann_renorm_table --input {input.genefam_eggnog} --output {output.genefam_eggnog_cpm} -s n --units cpm;
         humann_renorm_table --input {input.pathabun} --output {output.pathabun_cpm} -s n --units cpm; """


rule final_results_raw:
    input: 
        genefam_uniref_raw = config["output_dir"] + "/merged/merged_genefamilies_uniref_renamed.tsv",
        genefam_ko_raw = config["output_dir"] + "/merged/merged_genefamilies_ko_renamed.tsv",
        genefam_rxn_raw = config["output_dir"] + "/merged/merged_genefamilies_rxn_renamed.tsv",
        genefam_eggnog_raw = config["output_dir"] + "/merged/merged_genefamilies_eggnog_renamed.tsv",
        pathabun_raw = config["output_dir"] + "/merged/merged_pathabundance.tsv" 
    params: 
        output_dir=config["output_dir"] + "/final_results_raw"
    output:
        genefam_uniref_raw_s = config["output_dir"] + "/final_results_raw/merged_genefamilies_uniref_renamed_stratified.tsv",
        genefam_ko_raw_s = config["output_dir"] + "/final_results_raw/merged_genefamilies_ko_renamed_stratified.tsv",
        genefam_rxn_raw_s = config["output_dir"] + "/final_results_raw/merged_genefamilies_rxn_renamed_stratified.tsv",
        genefam_eggnog_raw_s = config["output_dir"] + "/final_results_raw/merged_genefamilies_eggnog_renamed_stratified.tsv",
        pathabun_raw_s = config["output_dir"] + "/final_results_raw/merged_pathabundance_stratified.tsv", 
        genefam_uniref_raw_us = config["output_dir"] + "/final_results_raw/merged_genefamilies_uniref_renamed_unstratified.tsv",
        genefam_ko_raw_us = config["output_dir"] + "/final_results_raw/merged_genefamilies_ko_renamed_unstratified.tsv",
        genefam_rxn_raw_us = config["output_dir"] + "/final_results_raw/merged_genefamilies_rxn_renamed_unstratified.tsv",
        genefam_eggnog_raw_us = config["output_dir"] + "/final_results_raw/merged_genefamilies_eggnog_renamed_unstratified.tsv",
        pathabun_raw_us = config["output_dir"] + "/final_results_raw/merged_pathabundance_unstratified.tsv" 

    conda: "humann3"
    shell:
        """ 
            humann_split_stratified_table --input {input.genefam_uniref_raw} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_ko_raw} --output  {params.output_dir};
            humann_split_stratified_table --input {input.genefam_rxn_raw} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_eggnog_raw} --output {params.output_dir};
            humann_split_stratified_table --input {input.pathabun_raw} --output {params.output_dir};

         """

rule final_results_relab:
    input: 
        genefam_uniref_relab = config["output_dir"] + "/merged/merged_genefamilies_uniref_renamed_relab.tsv",
        genefam_ko_relab = config["output_dir"] + "/merged/merged_genefamilies_ko_renamed_relab.tsv",
        genefam_rxn_relab = config["output_dir"] + "/merged/merged_genefamilies_rxn_renamed_relab.tsv",
        genefam_eggnog_relab = config["output_dir"] + "/merged/merged_genefamilies_eggnog_renamed_relab.tsv",
        pathabun_relab = config["output_dir"] + "/merged/merged_pathabundance_relab.tsv" 
    
    output:
        genefam_uniref_relab_s = config["output_dir"] + "/final_results_relab/merged_genefamilies_uniref_renamed_relab_stratified.tsv",
        genefam_ko_relab_s = config["output_dir"] + "/final_results_relab/merged_genefamilies_ko_renamed_relab_stratified.tsv",
        genefam_rxn_relab_s = config["output_dir"] + "/final_results_relab/merged_genefamilies_rxn_renamed_relab_stratified.tsv",
        genefam_eggnog_relab_s = config["output_dir"] + "/final_results_relab/merged_genefamilies_eggnog_renamed_relab_stratified.tsv",
        pathabun_relab_s = config["output_dir"] + "/final_results_relab/merged_pathabundance_relab_stratified.tsv", 
        genefam_uniref_relab_us = config["output_dir"] + "/final_results_relab/merged_genefamilies_uniref_renamed_relab_unstratified.tsv",
        genefam_ko_relab_us = config["output_dir"] + "/final_results_relab/merged_genefamilies_ko_renamed_relab_unstratified.tsv",
        genefam_rxn_relab_us = config["output_dir"] + "/final_results_relab/merged_genefamilies_rxn_renamed_relab_unstratified.tsv",
        genefam_eggnog_relab_us = config["output_dir"] + "/final_results_relab/merged_genefamilies_eggnog_renamed_relab_unstratified.tsv",
        pathabun_relab_us = config["output_dir"] + "/final_results_relab/merged_pathabundance_relab_unstratified.tsv" 

    params: 
        output_dir=config["output_dir"] + "/final_results_relab"

    conda: "humann3"
    shell:
        """ 
            humann_split_stratified_table --input {input.genefam_uniref_relab} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_ko_relab} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_rxn_relab} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_eggnog_relab} --output {params.output_dir};
            humann_split_stratified_table --input {input.pathabun_relab} --output {params.output_dir};
 
         """


rule final_results_cpm:
    input: 
        genefam_uniref_cpm = config["output_dir"] + "/merged/merged_genefamilies_uniref_renamed_cpm.tsv",
        genefam_ko_cpm = config["output_dir"] + "/merged/merged_genefamilies_ko_renamed_cpm.tsv",
        genefam_rxn_cpm = config["output_dir"] + "/merged/merged_genefamilies_rxn_renamed_cpm.tsv",
        genefam_eggnog_cpm = config["output_dir"] + "/merged/merged_genefamilies_eggnog_renamed_cpm.tsv",
        pathabun_cpm = config["output_dir"] + "/merged/merged_pathabundance_cpm.tsv" 
    
    output:
        genefam_uniref_cpm_s = config["output_dir"] + "/final_results_cpm/merged_genefamilies_uniref_renamed_cpm_stratified.tsv",
        genefam_ko_cpm_s = config["output_dir"] + "/final_results_cpm/merged_genefamilies_ko_renamed_cpm_stratified.tsv",
        genefam_rxn_cpm_s = config["output_dir"] + "/final_results_cpm/merged_genefamilies_rxn_renamed_cpm_stratified.tsv",
        genefam_eggnog_cpm_s = config["output_dir"] + "/final_results_cpm/merged_genefamilies_eggnog_renamed_cpm_stratified.tsv",
        pathabun_cpm_s = config["output_dir"] + "/final_results_cpm/merged_pathabundance_cpm_stratified.tsv", 
        genefam_uniref_cpm_us = config["output_dir"] + "/final_results_cpm/merged_genefamilies_uniref_renamed_cpm_unstratified.tsv",
        genefam_ko_cpm_us = config["output_dir"] + "/final_results_cpm/merged_genefamilies_ko_renamed_cpm_unstratified.tsv",
        genefam_rxn_cpm_us = config["output_dir"] + "/final_results_cpm/merged_genefamilies_rxn_renamed_cpm_unstratified.tsv",
        genefam_eggnog_cpm_us = config["output_dir"] + "/final_results_cpm/merged_genefamilies_eggnog_renamed_cpm_unstratified.tsv",
        pathabun_cpm_us = config["output_dir"] + "/final_results_cpm/merged_pathabundance_cpm_unstratified.tsv" 
    params: 
        output_dir=config["output_dir"] + "/final_results_cpm"

    conda: "humann3"
    shell:
        """ 
            humann_split_stratified_table --input {input.genefam_uniref_cpm} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_ko_cpm} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_rxn_cpm} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_eggnog_cpm} --output {params.output_dir};
            humann_split_stratified_table --input {input.pathabun_cpm} --output {params.output_dir};
         """
