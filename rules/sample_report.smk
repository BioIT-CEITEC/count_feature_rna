# REPORTING RULES
#

def merge_salmon_aln_input(wildcards):
    input = {}
    if config["salmon_align"]:
        input['salmon_align_tab'] = "qc_reports/{sample}/salmon/{sample}_aln/{sample}_aln.tsv"
    return input

rule merge_salmon_aln:
    input:  tsv=unpack(merge_salmon_aln_input),
    output: tsv="qc_reports/all_samples/salmon_aln.tsv"
    params: header=
    log:    "logs/all_samples/merge_salmon_aln.log"
    conda: "../wrappers/merge_salmon/env.yaml"
    script: "../wrappers/merge_salmon/script.py"

def merge_salmon_map_input(wildcards):
    input = {}
    if config["salmon_map"]:
        input['salmon_map_tab'] = "qc_reports/{sample}/salmon/{sample}_map/{sample}_map.tsv"
    return input

rule merge_salmon_map:
    input:  unpack(merge_salmon_map_input)
    output: tsv="qc_reports/all_samples/salmon_map.tsv"
    log:    "logs/all_samples/merge_salmon_map.log"
    conda: "../wrappers/merge_salmon/env.yaml"
    script: "../wrappers/merge_salmon/script.py"

def multiqc_report_input_files(wildcards):
    input = {}
    if config["feature_count"]:
        input['feature_count'] = "qc_reports/{sample}/feature_count/{sample}.feature_count.tsv"
    if config["RSEM"]:
        input['RSEM'] = "qc_reports/{sample}/RSEM/{sample}.genes.results"
    if config["salmon_align"]:
        input['salmon_align'] = "qc_reports/{sample}/salmon/{sample}_aln/{sample}.salmon_aln.sf",
        input['salmon_align_tab'] = "qc_reports/{sample}/salmon/{sample}_aln/{sample}_aln.tsv"
    if config["salmon_map"]:
        input['salmon_map'] = "qc_reports/{sample}/salmon/{sample}_map/{sample}.salmon_map.sf",
        input['salmon_map_tab'] = "qc_reports/{sample}/salmon/{sample}_map/{sample}_map.tsv"
    if config["kallisto"]:
        input['kallisto_h5'] = "qc_reports/{sample}/kallisto/{sample}.kallisto.h5",
        input['kallisto_tsv'] = "qc_reports/{sample}/kallisto/{sample}.kallisto.tsv"
    return input

def multiqc_report_input_paths(wildcards):
    input = {}
    if config["feature_count"]:
        input['feature_count'] = "qc_reports/*/feature_count/*"
    if config["RSEM"]:
        input['RSEM'] = "qc_reports/*/RSEM/*"
    if config["salmon_align"]:
        input['salmon_align'] = "qc_reports/*/salmon/*_aln/*"
        input['salmon_align_tab'] = "qc_reports/all_samples/salmon_aln.tsv"
    if config["salmon_map"]:
        input['salmon_map'] = "qc_reports/*/salmon/*_map/*"
    if config["kallisto"]:
        input['kallisto_h5'] = "qc_reports/*/kallisto/*"
    input_paths = " ".join(input.values())
    return input_paths

rule multiqc_report:
    input:  files=unpack(multiqc_report_input),
            paths=unpack(multiqc_report_paths)
    output: html="qc_reports/all_samples/multiqc_features.html"
    log:    "logs/all_samples/multiqc_features.log"
    params: multiqc_config = workflow.basedir+"/wrappers/multiqc_report/multiqc_config.txt",
            multiqc_path = "qc_reports/all_samples/"
    conda: "../wrappers/count_features_multiqc/env.yaml"
    script: "../wrappers/count_features_multiqc/script.py"


