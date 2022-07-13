# REPORTING RULES
#

def multiqc_report_input(wildcards):
    input = {}
    if config["feature_count"]:
        input['feature_count'] = "qc_reports/{sample}/feature_count/{sample}.feature_count.tsv"
    if config["RSEM"]:
        input['RSEM'] = "qc_reports/{sample}/RSEM/{sample}.genes.results"
    if config["salmon_align"]:
        input['salmon_align'] = "qc_reports/{sample}/salmon/{sample}_aln/{sample}.salmon_aln.sf"
    if config["salmon_map"]:
        input['salmon_map'] = "qc_reports/{sample}/salmon/{sample}_map/{sample}.salmon_map.sf"
    if config["kallisto"]:
        input['kallisto_h5'] = "qc_reports/{sample}/kallisto/{sample}.kallisto.h5",
        input['kallisto_tsv'] = "qc_reports/{sample}/kallisto/{sample}.kallisto.tsv"
        return input

rule multiqc_report:
    input:  unpack(multiqc_report_input)
    output: html="qc_reports/{sample}/multiqc_features.html"
    log:    "logs/{sample}/multiqc.log"
    params: multiqc_config = workflow.basedir+"/wrappers/multiqc_report/multiqc_config.txt",
            multiqc_path = "qc_reports/{sample}/"
    conda: "../wrappers/multiqc_report/env.yaml"
    script: "../wrappers/multiqc_report/script.py"


