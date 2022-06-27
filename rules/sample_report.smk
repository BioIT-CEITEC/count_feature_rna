# REPORTING RULES
#


def multiqc_report_input(wildcards):
    input = {}
    if wildcards.sample != "all_samples":
        if paired == "PE":
            input['raw_fastq_R1_report'] = "qc_reports/" + wildcards.sample + "/raw_fastqc/R1_fastqc.zip"
            input['raw_fastq_R2_report'] = "qc_reports/" + wildcards.sample + "/raw_fastqc/R2_fastqc.zip"
        else:
            input['raw_fastq_SE_report'] = "qc_reports/" + wildcards.sample + "/raw_fastqc/SE_fastqc.zip"
        if config["qc_qualimap_DNA"]:
            input['qc_qualimap_DNA'] = "qc_reports/{sample}/qc_qualimap_DNA/{sample}/qualimapReport.html"
        if config["qc_samtools"]:
            input['qc_samtools'] = "qc_reports/{sample}/qc_samtools/{sample}.idxstats.tsv"
        if config["qc_picard_DNA"]:
            input['qc_picard_DNA'] = "qc_reports/{sample}/qc_picard_DNA/picard.tsv"
        if config["qc_qualimap_RNA"]:
            input['qc_qualimap_RNA'] = "qc_reports/{sample}/qc_qualimap_RNA/{sample}/qualimapReport.html"
        if config["qc_picard_RNA"]:
            input['qc_picard_RNA'] = "qc_reports/{sample}/qc_picard_RNA/{sample}.npc.pdf"
        if config["feature_count"]:
            input['feature_count'] = "qc_reports/{sample}/feature_count/{sample}.feature_count.tsv"
        if config["qc_fastq_screen_RNA"]:
            input['qc_fastq_screen_RNA'] = expand("qc_reports/{sample}/qc_fastq_screen_RNA/{sample}{read_pair_tag}_screen.png",sample=sample_tab.sample_name,read_pair_tag=read_pair_tags)
        # if config["biobloom"]:
        #     input['biobloom'] = "qc_reports/{sample}/biobloom/{sample}.biobloom_summary.tsv"
        if config["qc_RSeQC_RNA"]:
            input['qc_RSeQC_RNA'] = "qc_reports/{sample}/qc_RSeQC_RNA/{sample}.RSeQC.read_distribution.txt"
        if config["qc_biotypes_RNA"]:
            input['qc_biotypes_RNA'] = "qc_reports/{sample}/qc_biotypes_RNA/{sample}.biotype_counts.txt"
        # if it's DNA or RNA
        if config["qc_qualimap_DNA"]:
            input['trim'] = expand("qc_reports/{sample}/trim_galore/trim_stats{read_pair_tag}.log",sample=sample_tab.sample_name,read_pair_tag=read_pair_tags)
        else:
            input['trim'] = expand("qc_reports/{sample}/trimmomatic/trim_stats.log",sample=sample_tab.sample_name)
    else:
        input['per_sample_reports'] = expand("qc_reports/{sample}/single_sample_alignment_report.html",sample=sample_tab.sample_name)
    return input


rule multiqc_report:
    input:  unpack(multiqc_report_input)
    output: html="qc_reports/{sample}/multiqc.html"
    log:    "logs/{sample}/multiqc.log"
    params: multiqc_config = workflow.basedir+"/wrappers/multiqc_report/multiqc_config.txt",
            multiqc_path = "qc_reports/{sample}/"
    conda: "../wrappers/multiqc_report/env.yaml"
    script: "../wrappers/multiqc_report/script.py"


def merge_single_sample_reports_input(wildcards):
    input = {}
    if config["qc_biotypes_RNA"]:
        input['biotype'] = expand("qc_reports/{sample}/qc_biotypes_RNA/{sample}.biotype_counts.txt",sample=sample_tab.sample_name)
    if config["qc_fastq_screen_RNA"]:
        input['fastq_screen'] = expand("qc_reports/{sample}/qc_fastq_screen_RNA/{sample}{read_pair_tag}_screen.pdf",sample=sample_tab.sample_name,read_pair_tag=read_pair_tags)
    if config["qc_dupradar_RNA"]:
        input['dupraxpbox'] = expand("qc_reports/{sample}/qc_dupradar_RNA/{sample}_duprateExpBoxplot.pdf",sample=sample_tab.sample_name)
        input['exphist'] = expand("qc_reports/{sample}/qc_dupradar_RNA/{sample}_expressionHist.pdf",sample=sample_tab.sample_name)
        input['dupraexpden'] = expand("qc_reports/{sample}/qc_dupradar_RNA/{sample}_duprateExpDens.pdf",sample=sample_tab.sample_name)
        input['multipergene'] = expand("qc_reports/{sample}/qc_dupradar_RNA/{sample}_multimapPerGene.pdf",sample=sample_tab.sample_name)
        input['readdist'] = expand("qc_reports/{sample}/qc_dupradar_RNA/{sample}_readDist.pdf",sample=sample_tab.sample_name)
    if config["chip_extra_qc"]:
        input['phantompeak'] = expand("qc_reports/{sample}/phantompeakqual/{sample}.{dups}.cross-correlation.pdf", sample=sample_tab.sample_name, dups=['no_dups'])
        input['phantompeak_dups'] = expand("qc_reports/{sample}/phantompeakqual/{sample}.{dups}.cross-correlation.pdf", sample=sample_tab.sample_name, dups=['keep_dups'])
    return input

rule merge_single_sample_reports:
    input:  unpack(merge_single_sample_reports_input)
    output: biotype_pdf = "qc_reports/all_samples/qc_biotype_RNA/biotype.pdf",
            fastq_screen_pdf = "qc_reports/all_samples/fastq_screen/fastq_screen.pdf",
            dupraxpbox_pdf = "qc_reports/all_samples/qc_dupradar_RNA/dupraxpbox.pdf",
            exphist_pdf = "qc_reports/all_samples/qc_dupradar_RNA/exphist.pdf",
            dupraexpden_pdf = "qc_reports/all_samples/qc_dupradar_RNA/dupraexpden.pdf",
            multipergene_pdf = "qc_reports/all_samples/qc_dupradar_RNA/multipergene.pdf",
            readdist_pdf = "qc_reports/all_samples/qc_dupradar_RNA/readdist.pdf",
            phantompeak_pdf="qc_reports/all_samples/phantompeakqual/cross-correlation.no_dups.pdf",
            phantompeak_dups_pdf="qc_reports/all_samples/phantompeakqual/cross-correlation.keep_dups.pdf",
    log:    "logs/all_samples/merge_single_sample_reports.log"
    conda: "../wrappers/merge_single_sample_reports/env.yaml"
    script: "../wrappers/merge_single_sample_reports/script.py"


def per_sample_alignment_report_input(wildcards):
    input = {}
    input['multiqc'] = "qc_reports/{sample}/multiqc.html"
    if paired == "PE":
        input['raw_fastq_R1_report'] = "qc_reports/{sample}/raw_fastqc/R1_fastqc.html"
        input['raw_fastq_R2_report'] = "qc_reports/{sample}/raw_fastqc/R2_fastqc.html"
    else:
        input['raw_fastq_SE_report'] = "qc_reports/{sample}/raw_fastqc/SE_fastqc.html"
    if config["qc_qualimap_DNA"]:
        input['qc_qualimap_DNA'] = "qc_reports/{sample}/qc_qualimap_DNA/{sample}/qualimapReport.html"
    if config["qc_picard_RNA"]:
        input['qc_picard_RNA'] = "qc_reports/{sample}/qc_picard_RNA/{sample}.npc.pdf"
    if config["qc_qualimap_RNA"]:
        input['qc_qualimap_RNA'] = "qc_reports/{sample}/qc_qualimap_RNA/{sample}/qualimapReport.html"
    if config["feature_count"]:
        input['feature_count'] = "qc_reports/{sample}/feature_count/{sample}.feature_count.tsv"
    if config["qc_fastq_screen_RNA"]:
        input['qc_fastq_screen_RNA'] = expand("qc_reports/" + wildcards.sample + "/qc_fastq_screen_RNA/" + wildcards.sample + "{read_pair_tag}_screen.png",read_pair_tag=read_pair_tags)
    # if config["biobloom"]:
    #     input['biobloom'] = "qc_reports/{sample}/biobloom/{sample}.biobloom_summary.tsv"
    if config["qc_biotypes_RNA"]:
        input['qc_biotypes_RNA'] = "qc_reports/{sample}/qc_biotypes_RNA/{sample}.biotype_counts.pdf"
    if config["RSEM"]:
        input['RSEM'] = "qc_reports/{sample}/RSEM/{sample}.genes.results"
    if config["salmon"]:
        input['salmon'] = "qc_reports/{sample}/salmon/{sample}.salmon.sf"
    if config["kallisto"]:
        input['kallisto_h5'] = "qc_reports/{sample}/kallisto/{sample}.kallisto.h5",
        input['kallisto_tsv'] = "qc_reports/{sample}/kallisto/{sample}.kallisto.tsv"
    if config["qc_RSeQC_RNA"]:
        input['qc_RSeQC_RNA'] = "qc_reports/{sample}/qc_RSeQC_RNA/{sample}.RSeQC.read_distribution.txt"
    if config["chip_extra_qc"]:
        input['samtools_contam'] = "qc_reports/{sample}/qc_samtools/{sample}.no_contam.flagstat.tsv",
        input['samtools_dups'] = "qc_reports/{sample}/qc_samtools/{sample}.no_dups.flagstat.tsv",
        input['no_dups_bam_cov'] = "mapped/{sample}.no_dups.bam.bigWig",
        input['phantompeak'] = "qc_reports/{sample}/phantompeakqual/{sample}.no_dups.cross-correlation.pdf",
        input['phantompeak_dups'] = "qc_reports/{sample}/phantompeakqual/{sample}.keep_dups.cross-correlation.pdf",
    return input

rule per_sample_alignment_report:
    input:  unpack(per_sample_alignment_report_input)
    output: sample_report = "qc_reports/{sample}/single_sample_alignment_report.html",
    params: sample_name = "{sample}",
            config = "./config.json",
            paired = paired,
    conda: "../wrappers/per_sample_alignment_report/env.yaml"
    script: "../wrappers/per_sample_alignment_report/script.Rmd"


def final_alignment_report_input(wildcards):
    input = {}
    input['all_sample_multiqc'] = "qc_reports/all_samples/multiqc.html"
    input['per_sample_reports'] = expand("qc_reports/{sample}/single_sample_alignment_report.html",sample=sample_tab.sample_name)
    if config["cross_sample_correlation"]:
        input['cross_sample_correlation'] = "qc_reports/all_samples/cross_sample_correlation/cross_sample_correlation.snps.html"
    if config["qc_biotypes_RNA"]:
        input['qc_biotypes_RNA'] = "qc_reports/all_samples/qc_biotype_RNA/biotype.pdf"
    if config["qc_fastq_screen_RNA"]:
        input['qc_fastq_screen_RNA'] = "qc_reports/all_samples/fastq_screen/fastq_screen.pdf"
    if config["qc_dupradar_RNA"]:
        input['dupraxpbox'] = "qc_reports/all_samples/qc_dupradar_RNA/dupraxpbox.pdf"
        input['exphist'] = "qc_reports/all_samples/qc_dupradar_RNA/exphist.pdf"
        input['dupraexpden'] = "qc_reports/all_samples/qc_dupradar_RNA/dupraexpden.pdf"
        input['multipergene'] = "qc_reports/all_samples/qc_dupradar_RNA/multipergene.pdf"
        input['readdist'] = "qc_reports/all_samples/qc_dupradar_RNA/readdist.pdf"
    if config["chip_extra_qc"]:
        input['phantompeak'] = "qc_reports/all_samples/phantompeakqual/cross-correlation.no_dups.pdf",
        input['phantompeak_dups'] = "qc_reports/all_samples/phantompeakqual/cross-correlation.keep_dups.pdf",
        input['corr_heatmap'] = "qc_reports/all_samples/deeptools/correlation_heatmap.no_dups.pdf",
        input['corr_heatmap_dups'] = "qc_reports/all_samples/deeptools/correlation_heatmap.keep_dups.pdf",
        input["fingerprint"] = "qc_reports/all_samples/deeptools/fingerprint.no_dups.pdf",
        input["fingerprint_dups"] = "qc_reports/all_samples/deeptools/fingerprint.keep_dups.pdf",
    return input

rule final_alignment_report:
    input:  unpack(final_alignment_report_input)
    output: html = "qc_reports/final_alignment_report.html"
    params: sample_name = sample_tab.sample_name,
            config = "./config.json"
    conda: "../wrappers/final_alignment_report/env.yaml"
    script: "../wrappers/final_alignment_report/script.Rmd"
