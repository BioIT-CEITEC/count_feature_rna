def multiqc_report_input_files(wildcards):
    input = {}
    if config["featureCount_exon"]:
        input["featureCount_exon"] = expand("qc_reports/{sample}/featureCount_exon/{sample}.featureCount_exon.tsv", sample = sample_tab.sample_name)
    if config["featureCount_gene"]:
        input["featureCount_gene"] = expand("qc_reports/{sample}/featureCount_gene/{sample}.featureCount_gene.tsv",sample=sample_tab.sample_name)
    if config["featureCount_transcript"]:
        input["featureCount_transcript"] = expand("qc_reports/{sample}/featureCount_transcript/{sample}.featureCount_transcript.tsv",sample=sample_tab.sample_name)
    if config["featureCount_3pUTR"]:
        input["featureCount_3pUTR"] = expand("qc_reports/{sample}/featureCount_3pUTR/{sample}.featureCount_3pUTR.tsv",sample=sample_tab.sample_name)
    if config["featureCount_5pUTR"]:
        input["featureCount_5pUTR"] = expand("qc_reports/{sample}/featureCount_5pUTR/{sample}.featureCount_5pUTR.tsv",sample=sample_tab.sample_name)
    if config["RSEM"]:
        input["RSEM"] = expand("qc_reports/{sample}/RSEM/{sample}.genes.results", sample = sample_tab.sample_name)
    if config["salmon_align"]:
        input["salmon_align"] = expand("qc_reports/{sample}/salmon_aln/{sample}.salmon_aln.sf", sample = sample_tab.sample_name)
        input["salmon_align_tab"] = expand("qc_reports/{sample}/salmon_aln/{sample}_aln.tsv", sample = sample_tab.sample_name)
    if config["salmon_map"]:
        input["salmon_map"] = expand("qc_reports/{sample}/salmon_map/{sample}.salmon_map.sf", sample = sample_tab.sample_name)
        input["salmon_map_tab"] = expand("qc_reports/{sample}/salmon_map/{sample}_map.tsv", sample = sample_tab.sample_name)
    if config["kallisto"]:
        input["kallisto_h5"] = expand("qc_reports/{sample}/kallisto/{sample}.kallisto.h5", sample = sample_tab.sample_name)
        input["kallisto_tsv"] = expand("qc_reports/{sample}/kallisto/{sample}.kallisto.tsv", sample = sample_tab.sample_name)
    return input

rule multiqc_report:
    input:  unpack(multiqc_report_input_files),
    output: html="qc_reports/all_samples/multiqc_features.html"
    log:    "logs/all_samples/multiqc_features.log"
    params: multiqc_config = workflow.basedir+"/wrappers/count_features_multiqc/multiqc_config.txt",
            multiqc_path = "qc_reports/all_samples/"
    conda: "../wrappers/count_features_multiqc/env.yaml"
    script: "../wrappers/count_features_multiqc/script.py"

rule featureCount_exon:
     input:  bam = "mapped/{sample}.bam",
             gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
     output: feature_count = "qc_reports/{sample}/featureCount_exon/{sample}.featureCount_exon.tsv"
     log:    "logs/{sample}/featureCount_exon.log"
     threads: 10
     resources:  mem = 10
     params: count_over = "exon",
             paired = paired,
             strandness = config["strandness"],
     conda:  "../wrappers/feature_count/env.yaml"
     script: "../wrappers/feature_count/script.py"

rule featureCount_gene:
     input:  bam = "mapped/{sample}.bam",
             gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
     output: feature_count = "qc_reports/{sample}/featureCount_gene/{sample}.featureCount_gene.tsv"
     log:    "logs/{sample}/featureCount_gene.log"
     threads: 10
     resources:  mem = 10
     params: count_over = "gene",
             paired = paired,
             strandness = config["strandness"],
     conda:  "../wrappers/feature_count/env.yaml"
     script: "../wrappers/feature_count/script.py"

rule featureCount_transcript:
     input:  bam = "mapped/{sample}.bam",
             gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
     output: feature_count = "qc_reports/{sample}/featureCount_transcript/{sample}.featureCount_transcript.tsv"
     log:    "logs/{sample}/featureCount_transcript.log"
     threads: 10
     resources:  mem = 10
     params: count_over = "transcript",
             paired = paired,
             strandness = config["strandness"],
     conda:  "../wrappers/feature_count/env.yaml"
     script: "../wrappers/feature_count/script.py"

rule featureCount_3pUTR:
     input:  bam = "mapped/{sample}.bam",
             gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
     output: feature_count = "qc_reports/{sample}/featureCount_3pUTR/{sample}.featureCount_3pUTR.tsv"
     log:    "logs/{sample}/featureCount_3pUTR.log"
     threads: 10
     resources:  mem = 10
     params: count_over = "three_prime_utr",
             paired = paired,
             strandness = config["strandness"],
     conda:  "../wrappers/feature_count/env.yaml"
     script: "../wrappers/feature_count/script.py"

rule featureCount_5pUTR:
     input:  bam = "mapped/{sample}.bam",
             gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
     output: feature_count = "qc_reports/{sample}/featureCount_5pUTR/{sample}.featureCount_5pUTR.tsv"
     log:    "logs/{sample}/featureCount_5pUTR.log"
     threads: 10
     resources:  mem = 10
     params: count_over = "five_prime_utr",
             paired = paired,
             strandness = config["strandness"],
     conda:  "../wrappers/feature_count/env.yaml"
     script: "../wrappers/feature_count/script.py"

rule RSEM:
    input:  bam = "mapped/{sample}.bam",
            transcriptome = "mapped/transcriptome/{sample}.transcriptome.bam",
            rsem_index = expand("{ref_dir}/index/RSEM/{ref}.idx.fa",ref_dir=reference_directory,ref=config["reference"])[0],
    output: rsem_out = "qc_reports/{sample}/RSEM/{sample}.genes.results"
    log:    "logs/{sample}/RSEM.log"
    threads: 5
    resources:  mem = 10
    params: paired = paired,
            strandness = config["strandness"],
    conda:  "../wrappers/RSEM/env.yaml"
    script: "../wrappers/RSEM/script.py"


def salmon_kallisto_input(wildcards):
    preprocessed = "cleaned_fastq"
    input = {}
    if not config["is_paired"]:
        input['r1'] = os.path.join(preprocessed,"{sample}.fastq.gz")
    else:
        input['r1'] = os.path.join(preprocessed,"{sample}_R1.fastq.gz")
        input['r2'] = os.path.join(preprocessed,"{sample}_R2.fastq.gz")
    return input

rule Salmon_align:
    input:  bam = "mapped/transcriptome/{sample}.transcriptome.bam",
            cds = expand("{ref_dir}/seq/{ref}.cds.fa",ref_dir=reference_directory,ref=config["reference"])[0],
    output: sf = "qc_reports/{sample}/salmon_aln/{sample}.salmon_aln.sf",
            tsv= "qc_reports/{sample}/salmon_aln/{sample}_aln.tsv",
    log:    "logs/{sample}/salmon_aln.log"
    threads: 5
    resources:  mem = 10
    params: prefix = "qc_reports/{sample}/salmon_aln",
            lib_type = config["lib_type"],
            sample_name= "{sample}_aln",
            info="qc_reports/{sample}/salmon_aln/aux_info/meta_info.json",
    conda:  "../wrappers/Salmon_align/env.yaml"
    script: "../wrappers/Salmon_align/script.py"

rule Salmon_map:
    input:  unpack(salmon_kallisto_input),
            index = expand("{ref_dir}/index/Salmon",ref_dir=reference_directory,ref=config["reference"])[0],
    output: sf = "qc_reports/{sample}/salmon_map/{sample}.salmon_map.sf",
            tsv = "qc_reports/{sample}/salmon_map/{sample}_map.tsv",
    log:    "logs/{sample}/salmon_map.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "qc_reports/{sample}/salmon_map",
            lib_type = config["lib_type"],
            sample_name = "{sample}_map",
            info = "qc_reports/{sample}/salmon_map/aux_info/meta_info.json",
            gcbias = config["gcbias"],
            numGibbsSamples = config["numGibbsSamples"],
            paired = paired,
    conda: "../wrappers/Salmon_map/env.yaml"
    script: "../wrappers/Salmon_map/script.py"

rule Kallisto:
    input:  unpack(salmon_kallisto_input),
            index = expand("{ref_dir}/index/Kallisto",ref_dir=reference_directory,ref=config["reference"])
    output: h5 = "qc_reports/{sample}/kallisto/{sample}.kallisto.h5",
            tsv = "qc_reports/{sample}/kallisto/{sample}.kallisto.tsv"
    log:    "logs/{sample}/kallisto.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "qc_reports/{sample}/kallisto",
            samplelog = "qc_reports/{sample}/kallisto/{sample}.kallisto.log",
            paired = paired,
    conda: "../wrappers/Kallisto/env.yaml"
    script: "../wrappers/Kallisto/script.py"
