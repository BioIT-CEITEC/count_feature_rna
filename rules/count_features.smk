def multiqc_report_input_files(wildcards):
    input = {}
    if config["feature_count"]:
        input["feature_count"] = expand("qc_reports/{sample}/feature_count/{sample}.feature_count.tsv", sample = sample_tab.sample_name)
    if config["RSEM"]:
        input["RSEM"] = expand("qc_reports/{sample}/RSEM/{sample}.genes.results", sample = sample_tab.sample_name)
    if config["salmon_align"]:
        input["salmon_align"] = expand("qc_reports/{sample}/salmon/{sample}_aln/{sample}.salmon_aln.sf", sample = sample_tab.sample_name)
        input["salmon_align_tab"] = expand("qc_reports/{sample}/salmon/{sample}_aln/{sample}_aln.tsv", sample = sample_tab.sample_name)
    if config["salmon_map"]:
        input["salmon_map"] = expand("qc_reports/{sample}/salmon/{sample}_map/{sample}.salmon_map.sf", sample = sample_tab.sample_name)
        input["salmon_map_tab"] = expand("qc_reports/{sample}/salmon/{sample}_map/{sample}_map.tsv", sample = sample_tab.sample_name)
    if config["kallisto"]:
        input["kallisto_h5"] = expand("qc_reports/{sample}/kallisto/{sample}.kallisto.h5", sample = sample_tab.sample_name)
        input["kallisto_tsv"] = expand("qc_reports/{sample}/kallisto/{sample}.kallisto.tsv", sample = sample_tab.sample_name)
    return input

rule multiqc_report:
    input:  unpack(multiqc_report_input_files),
    output: html="qc_reports/all_samples/multiqc_features.html"
    log:    "logs/all_samples/multiqc_features.log"
    params: multiqc_config = workflow.basedir+"/wrappers/multiqc_report/multiqc_config.txt",
            multiqc_path = "qc_reports/all_samples/"
    conda: "../wrappers/count_features_multiqc/env.yaml"
    script: "../wrappers/count_features_multiqc/script.py"

rule feature_count:
     input:  bam = "mapped/{sample}.bam",
             gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
     output: feature_count = "qc_reports/{sample}/feature_count/{sample}.feature_count.tsv"
     log:    "logs/{sample}/feature_count.log"
     threads: 10
     resources:  mem = 10
     params: count_over = config["count_over"],
             paired = paired,
             strandness = config["strandness"],
     conda:  "../wrappers/feature_count/env.yaml"
     script: "../wrappers/feature_count/script.py"


rule RSEM:
    input:  bam = "mapped/{sample}.bam",
            transcriptome = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam",
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
    input:  bam = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam",
            cds = expand("{ref_dir}/seq/{ref}.cds.fa",ref_dir=reference_directory,ref=config["reference"])[0],
    output: sf = "qc_reports/{sample}/salmon/{sample}_aln/{sample}.salmon_aln.sf",
            tsv= "qc_reports/{sample}/salmon/{sample}_aln/{sample}_aln.tsv",
    log:    "logs/{sample}/salmon_align.log"
    threads: 5
    resources:  mem = 10
    params: prefix = "qc_reports/{sample}/salmon/{sample}_aln",
            lib_type = config["lib_type"],
            sample_name= "{sample}_aln",
            info="qc_reports/{sample}/salmon/{sample}_aln/aux_info/meta_info.json",
    conda:  "../wrappers/Salmon_align/env.yaml"
    script: "../wrappers/Salmon_align/script.py"

rule Salmon_map:
    input:  unpack(salmon_kallisto_input),
            index = expand("{ref_dir}/index/Salmon",ref_dir=reference_directory,ref=config["reference"])[0],
    output: sf = "qc_reports/{sample}/salmon/{sample}_map/{sample}.salmon_map.sf",
            tsv = "qc_reports/{sample}/salmon/{sample}_map/{sample}_map.tsv",
    log:    "logs/{sample}/salmon_map.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "qc_reports/{sample}/salmon/{sample}_map",
            lib_type = config["lib_type"],
            sample_name = "{sample}_map",
            info = "qc_reports/{sample}/salmon/{sample}_map/aux_info/meta_info.json",
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

