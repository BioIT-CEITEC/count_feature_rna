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

rule Salmon_bam:
    input:  bam = "mapped/{sample}.bam",
            ref_transcriptome = expand("{ref_dir}/index/RSEM/{ref}.idx.fa",ref_dir=reference_directory,ref=config["reference"])[0],
    output: sf = "qc_reports/{sample}/salmon_bam/{sample}.salmon_bam.sf",
    log:    "logs/{sample}/salmon_bam.log"
    threads: 5
    resources:  mem = 10
    params: prefix = "qc_reports/{sample}/salmon_bam/{sample}",
            lib_type = config["lib_type"],
    conda:  "../wrappers/Salmon_bam/env.yaml"
    script: "../wrappers/Salmon_bam/script.py"

rule Salmon_fastq:
    input:  unpack(salmon_kallisto_input),
            index = expand("{ref_dir}/index/Salmon_fastq",ref_dir=reference_directory,ref=config["reference"])
    output: sf = "qc_reports/{sample}/salmon_fastq/{sample}.salmon_fastq.sf",
    log:    "logs/{sample}/salmon_fastq.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "qc_reports/{sample}/salmon_fastq/{sample}",
            lib_type = config["lib_type"],
            gcbias = config["gcbias"],
            numGibbsSamples = config["numGibbsSamples"],
            paired = paired,
    conda: "../wrappers/Salmon_fastq/env.yaml"
    script: "../wrappers/Salmon_fastq/script.py"

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

