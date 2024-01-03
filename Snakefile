import pandas as pd
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

read_pair_tags = BR.set_read_pair_qc_tags()
pair_tag = BR.set_read_pair_tags()
paired = BR.set_paired_tags()

config = BR.load_organism()

# DNA parameters processing
#
if not "lib_ROI" in config:
    config["lib_ROI"] = "rna"

# RNA parameters processing
#
if not "strandness" in config:
    config["strandness"] = "unstr"

if not "featureCount" in config:
    config["featureCount"] = False

if not "RSEM" in config:
    config["RSEM"] = False

if not "kallisto" in config:
    config["kallisto"] = False

if not "salmon_map" in config:
    config["salmon_map"] = False

if not "salmon_align" in config:
    config["salmon_map"] = False

## featureCount count_over option
if not "count_over" in config:
    config["count_over"] = "exon"

count_over_list = config['count_over'].split(",")

## Salmon parameters
if not "lib_type" in config:
    config["lib_type"] = "A"

if not "gcbias" in config:
    config["gcbias"] = True

if not "numGibbsSamples" in config:
    config["numGibbsSamples"] = "20"

# ChIP-seq parameters processing
#
if not "effective_genome_size" in config:
    config["effective_genome_size"] = "unk"

if not "fragment_length" in config:
    config["fragment_length"] = "unk"
    
if not "summary_correlation_method" in config:
    config["summary_correlation_method"] = "spearman"

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "R1|R2|SE",
    pair_tags = "|".join(pair_tag),
    count_over_list = "exon|gene|transcript|three_prime_UTR|five_prime_UTR"

##### Target rules #####

rule all:
    input:  "qc_reports/all_samples/multiqc_features.html"

##### Modules #####

include: "rules/count_features.smk"
