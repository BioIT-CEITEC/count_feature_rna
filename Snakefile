import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_TMPD_PATH = "./tmp/"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

# DNA parameteres processing
#
if not "lib_ROI" in config:
    config["lib_ROI"] = "wgs"

# RNA parameteres processing
#
if not "strandness" in config:
    config["strandness"] = "unstr"

if not "count_over" in config:
    config["count_over"] = "exon"

if not "feature_count" in config:
    config["feature_count"] = False

if not "RSEM" in config:
    config["RSEM"] = False

if not "kallisto" in config:
    config["kallisto"] = False

if not "salmon" in config:
    config["salmon"] = False

if not "salmon_map" in config:
    config["salmon_map"] = False

if not "salmon_align" in config:
    config["salmon_map"] = False
    
# ChIP-seq parameters processing
#
if not "effective_genome_size" in config:
    config["effective_genome_size"] = "unk"

if not "fragment_length" in config:
    config["fragment_length"] = "unk"
    
if not "summary_correlation_method" in config:
    config["summary_correlation_method"] = "spearman"

# Reference processing
#
if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()
config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
if len(config["species_name"].split(" (")) > 1:
    config["species"] = config["species_name"].split(" (")[1].replace(")","")


##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if not config["is_paired"]:
    read_pair_tags = ["SE"]
    pair_tag = [""]
    paired = "SE"
else:
    read_pair_tags = ["_R1","_R2"]
    pair_tag = ["_R1","_R2"]
    paired = "PE"

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "R1|R2|SE",
    pair_tags = "|".join(pair_tag),

##### Target rules #####

rule all:
    input:  "qc_reports/all_samples/multiqc_features.html"

##### Modules #####

include: "rules/count_features.smk"
