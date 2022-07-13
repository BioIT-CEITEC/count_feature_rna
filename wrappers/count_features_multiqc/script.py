######################################
# wrapper for rule: count_features_multiqc
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: count_features_multiqc \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

multiqc_search_paths = "./qc_reports/*/feature_count/*" + " ./qc_reports/{sample}/RSEM/*" + " ./qc_reports/{sample}/salmon_bam/*" + " ./qc_reports/{sample}/salmon_fastq/*" + " ./qc_reports/{sample}/kallisto/*"

command = "multiqc -f -n " + snakemake.output.html + " " + multiqc_search_paths + \
              " --cl_config \"{{read_count_multiplier: 0.001, read_count_prefix: 'K', read_count_desc: 'thousands' }}\" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

