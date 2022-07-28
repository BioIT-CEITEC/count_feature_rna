######################################
# wrapper for rule: count_features_multiqc
######################################
import os.path
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

multiqc_search_paths = " ".join([os.path.dirname(i) for i in snakemake.input])

command = "multiqc -f -n " + snakemake.output.html + " " + multiqc_search_paths + \
              " --config " + snakemake.params.multiqc_config + " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

