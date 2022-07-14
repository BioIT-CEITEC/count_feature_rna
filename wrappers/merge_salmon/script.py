######################################
# wrapper for rule: merge_salmon
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: merge_salmon \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.output.tsv)+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

merge_salmon = os.path.abspath(os.path.dirname(__file__))+ "/merge_tables.R"

command = "(time Rscript "+merge_salmon+" "+snakemake.input+" "+ snakemake.output.tsv + ") >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
