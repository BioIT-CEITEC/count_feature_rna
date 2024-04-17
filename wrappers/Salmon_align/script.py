######################################
# wrapper for rule: Salmon_align
######################################
import os
import subprocess
import glob
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: Salmon_align \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "salmon quant -p " + str(snakemake.threads) + \
               " -t " + snakemake.input.cds + \
               " -l " + snakemake.params.lib_type + \
               " -a " + snakemake.input.bam + \
               " -o " + snakemake.params.prefix

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "cat " + snakemake.params.prefix + "/logs/salmon_quant.log >> " +log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + snakemake.params.prefix + "/quant.sf " + snakemake.output.sf + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

parse_salmon = os.path.abspath(os.path.dirname(__file__))+ "/parse_table.R"

command = "(time Rscript " + parse_salmon + " " + \
          snakemake.params.info + " " + \
          snakemake.params.sample_name + " " + \
          snakemake.output.tsv + ") >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: parsing results\n")
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
