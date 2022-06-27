######################################
# wrapper for rule: Kallisto
######################################
import os
import subprocess
import glob
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: Kallisto \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if snakemake.params.paired == "SE":
    input_fastqs = "--single -l 200 -s 25 " + snakemake.input.r1
else:
    input_fastqs = snakemake.input.r1 + " " + snakemake.input.r2

command = "kallisto quant -t " + str(snakemake.threads) + \
               " -i " + str(snakemake.input.index) + \
               " -o " + snakemake.params.prefix + \
               " --seed 12345 " + input_fastqs + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)

command = "mv " + snakemake.params.prefix + "/abundance.h5 " + snakemake.output.h5 + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + snakemake.params.prefix + "/abundance.tsv " + snakemake.output.tsv + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "cp " + log_filename + " " + snakemake.params.samplelog
shell(command)
