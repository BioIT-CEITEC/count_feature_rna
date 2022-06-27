######################################
# wrapper for rule: Salmon_fastq
######################################
import os
import subprocess
import glob
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: Salmon_fastq \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if snakemake.params.paired == "SE":
    input_fastqs = " -r " + snakemake.input.r1
else:
    input_fastqs = " -1 " + snakemake.input.r1 + " -2 " + snakemake.input.r2

if snakemake.params.gcbias == True:
    salmon_gcbias = " --gcBias "
else:
    salmon_gcbias = ""

command = "salmon quant -p " + str(snakemake.threads) + \
               " -i " + str(snakemake.input.index) + \
               " -l " + snakemake.params.lib_type + input_fastqs + \
               " -o " + snakemake.params.prefix + \
               " --numGibbsSamples " + snakemake.params.numGibbsSamples + salmon_gcbias + \
               " --validateMappings >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)

command = "mv " + snakemake.params.prefix + "/quant.sf " + snakemake.output.sf + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
