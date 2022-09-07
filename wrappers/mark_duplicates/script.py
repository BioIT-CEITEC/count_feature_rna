######################################
# wrapper for rule: mark_duplicates
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: mark_duplicates \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.params.mtx)+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if snakemake.params.paired == True:
    extra = " --paired"
else:
    extra = ""

# if snakemake.params.mark_duplicates == True:
#     os.makedirs(os.path.dirname(snakemake.params.mtx), exist_ok=True)
#
#     command = "samtools sort -@" + str(snakemake.threads) + " -o " + snakemake.params.temp_bam + " " + snakemake.input.transcriptome_bam
#     f = open(log_filename, 'at')
#     f.write("## COMMAND: " + command + "\n")
#     f.close()
#     shell(command)
#
#     command = "samtools index -@" + str(snakemake.threads) + " " + snakemake.params.temp_bam
#     f = open(log_filename, 'at')
#     f.write("## COMMAND: " + command + "\n")
#     f.close()
#     shell(command)
#
#     if snakemake.params.UMI == "no_umi" or snakemake.params.umi_usage == "no":
#
#         command = "mv " + snakemake.input.transcriptome_bam + " " + snakemake.output.bam
#         f = open(log_filename, 'at')
#         f.write("## COMMAND: "+command+"\n")
#         f.close()
#         shell(command)
#
#     else:
#
#         command = "umi_tools dedup -I " + snakemake.params.temp_bam + " -S " + snakemake.params.ddup_bam + " --log " + snakemake.params.mtx \
#         + extra + " --extract-umi-method=read_id --umi-separator='_' --method=directional --edit-distance-threshold=0 --spliced-is-unique --multimapping-detection-method=NH"
#         f = open(log_filename, 'at')
#         f.write("## COMMAND: "+command+"\n")
#         f.close()
#         shell(command)
#
#     command = "rm " + snakemake.params.temp_bam + "*"
#     f = open(log_filename, 'at')
#     f.write("## COMMAND: " + command + "\n")
#     f.close()
#     shell(command)
#
#     command = "samtools view -H " + snakemake.input.transcriptome_bam + " > " + snakemake.input.transcriptome_bam.replace("not_markDups.transcriptome.bam", "header.txt")
#     f = open(log_filename, 'at')
#     f.write("## COMMAND: " + command + "\n")
#     f.close()
#     shell(command)
#
#
#
#     if snakemake.params.keep_not_markDups_bam == False:
#         command = "rm " + snakemake.input.bam
#         f = open(log_filename, 'at')
#         f.write("## COMMAND: " + command + "\n")
#         f.close()
#         shell(command)
#
#         command = "rm " + snakemake.input.bam + ".bai"
#         f = open(log_filename, 'at')
#         f.write("## COMMAND: " + command + "\n")
#         f.close()
#         shell(command)
#
# else:

    shell("mv -T " + snakemake.input.transcriptome_bam + " " + snakemake.output.bam)
