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

sorted_bam = snakemake.input.transcriptome_bam.replace("not_markDups.transcriptome.bam","sorted.bam")
ddup_bam = snakemake.input.transcriptome_bam.replace("not_markDups.transcriptome.bam","ddup.bam")
header_sam = snakemake.input.transcriptome_bam.replace("not_markDups.transcriptome.bam","header.sam")
ddup_txt = snakemake.input.transcriptome_bam.replace("not_markDups.transcriptome.bam","ddup.txt")
ddup_sam = snakemake.input.transcriptome_bam.replace("not_markDups.transcriptome.bam","ddup.sam")

if snakemake.params.mark_duplicates == True:
    os.makedirs(os.path.dirname(snakemake.params.mtx), exist_ok=True)

    command = "samtools sort -@" + str(snakemake.threads) + " -o " + sorted_bam + " " + snakemake.input.transcriptome_bam
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "samtools index -@" + str(snakemake.threads) + " " + sorted_bam
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    if snakemake.params.UMI == "no_umi" or snakemake.params.umi_usage == "no":

        command = "export LD_BIND_NOW=1"
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

        command = "picard MarkDuplicates INPUT=" + sorted_bam + " OUTPUT=" + ddup_bam + \
                  " METRICS_FILE=" + snakemake.params.mtx + " REMOVE_DUPLICATES=" + str(snakemake.params.rmDup) + \
                  " ASSUME_SORTED=true PROGRAM_RECORD_ID=null VALIDATION_STRINGENCY=LENIENT" + \
                  " -Xmx" + str(snakemake.resources.mem) + "g 2>> " + log_filename + " "
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

    else:

        command = "umi_tools dedup -I " + sorted_bam + " -S " + ddup_bam + \
                  " --log " + snakemake.params.mtx + extra + \
                  " --extract-umi-method=read_id --umi-separator='_' --method=directional --edit-distance-threshold=0" + \
                  " --spliced-is-unique --multimapping-detection-method=NH"
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

    command = "rm " + sorted_bam + "*"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "samtools view -H " + snakemake.input.transcriptome_bam + " > " + header_sam
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "samtools view -@ " + str(snakemake.threads) + " " + ddup_bam + \
              " | cut -f 1 | sort | uniq > " + ddup_txt
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "samtools view -@ " + str(snakemake.threads) + " " + snakemake.input.transcriptome_bam + \
              " | grep -F -f " + ddup_txt + \
              " > " + ddup_sam
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "cat " + header_sam + " " + ddup_sam + \
              " | samtools view -@ " + str(snakemake.threads) + " -b -o " + snakemake.output.bam
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "rm " + header_sam + " " + ddup_txt + " " + ddup_sam
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    if snakemake.params.keep_not_markDups_bam == False:
        command = "rm " + snakemake.input.transcriptome_bam
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

else:

    command = "mv -T " + snakemake.input.transcriptome_bam + " " + snakemake.output.bam
    f = open(log_filename, 'at')
    f.write("## No markduplicate was requested" + "\n")
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
