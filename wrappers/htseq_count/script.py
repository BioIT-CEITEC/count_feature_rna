######################################
# wrapper for rule: htseq_count
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)
summary_file = snakemake.output.htseq_count+".summary"

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: gene_counts_HTSeqCounts \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.output.htseq_count)+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if snakemake.params.paired == "PE" :
    extra_flags_feature = " -r pos "
else:
    extra_flags_feature = ""

if snakemake.params.strandness == "fwd":
    extra_flags_feature += " -s yes"
elif snakemake.params.strandness == "rev":
    extra_flags_feature += " -s reverse"
else:
    extra_flags_feature += " -s no"

command = "htseq-count -t "+snakemake.params.count_over+" -i gene_id "+extra_flags_feature+" -m "+snakemake.params.mode+" --nonunique="+snakemake.params.nonunique+" -c "+snakemake.output.htseq_count+" "+snakemake.input.bam+" "+snakemake.input.gtf+" >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "sed -i 's|mapped/|mapped/" + snakemake.params.count_over + "_|' " + summary_file
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
