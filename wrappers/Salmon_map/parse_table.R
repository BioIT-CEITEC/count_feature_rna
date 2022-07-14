library(data.table)
library(jsonlite)

run_all <- function(args){

  info_file <- args[1]
  sample_name <- args[2]
  output <- args[3]

  salmon.tbl<-as.data.table(read_json(info_file))
  salmon.tbl<-salmon.tbl[, Samples:=sample_name]
  salmon.tbl<-unique(salmon.tbl[, .(Samples, Total=num_processed, mapped=num_mapped, unmapped=num_processed-num_mapped)])
  fwrite(file=output,x=salmon.tbl,sep="\t")

}

# run as Rscript
#script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)