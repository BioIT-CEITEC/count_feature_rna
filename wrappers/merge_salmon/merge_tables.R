library(data.table)

run_all <- function(args){

  list_file <- args[1]
  output <- args[2]

  salmon.tbl<-fread(list_file[1])
  for(i in 2:length(list_file)){
    salmon.tbl.nxt<-fread(list_file[i])
    salmon.tbl<-rbind(salmon.tbl,salmon.tbl.nxt)
  }

  fwrite(file=output,x=salmon.tbl,sep="\t")

}

# run as Rscript
#script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)