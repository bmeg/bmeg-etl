suppressMessages(library(optparse))

option_list = list(
  make_option(c("-s", "--source"), type="character", default=NULL, help="\"source\" directory where pharmacoGx objects will be stored. Will create a new directory doesn't exist yet", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$source)){
  print_help(opt_parser)
  stop("Missing parameter: source storage directory", call.=FALSE)
}
sourceDir = opt$source
cat(c("Source directory:", sourceDir, "\n"))

suppressMessages(library(PharmacoGx))
cat("library loaded\n")

logLoc = paste(sourceDir, "downloadLog.tsv", sep="/")
if(file.exists(logLoc)){
  log = read.table(logLoc, header=TRUE, sep="\t", comment.char="")
  failed = log[log$status == "started","name"]
  for(fname in failed){
    path = paste(sourceDir, fname, sep="/")
    cat("Removing failed download:", path, "\n")
    log = log[log$name != fname,]
    file.remove(path)
  }
  cat("Completed jobs:\n")
  completed = log$name
  cat(completed, "\n")
}else{
  cat("No log file detected, starting from scratch\n")
  log = NULL
  completed = NULL  
}  

names = availablePSets()$`PSet Name`
remaining = setdiff(names, completed)

for(setName in remaining){
  cat("Loading", setName, "\n")
  write.table(rbind(log, list(name=setName, status = "started")), file=logLoc, sep = "\t", row.names = FALSE, na="")
  PSET = downloadPSet(setName, saveDir=sourceDir, pSetFileName=setName, timeout=2400)
  cat(setName, "loaded successfully!\n")
  log = rbind(log, list(name=setName, status = "completed"))
}
cat("Done\n")
