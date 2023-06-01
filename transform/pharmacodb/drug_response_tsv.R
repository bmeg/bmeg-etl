suppressMessages(library(optparse))
suppressMessages(library(stringr))

#treatmentResponse@info: inhibitor, sampleid, inhibitor_panel, time_of_read, treatmentid, nbr.conc.tested (whatever that means) 
#treatmentResponse@raw: dose and viability for each aliquot
#treatmentResponse@profiles: final stats

#arg handling
option_list = list(
  make_option(c("-o", "--out"), type="character", default=NULL, help="tsv output destination", metavar="character"),
  make_option(c("-s", "--source"), type="character", default=NULL, help="\"source\" directory where pharmacoGx objects will be stored. Will create a new directory doesn't exist yet", metavar="character"),
  make_option(c("-r", "--resume"), type="character", default=NULL, help="tsv output from previous run to resume", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Missing parameter: tsv output destination", call.=FALSE)
}
outFile = opt$out
cat(c("Output:", outFile, "\n"))

if (is.null(opt$source)){
  print_help(opt_parser)
  stop("Missing parameter: source storage directory", call.=FALSE)
}
sourceDir = opt$source
cat(c("Source directory:", sourceDir, "\n"))

suppressMessages(library(PharmacoGx))

friendlyRbind <- function(x, y, id){
  #allows the merging of two tables with some shared column names and no identical row names
  if(ncol(x)==0){
    return(y)
  }
  if(ncol(y)==0){
    return(x)
  }
  sharedCols = intersect(colnames(x),colnames(y))
  if( length(sharedCols)<1){
    stop("Improper use of friendlyRbind, no overlap between x and y column names")
  }
  xexclusive = x[,c(id, setdiff(colnames(x), sharedCols))]
  yexclusive = y[,c(id, setdiff(colnames(y), sharedCols))]
  result = rbind(x[,sharedCols],y[,sharedCols])
  if(ncol(x)>length(sharedCols)){
    result = merge(x = xexclusive, y = result, by=id, all=TRUE)
  }
  if(ncol(y)>length(sharedCols)){
    result = merge(x = result, y = yexclusive, by=id, all=TRUE)
  }
  return(result)
}


cat("library loaded\n")

names = availablePSets()$`PSet Name`

CASE_IDS = c('Cellosaurus.Accession.id', 'PatientId', 'depmap_id', 'Cell line', 'Cosmic ID', 'sampleid','MODEL')

if (is.null(opt$resume)){
  cat("No resume value given, starting from the top")
  big = data.frame()
  remaining = names
}else{
  cat("Resuming from file", opt$resume)
  big = read.table(opt$resume, header=TRUE, sep="\t", comment.char="")
  done = unique(big$project)
  remaining = setdiff(names, done)
}
cat(",", as.character(length(remaining)), "sets to process\n")

for(setName in remaining){
  cat("Loading", setName, "\n")

  PSET = downloadPSet(setName, saveDir=sourceDir, pSetFileName=setName, timeout=1200)

  cat("Processing PSet:\n")
  info = PSET@treatmentResponse$info
  raw = PSET@treatmentResponse$raw
  profiles = PSET@treatmentResponse$profiles
  colnames(profiles)=paste("profileInfo_", toupper(colnames(profiles)), sep="")
  sampleFields = names(PSET@sample)
  applicableFields = intersect(CASE_IDS, sampleFields)
  if(length(applicableFields >= 2)){
    cat("Applicable sample/case ID fields include ", applicableFields, "\n")
    samples = PSET@sample[applicableFields]
  }else{
    cat("Not enough fields found. Using the first two columns")
    samples = PSET@sample[1:2]
  }

  treatments = PSET@curation$treatment

  cat(" -converting raw data array to data.frames\n")  
  doses = as.data.frame(raw[,,"Dose"])
  colnames(doses) = paste("rawDdose",1:ncol(doses),sep="")
  doses$responseID = paste(setName,rownames(doses),sep="/")
  viabilities = as.data.frame(raw[,,"Viability"])
  colnames(viabilities) = paste("rawVdose",1:ncol(viabilities),sep="")
  viabilities$responseID = paste(setName,rownames(viabilities),sep="/")
  dv = merge(doses, viabilities, by="responseID") 
  secondii = merge(info,profiles, by.x='row.names', by.y='row.names')

  cat(" -renaming columns\n")  
  samples$sampleid = rownames(samples)
  cat("Treatment names before:", names(treatments), "\n")
  treatments = treatments[,1:2]
  names(treatments)[which(names(treatments)!="unique.treatmentid")] = "projecttreatmentid"
  cat("Treatment names after:", names(treatments), "\n")
  
  treatments = treatments[which(treatments$projecttreatmentid != ""),]
  twoPart = treatments[str_detect(treatments$projecttreatmentid,"///"),]
  treatments = treatments[!(str_detect(treatments$projecttreatmentid,"///")),]
  cat("Length of twoPart:", nrow(twoPart), "\n")
  if(nrow(twoPart) >0){
    split = str_split(twoPart$projecttreatmentid,"///", simplify=TRUE)
    firstHalf = twoPart
    firstHalf$projecttreatmentid = split[,1]
    secondHalf = twoPart
    secondHalf$projecttreatmentid = split[,2]
    append = rbind(firstHalf,secondHalf)
    treatments = rbind(treatments, append)
  }

  secondi = merge(secondii,samples, by="sampleid")
  second = merge(secondi,treatments, by.x="treatmentid", by.y="projecttreatmentid", all.x = TRUE)

  second$project = setName
  colnames(second)[which(names(second)=="Row.names")] = "experimentID"
  second$responseID = paste(second$project, second$experimentID, sep = "/")
  final = merge(second, dv, by="responseID")

  cat("Appending to big data.frame\n")  
  big = friendlyRbind(big,final,"responseID")

  cat("Saving to files\n")
}

cat("done\n")

