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

#output location is required
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Missing parameter: tsv output destination", call.=FALSE)
}
outFile = opt$out
cat(c("Output:", outFile, "\n"))

#storage location for pharmacogx objects is required
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

#finding where we left off if this program was run before
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

#Processing each PSet remaining
for(setName in remaining){
  cat("Loading", setName, "\n")

  PSET = downloadPSet(setName, saveDir=sourceDir, pSetFileName=setName, timeout=1200)

  cat("Processing PSet:\n")
  #Pulling tables from the pharmacoGx object. info = sample and drug info for each experiment, raw = raw dose and viability data, profiles = drug response profiles, samples = sample info, treatments = standardized treatment names
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

  cat("Rearranging data.frames\n")  
  #raw array is split into two frames and changing the column names
  doses = as.data.frame(raw[,,"Dose"])
  colnames(doses) = paste("rawDdose",1:ncol(doses),sep="")
  doses$responseID = paste(setName,rownames(doses),sep="/")
  viabilities = as.data.frame(raw[,,"Viability"])
  colnames(viabilities) = paste("rawVdose",1:ncol(viabilities),sep="")
  viabilities$responseID = paste(setName,rownames(viabilities),sep="/")

  #merging the dose and viabilities tables together
  dv = merge(doses, viabilities, by="responseID") 

  #merging the info and drug response profiles together
  secondii = merge(info,profiles, by.x='row.names', by.y='row.names')

  #re-formatting the treatment data frame so it is easier to merge with the info/profile table
  samples$sampleid = rownames(samples)
  treatments = treatments[,1:2]
  names(treatments)[which(names(treatments)!="unique.treatmentid")] = "projecttreatmentid"
  treatments = treatments[which(treatments$projecttreatmentid != ""),]
  twoPart = treatments[str_detect(treatments$projecttreatmentid,"///"),]
  treatments = treatments[!(str_detect(treatments$projecttreatmentid,"///")),]
  if(nrow(twoPart) >0){
    split = str_split(twoPart$projecttreatmentid,"///", simplify=TRUE)
    firstHalf = twoPart
    firstHalf$projecttreatmentid = split[,1]
    secondHalf = twoPart
    secondHalf$projecttreatmentid = split[,2]
    append = rbind(firstHalf,secondHalf)
    treatments = rbind(treatments, append)
  }

  #merging the info/pofile table with the samples and treatments tables
  secondi = merge(secondii,samples, by="sampleid")
  second = merge(secondi,treatments, by.x="treatmentid", by.y="projecttreatmentid", all.x = TRUE)

  #adding project, experimentID and responseID to the info/profile/samples/treatments table
  second$project = setName
  colnames(second)[which(names(second)=="Row.names")] = "experimentID"
  second$responseID = paste(second$project, second$experimentID, sep = "/")

  #merging the info/profile/samples/treatments table with the raw dose/viability table
  final = merge(second, dv, by="responseID")

  cat("Appending to big data.frame\n")  
  #adding the data.frame for this set to the big one that represents all datasets on pharmacodb
  big = friendlyRbind(big,second,"responseID")

  cat("Saving to file\n")
  #backing up to a file
  write.table(big, file=outFile, sep = "\t", row.names = FALSE, na="")
}

cat("done\n")

