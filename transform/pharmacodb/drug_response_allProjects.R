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

friendlyRBind = function(x, y){
  xDiff = setdiff(names(x),names(y))
  yDiff = setdiff(names(y),names(x))
  for(xCol in xDiff){
    y[[xCol]] = NA
  }
  for(yCol in yDiff){
    x[[yCol]] = NA
  }
  return(rbind(x,y))
} 


cat("library loaded\n")


names = availablePSets()$`PSet Name`

EXCLUDE = c('GDSC_2020(v1-8.2)','GBM_scr2')

#finding where we left off if this program was run before
if (is.null(opt$resume)){
  cat("No resume value given, starting from the top")
  big = data.frame()
  bigSample = data.frame()
  bigTreatment = data.frame()
  remaining = names
}else{
  cat("Resuming from file", opt$resume)
  big = read.table(opt$resume, header=TRUE, sep="\t", comment.char="")
  done = unique(big$project)
  remaining = setdiff(names, done)
}
remaining = setdiff(remaining, EXCLUDE)
cat(",", as.character(length(remaining)), "sets to process\n")

#Processing each PSet remaining
for(setName in remaining){
  cat("Loading", setName, "\n")
  setNameShort = gsub("_.*", "", setName) 

  PSET = downloadPSet(setName, saveDir=sourceDir, pSetFileName=setName, timeout=1200)

  cat("Processing PSet:\n")
  #Pulling tables from the pharmacoGx object. info = sample and drug info for each experiment, raw = raw dose and viability data, profiles = drug response profiles, samples = sample info, treatments = standardized treatment names
  info = PSET@treatmentResponse$info
  raw = PSET@treatmentResponse$raw
  profiles = PSET@treatmentResponse$profiles
  colnames(profiles)=paste("profileInfo_", toupper(colnames(profiles)), sep="")
  samples = PSET@sample
  treatments = PSET@curation$treatment

  cat("Rearranging data.frames\n")  
  #raw array is split into two frames and changing the column names
  doses = as.data.frame(raw[,,"Dose"])
  colnames(doses) = paste("rawDdose",1:ncol(doses),sep="")
  doses$responseID = paste(setNameShort,rownames(doses),sep="/")
  viabilities = as.data.frame(raw[,,"Viability"])
  colnames(viabilities) = paste("rawVdose",1:ncol(viabilities),sep="")
  viabilities$responseID = paste(setNameShort,rownames(viabilities),sep="/")

  #merging the dose and viabilities tables together
  dv = merge(doses, viabilities, by="responseID") 
  #merging the info and drug response profiles together
  secondii = merge(info,profiles, by.x='row.names', by.y='row.names')
  #re-formatting the treatment data frame so it is easier to merge with the info/profile table
  samples$sampleid = rownames(samples)
  
  cat(names(treatments),"\n")
  treatments = treatments[,1:2]
  names(treatments) = c("UNIQUEtreatmentid", "PROJECTtreatmentid")
  treatments = treatments[which(treatments$PROJECTtreatmentid != ""),]
  twoPart = treatments[str_detect(treatments$PROJECTtreatmentid,"///"),]
  treatments = treatments[!(str_detect(treatments$PROJECTtreatmentid,"///")),]
  if(nrow(twoPart) >0){
    split = str_split(twoPart$PROJECTtreatmentid,"///", simplify=TRUE)
    firstHalf = twoPart
    firstHalf$PROJECTtreatmentid = split[,1]
    secondHalf = twoPart
    secondHalf$PROJECTtreatmentid = split[,2]
    append = rbind(firstHalf,secondHalf)
    treatments = rbind(treatments, append)
  }
  
  #merging the info/pofile table with the samples and treatments tables
  #secondi = merge(secondii,samples, by="sampleid")
  second = merge(secondii,treatments, by.x="treatmentid", by.y="PROJECTtreatmentid", all.x = TRUE)
  #adding project, experimentID and responseID to the info/profile/samples/treatments table
  second$project = setNameShort
  samples$project = setNameShort
  colnames(second)[which(names(second)=="Row.names")] = "experimentID"
  second$responseID = paste(second$project, second$experimentID, sep = "/")

  #merging the info/profile/samples/treatments table with the raw dose/viability table
  final = merge(second, dv, by="responseID")

  cat("Appending to big data.frame\n")  
  #adding the data.frame for this set to the big one that represents all datasets on pharmacodb
  if(length(intersect(names(big), names(final))) == 0){
    big = final
    bigSample = samples
    bigTreatment = treatments
  }else{  
    big = friendlyRBind(big,final)
    bigSample = friendlyRBind(bigSample, samples)
    bigTreatment = rbind(bigTreatment, treatments)
  }
}
cat("Saving to file\n")
#backing up to a file
write.table(big, file=outFile, sep = "\t", row.names = FALSE, na="", qmethod="d")
write.table(bigSample, file=paste(outFile, "samples.tsv", sep="."), sep = "\t", row.names = FALSE, na="", qmethod="d")
write.table(bigTreatment, file=paste(outFile, "treatments.tsv", sep="."), sep = "\t", row.names = FALSE, na="", qmethod="d")
cat("done\n")
