suppressMessages(library(optparse))
suppressMessages(library(stringr))

#treatmentResponse@info: inhibitor, sampleid, inhibitor_panel, time_of_read, treatmentid, nbr.conc.tested (whatever that means) 
#treatmentResponse@raw: dose and viability for each aliquot
#treatmentResponse@profiles: final stats

#arg handling
option_list = list(
  make_option(c("-o", "--out"), type="character", default=NULL, help="tsv output directory", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, help="file for pharmacoGx to run on", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);

#output location is required
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Missing parameter: tsv output destination", call.=FALSE)
}
outDir = opt$out
cat(c("Output directory:", outDir, "\n"))

#storage location for pharmacogx objects is required
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Missing parameter: input file", call.=FALSE)
}
sourcePath = opt$input

#Processing each PSet remaining
cat("Loading", sourcePath, "\n")

PSET = readRDS(sourcePath)

cat("Processing PSet:\n")
#Pulling tables from the pharmacoGx object. info = sample and drug info for each experiment, raw = raw dose and viability data, profiles = drug response profiles, samples = sample info, treatments = standardized treatment names
info = PSET@treatmentResponse$info
raw = PSET@treatmentResponse$raw
profiles = PSET@treatmentResponse$profiles
colnames(profiles)=paste("profileInfo_", toupper(colnames(profiles)), sep="")
samples = PSET@sample
treatments = PSET@curation$treatment

setNameShort="blaa"

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
second = merge(info,profiles, by.x='row.names', by.y='row.names')
#re-formatting the treatment data frame so it is easier to merge with the info/profile table
samples$sampleid = rownames(samples)

#cat(names(treatments),"\n")
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
#second = merge(secondii,treatments, by.x="treatmentid", by.y="PROJECTtreatmentid", all.x = TRUE)
#adding project, experimentID and responseID to the info/profile/samples/treatments table
second$project = setNameShort
samples$project = setNameShort
colnames(second)[which(names(second)=="Row.names")] = "experimentID"
second$responseID = paste(second$project, second$experimentID, sep = "/")

#merging the info/profile/samples/treatments table with the raw dose/viability table
final = merge(second, dv, by="responseID")

cat("Appending to big data.frame\n")  
#adding the data.frame for this set to the big one that represents all datasets on pharmacodb
#if(length(intersect(names(big), names(final))) == 0){
#  big = final
#  bigSample = samples
#  bigTreatment = treatments
#}else{  
#  big = friendlyRBind(big,final)
#  bigSample = friendlyRBind(bigSample, samples)
#  bigTreatment = rbind(bigTreatment, treatments)
#}

cat("Saving to file\n")
if (!(endsWith(outDir,'/'))){
  outDir = paste0(outDir, '/')
}
#backing up to a file
write.table(final, file=paste0(outDir, "response.tsv"), sep = "\t", row.names = FALSE, na="", qmethod="d")
write.table(samples, file=paste0(outDir, "samples.tsv"), sep = "\t", row.names = FALSE, na="", qmethod="d")
write.table(treatments, file=paste0(outDir, "treatments.tsv"), sep = "\t", row.names = FALSE, na="", qmethod="d")
cat("done\n")
