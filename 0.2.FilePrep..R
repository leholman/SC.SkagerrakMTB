
#############################################
####===Data prep before analysis Skag.MTB===####
####==== Luke E. Holman====15.09.2023====####
#############################################

####====0.0 Packages & Parameters====####
library("Biostrings")

#Set some variables 
minreads <- 2
items <- NULL

#Set the seed 
set.seed("123456")

#Read in metadata
metadata<-read.csv("metadata.csv") 
metadata$rep <- gsub(".*-([0-9])$","\\1",metadata$SampleID)



####====1.0 Make a raw dataset to play with/visualise====####


EUK.data<- read.csv("rawData/EUK.raw.names.csv",row.names = 1)

EUK.tax <- read.csv("taxonomy/EUK.parsed.csv",row.names = 1)

colnames(EUK.data)

EUK.all <- cbind(EUK.data,EUK.tax[match(rownames(EUK.data),EUK.tax$OTU),])

sampleIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(EUK.data))

unique(sampleIndex)

EUK.nReps <- data.frame(matrix(0,nrow = length(EUK.data[,1]),ncol=length(unique(sampleIndex))))

colnames(EUK.nReps) <- unique(sampleIndex)
rownames(EUK.nReps) <- rownames(EUK.data)

EUK.binary <-EUK.data
EUK.binary[EUK.binary<1] <- 0
EUK.binary[EUK.binary>0] <- 1
binaryIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(EUK.binary))

for (column in 1:length(EUK.nReps[1,])){
  EUK.nReps[,column] <- rowSums(EUK.binary[,binaryIndex %in% colnames(EUK.nReps)[column]])
  }

EUK.all.nReps <- cbind(EUK.nReps,EUK.tax[match(rownames(EUK.binary),EUK.tax$OTU),])

write.csv(EUK.all,"rawdata/EUK.raw.wTAX.csv")
write.csv(EUK.all.nReps,"rawdata/EUK.raw.nReps.wTAX.csv")



RIZ.data<- read.csv("rawData/RIZ.raw.names.csv",row.names = 1)

RIZ.tax <- read.csv("taxonomy/RIZ.parsed.csv",row.names = 1)

colnames(RIZ.data)

RIZ.all <- cbind(RIZ.data,RIZ.tax[match(rownames(RIZ.data),RIZ.tax$OTU),])

sampleIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(RIZ.data))

unique(sampleIndex)

RIZ.nReps <- data.frame(matrix(0,nrow = length(RIZ.data[,1]),ncol=length(unique(sampleIndex))))

colnames(RIZ.nReps) <- unique(sampleIndex)
rownames(RIZ.nReps) <- rownames(RIZ.data)

RIZ.binary <-RIZ.data
RIZ.binary[RIZ.binary<1] <- 0
RIZ.binary[RIZ.binary>0] <- 1
binaryIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(RIZ.binary))

for (column in 1:length(RIZ.nReps[1,])){
  if (is.vector(RIZ.binary[,binaryIndex %in% colnames(RIZ.nReps)[column]])){
    RIZ.nReps[,column]  <- RIZ.binary[,binaryIndex %in% colnames(RIZ.nReps)[column]]}else{
  RIZ.nReps[,column] <- rowSums(RIZ.binary[,binaryIndex %in% colnames(RIZ.nReps)[column]])}
}


RIZ.all.nReps <- cbind(RIZ.nReps,RIZ.tax[match(rownames(RIZ.binary),RIZ.tax$OTU),])

write.csv(RIZ.all,"rawdata/RIZ.raw.wTAX.csv")
write.csv(RIZ.all.nReps,"rawdata/RIZ.raw.nReps.wTAX.csv")



####====2.0 Now let's try and do some filtering and make some clean datasets====####

# first lets write this little function to help us collapse replicates 
NrepsMaker <- function(INdataframe,vector){
  ##write these checks
  #check the dataframe is a dataframe
  if(!is.data.frame(INdataframe)){stop("Input dataframe doesn't look like a dataframe")}
  #check the vector is a vector
  if(!is.vector(vector)){stop("Input vector doesn't look like a vector")}
  #check the dataframe contains the vector
  ## TO DO
  #make a new dataframe to captuire the output
  newDataFrame <- data.frame(matrix(0,nrow = length(INdataframe[,1]),ncol=length(unique(vector))))
  #name stuff
  colnames(newDataFrame) <- unique(vector)
  rownames(newDataFrame) <- rownames(INdataframe)
  #make it binary
  INdataframe[INdataframe<1] <- 0
  INdataframe[INdataframe>0] <- 1
  #loop over all the samples with replicates, summing according to the vector
  for (column in 1:length(newDataFrame[1,])){
    # this if statement checks in case there is only one replicate remaining (this sometimes happens in controls)
    if(is.vector(INdataframe[,vector %in% colnames(newDataFrame)[column]])){
      newDataFrame[,column]  <- INdataframe[,vector %in% colnames(newDataFrame)[column]]}else{
        newDataFrame[,column] <- rowSums(INdataframe[,vector %in% colnames(newDataFrame)[column]])}
  }
  return(newDataFrame)
}




# we cut the metadata here becuase it is strangely formatted
metadata <- metadata[0:120,]


for (dataset in list.files("rawdata",pattern="....raw.names.csv")){

datasetname <- substr(dataset,1,3)  
  
indata <- read.csv(paste0("rawdata/",datasetname,".raw.names.csv"),row.names = 1)

expSamples <- indata[,na.omit(match(gsub(",|-",".",sort(metadata$SampleID[metadata$SampleType=="Experimental"])),substring(colnames(indata),5)))]

ctlSamples <-indata[,na.omit(match(gsub(",|-",".",sort(metadata$SampleID[metadata$SampleType=="ControlN"])),substring(colnames(indata),5)))]

#Filter 1 - minimum number of reads for any ID
expSamples[expSamples< minreads] <- 0
expSamples <- expSamples[rowSums(expSamples) > 0,]


#####NOTE - we are not running this filter as we want to be V sensitive! 
#Filter 2 - within samples OTU must appear in more than one sample (this works because there are lots of reps per site and sample)
#filtersam <- expSamples
#filtersam[filtersam>0 ] <- 1
#filtersam <-filtersam[rowSums(filtersam) > 1,]
#expSamples <- expSamples[rownames(expSamples) %in% rownames(filtersam),]

#Filter 3 -Maximum value in neg = 0 value in samples

controlsCONTAM <- ctlSamples[rowSums(ctlSamples) > 0,]
for (contamOTU in 1:length(controlsCONTAM[,1])){
  loopOTU <- row.names(controlsCONTAM[contamOTU,])
  loopMax <- max(as.numeric(controlsCONTAM[contamOTU,]))
  #loopSum <- sum(as.numeric(controlsCONTAM[contamOTU,]))
  if (any(is.na(expSamples[loopOTU,]))){next}
  expSamples[loopOTU,expSamples[loopOTU,]<loopMax] <- 0
  print(paste("Cleaning contaminants",contamOTU))
}


##### make a version of the data with Nreps
expSamplesNreps <- NrepsMaker(expSamples,gsub("(^.*)[.][0-9]$","\\1",colnames(expSamples)))


#Reattach taxonomy and ASVs

rawSeqs <- as.character(readDNAStringSet(paste0("rawdata/ASVs/",datasetname,".DADA2.ASVs.fasta")))
Assignments <- read.csv(paste0("taxonomy/",datasetname,".parsed.csv"),row.names = 1)

CleanedOutput <- cbind(expSamples,
                       unname(rawSeqs)[match(row.names(expSamples),names(rawSeqs))],
                       Assignments[match(row.names(expSamples),Assignments$OTU),])
dir.create("cleanedData",showWarnings = F)
write.csv(CleanedOutput,paste0("cleanedData/clean.",dataset,".csv"))

CleanedNrepsOutput <- cbind(expSamplesNreps,
                       unname(rawSeqs)[match(row.names(expSamplesNreps),names(rawSeqs))],
                       Assignments[match(row.names(expSamplesNreps),Assignments$OTU),])

write.csv(CleanedNrepsOutput,paste0("cleanedData/clean.",dataset,".Nreps.csv"))


}


####====3.0 ====#### BASEMENT

## Lets write a function to collapse a dataframe by a vector so it makes a number of positive reps dataframe


sampleIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(EUK.data))

unique(sampleIndex)

EUK.nReps <- data.frame(matrix(0,nrow = length(EUK.data[,1]),ncol=length(unique(sampleIndex))))

colnames(EUK.nReps) <- unique(sampleIndex)
rownames(EUK.nReps) <- rownames(EUK.data)

EUK.binary <-EUK.data
EUK.binary[EUK.binary<1] <- 0
EUK.binary[EUK.binary>0] <- 1
binaryIndex <- gsub("(^.*)[.][0-9]$","\\1",colnames(EUK.binary))

for (column in 1:length(EUK.nReps[1,])){
  EUK.nReps[,column] <- rowSums(EUK.binary[,binaryIndex %in% colnames(EUK.nReps)[column]])
}

EUK.all.nReps <- cbind(EUK.nReps,EUK.tax[match(rownames(EUK.binary),EUK.tax$OTU),])

write.csv(EUK.all,"rawdata/EUK.raw.wTAX.csv")
write.csv(EUK.all.nReps,"rawdata/EUK.raw.nReps.wTAX.csv")


NrepsMaker <- function(INdataframe,vector){
  ##write these checks
  #check the dataframe is a dataframe
  if(!is.data.frame(INdataframe)){stop("Input dataframe doesn't look like a dataframe")}
  #check the vector is a vector
  if(!is.vector(vector)){stop("Input vector doesn't look like a vector")}
  #check the dataframe contains the vector
  ## TO DO
  #make a new dataframe to captuire the output
  newDataFrame <- data.frame(matrix(0,nrow = length(INdataframe[,1]),ncol=length(unique(vector))))
  #name stuff
  colnames(newDataFrame) <- unique(vector)
  rownames(newDataFrame) <- rownames(INdataframe)
  #make it binary
  INdataframe[INdataframe<1] <- 0
  INdataframe[INdataframe>0] <- 1
  #loop over all the samples with replicates, summing according to the vector
  for (column in 1:length(newDataFrame[1,])){
    # this if statement checks in case there is only one replicate remaining (this sometimes happens in controls)
    if(is.vector(INdataframe[,vector %in% colnames(newDataFrame)[column]])){
      newDataFrame[,column]  <- INdataframe[,vector %in% colnames(newDataFrame)[column]]}else{
        newDataFrame[,column] <- rowSums(INdataframe[,vector %in% colnames(newDataFrame)[column]])}
  }
  return(newDataFrame)
}

RIZ.

test <- NrepsMaker(EUK.data,gsub("(^.*)[.][0-9]$","\\1",colnames(EUK.data)))




