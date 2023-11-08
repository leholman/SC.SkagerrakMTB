#############################################
####= SeaChange Skagerrak Core Comparison=####
####==== Luke E. Holman====25.09.2023====####
#############################################

####====0.0 Packages====####
library(tidyverse)
library(dplyr)
library(maditr)
library(vegan)
library(breakaway)

#### METABARCODING ####

euk <- read.csv("cleanedData/clean.EUK.raw.names.csv.csv",row.names = 1)
euk.Nreps <- read.csv("cleanedData/clean.EUK.raw.names.csv.Nreps.csv",row.names = 1)
dates <-read.csv("metadataAge.csv")
dates$sampleID <- gsub("\\.","-",gsub(".*(MD9.\\d{4})","\\1", dates$Sample))
taxPR2 <-read.csv("taxonomy/tax.PR2.csv",row.names = 1)


#Correct col names
colnames(euk) <- c(gsub("\\.","-",gsub(".*(MD9.\\d{4})","\\1",colnames(euk[,1:88]))),colnames(euk[,89:99]))
colnames(euk.Nreps) <- c(gsub("\\.","-",gsub(".*(MD9.\\d{4})","\\1",colnames(euk.Nreps[,1:11]))),colnames(euk.Nreps[,12:22]))

#Alpha diversity 

pdf("figures/rarefaction.big.pdf",width = 20,height = 13)
rarecurve(t(euk[,1:88]),label=FALSE,step = 1000)
ordilabel(cbind(rowSums(t(euk[,1:88])), specnumber(t(euk[,1:88])))+5, labels=rownames(t(euk[,1:88])),cex=0.5, border = NA,fill=NULL,col="Darkred")
dev.off()
pdf("figures/rarefaction.small.pdf",width = 9,height = 6.5)
rarecurve(t(euk[,1:88]),label=FALSE,step = 1000)
ordilabel(cbind(rowSums(t(euk[,1:88])), specnumber(t(euk[,1:88])))+5, labels=rownames(t(euk[,1:88])),cex=0.5, border = NA,fill=NULL,col="Darkred")
dev.off()

##ASV richness (blind)

euk.Nreps.high.binary.1rep <- euk.Nreps[,1:11]
euk.Nreps.high.binary.1rep[euk.Nreps.high.binary.1rep>0.5] <- 1

euk.Nreps.high.binary.3rep <- euk.Nreps[,1:11]
euk.Nreps.high.binary.3rep[euk.Nreps.high.binary.3rep<3] <- 0
euk.Nreps.high.binary.3rep[euk.Nreps.high.binary.3rep>2.5] <- 1

euk.Nreps.high.binary.8rep <- euk.Nreps[,1:11]
euk.Nreps.high.binary.8rep[euk.Nreps.high.binary.8rep<8] <- 0
euk.Nreps.high.binary.8rep[euk.Nreps.high.binary.8rep>7.5] <- 1




pdf("figures/richness.sample.pdf",width = 8,height = 5)
par(mar=c(6.1,4.1,1.1,1.1),mfrow=c(1,3))
plot(colSums(euk.Nreps.high.binary.1rep),pch=16,cex=1.5,xaxt="n",ylab="ASV Richness",xlab="")
axis(1,at=1:11,labels = colnames(euk.Nreps.high.binary.1rep),las=2,cex.axis=1)
plot(colSums(euk.Nreps.high.binary.3rep),pch=16,cex=1.5,xaxt="n",ylab="ASV Richness",xlab="")
axis(1,at=1:11,labels = colnames(euk.Nreps.high.binary.3rep),las=2,cex.axis=1)
plot(colSums(euk.Nreps.high.binary.8rep),pch=16,cex=1.5,xaxt="n",ylab="ASV Richness",xlab="")
axis(1,at=1:11,labels = colnames(euk.Nreps.high.binary.8rep),las=2,cex.axis=1)
dev.off()

pdf("figures/richness.date.pdf",width = 8,height = 5)
par(mar=c(5.1,4.1,1.1,1.1),mfrow=c(1,3))
plot(dates$Median[match(colnames(euk.Nreps.high.binary.1rep),dates$sampleID)],
     colSums(euk.Nreps.high.binary.1rep),
     pch=16,cex=1.5,ylab="ASV Richness (1rep)",xlab="CalYrBP")
plot(dates$Median[match(colnames(euk.Nreps.high.binary.3rep),dates$sampleID)],
     colSums(euk.Nreps.high.binary.3rep),
     pch=16,cex=1.5,ylab="ASV Richness (3rep)",xlab="CalYrBP")
plot(dates$Median[match(colnames(euk.Nreps.high.binary.8rep),dates$sampleID)],
     colSums(euk.Nreps.high.binary.8rep),
     pch=16,cex=1.5,ylab="ASV Richness (8rep)",xlab="CalYrBP")
dev.off()


### Here we are estimating richness using breakaway which uses this approach - https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12332
###since our metabarcoding data will only ever amplify a fraction of total biodiversity 
###therefore these are useless as absolute estimates of richness but useful as relative measures of richness along the core

#we use the euk dataset as it contains read data 
richEst <- breakaway(euk[,1:88])


richEstimate <- unlist(lapply(richEst,FUN = function(x){x[["estimate"]]}))
richEstimateCIlwr <- unlist(lapply(richEst,FUN = function(x){x[["ci"]][1]}))
richEstimateCIupr <- unlist(lapply(richEst,FUN = function(x){x[["ci"]][2]}))



pdf("figures/richness.freqEst.pdf",width = 8,height = 5)
plot(dates$Median[match(as.factor(substr(names(richEstimate),1,8)),dates$sampleID)],
     richEstimate,pch=16,
     xlab="CalYrBP",
     ylab="EstimatedRichness")
dev.off()


## Are breakaway estimates and raw richness correlated?

compRich <- data.frame("breakaway"=tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8)),
                       "Richness.8rep"=colSums(euk.Nreps.high.binary.8rep),
                       "Richness.3rep"=colSums(euk.Nreps.high.binary.3rep),
                       "Richness.1rep"=colSums(euk.Nreps.high.binary.1rep))

m1 <- lm(colSums(euk.Nreps.high.binary.1rep)~tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8)))
m2 <- lm(colSums(euk.Nreps.high.binary.3rep)~tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8)))
m3 <- lm(colSums(euk.Nreps.high.binary.8rep)~tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8)))



pdf("figures/breakaway-richness.comp.pdf",width = 8,height = 6.5)
plot(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8)),
     colSums(euk.Nreps.high.binary.1rep),pch=16,cex=1.5,col="lightblue",ylim=c(0,2200),xlim=c(0,1100),
     xlab="Breakaway Estimate",ylab="Observed Richness (8rep)")
abline(m1,col="lightblue",lwd=1.5)
points(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8)),
     colSums(euk.Nreps.high.binary.3rep),pch=16,cex=1.5,col="dodgerblue")
abline(m2,col="dodgerblue",lwd=1.5)
points(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8)),
     colSums(euk.Nreps.high.binary.8rep),pch=16,cex=1.5,col="darkblue")
abline(m3,col="darkblue",lwd=1.5)

legend("topleft",border=NA,legend=c("1 rep","3 rep","8 rep"),col=c("lightblue","dodgerblue","darkblue"),pch=16,pt.cex=2)

text(800,1800,round(summary(lm(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8))~colSums(euk.Nreps.high.binary.1rep)))$adj.r.squared,4),col="lightblue",adj=0)
text(800,700,round(summary(lm(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8))~colSums(euk.Nreps.high.binary.3rep)))$adj.r.squared,4),col="dodgerblue",adj=0)
text(800,200,round(summary(lm(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8))~colSums(euk.Nreps.high.binary.8rep)))$adj.r.squared,4),col="darkblue",adj=0)

dev.off()

summary(lm(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8))~colSums(euk.Nreps.high.binary.1rep)))$adj.r.squared
summary(lm(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8))~colSums(euk.Nreps.high.binary.3rep)))$adj.r.squared
summary(lm(tapply(richEstimate,FUN=mean,INDEX = substr(names(richEstimate),1,8))~colSums(euk.Nreps.high.binary.8rep)))$adj.r.squared

##ToDo List
Alpha
-ASV richness
-Bacterial richness
-Non-metazoans, non-bacteria - protists and others
-Metazoan 

Beta (above datasets)

Taxonomic changes

High Conf Metazoan Comparison


Alpha 
Richness of unfiltered reads
-kraken2
-metaPhylan
-Kmer something?

Bacterial richness
-metaPhylan

Non-metazoans (protists) 
-SMAGs mapping 

Metzoa
-nt/refseq mapping


Beta
Raw seqs 
-simka












  





#plot high quality assignments
euk.Nreps.high <- euk.Nreps[!is.na(euk.Nreps$assignmentQual),]
euk.Nreps.high <- euk.Nreps.high[ euk.Nreps.high$assignmentQual=="High-MH"|euk.Nreps.high$assignmentQual=="High",]
euk.Nreps.high <- euk.Nreps.high[euk.Nreps.high$superkingdom !="Bacteria" ,]
euk.Nreps.high <- euk.Nreps.high[euk.Nreps.high$superkingdom !="Archaea" ,]

sort(table(euk.Nreps$phylum[euk.Nreps$assignmentQual=="High" | euk.Nreps$assignmentQual=="High-MH"]))

#plot richness

euk.Nreps.high.binary <- euk.Nreps.high
euk.Nreps.high.binary[euk.Nreps.high.binary<2] <- 0
euk.Nreps.high.binary[euk.Nreps.high.binary>1] <- 1
pdf("figures/richness.pdf",width = 6.5,height = 5)
par(mar=c(5.1,4.1,1.1,1.1))
plot(colSums(euk.Nreps.high.binary[,1:11]),pch=16,cex=1.5,xaxt="n",ylab="ASV Richness",xlab="")
axis(1,at=1:11,labels = colnames(euk.Nreps.high.binary[,1:11]),las=2,cex.axis=1)
dev.off()


#Little beta div plot 
out <- metaMDS(vegdist(t(euk.Nreps.high[1:11]))) 


pdf("figures/betadiv.pdf",height = 6,width = 6)
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1))
text(out$points[,1],out$points[,2]+0.02,labels = metadata$Mean[match(rownames(out$points),metadata$SampleID.2)])
dev.off()

# pull out the following taxa






###### METAGENOMICS #####

metadata$Mean[match(rownames(out$points),metadata$SampleID.2)]


####====0.1 Load data====####
r100.metazoa <- read.csv("rawdata/metaG.r100.metazoa.csv",row.names = 1)


metadata <- read.csv("metadataAge.csv")
metadata$SampleID.2 <- gsub(".*(MD9-\\d{4}).*","\\1",metadata$Sample)

####====1.0 taxonomic subsets ====####

## We want a bubble plot of very open (1 read any damage) and very closed (100 reads damaged/undamaged) taxonomy 

# Lets start with metazoa

inputdata <- r100.metazoa

inputdata$sample2 <- factor(inputdata$sample2,levels= sort(unique(inputdata$sample2),decreasing = TRUE))
inputdata.genus <- inputdata[inputdata$tax_rank=="genus",]
genus.inputdata.wide <- dcast(inputdata.genus, tax_name ~ sample2, value.var="N_reads")

genus.inputdata.wide.m <- as.matrix(genus.inputdata.wide[,2:length(genus.inputdata.wide[1,])])
rownames(genus.inputdata.wide.m) <- genus.inputdata.wide$tax_name 
genus.inputdata.wide.m[is.na(genus.inputdata.wide.m)] <- 0


inputdata.genus$ancient <- rep(1,length(inputdata.genus$sample))
inputdata.genus$ancient[inputdata.genus$D_max>0.1 & inputdata.genus$lambda_LR>2] <- 2


metadata$Mean[match(as.character(inputdata.genus$sample2),metadata$SampleID.2)]


pdf("figures/metazoa.r100.pdf",height = 6,width = 5)
par(mar=c(7.1,7.1,1.1,1.1))
plot(as.numeric(as.factor(inputdata.genus$tax_name)),
     metadata$Mean[match(as.character(inputdata.genus$sample2),metadata$SampleID.2)],
     xlab="",
     ylab="Years Cal BP",
     ylim=c(8500,0),
     xaxt="n",
     frame.plot=TRUE,
     pch=16,
     col=inputdata.genus$ancient,
     cex=as.numeric(as.character(cut(inputdata.genus$N_reads,breaks=c(50,100,200,300,500,1000,3000,100000000),labels=c(0,0.2,0.4,0.8,1.2,1.6,1.8)))))
axis(1,at=1:length(levels(as.factor(inputdata.genus$tax_name))),labels = levels(as.factor(inputdata.genus$tax_name)),las=2,cex.axis=0.5)
dev.off()













