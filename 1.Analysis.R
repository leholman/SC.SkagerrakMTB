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
library(RColorBrewer)
library(scales)

#### METABARCODING ####

euk <- read.csv("cleanedData/clean.EUK.raw.names.csv.csv",row.names = 1)
euk.Nreps <- read.csv("cleanedData/clean.EUK.raw.names.csv.Nreps.csv",row.names = 1)
dates <-read.csv("metadataAge.csv")
dates$sampleID <- gsub("\\.","-",gsub(".*(MD9.\\d{4})","\\1", dates$Sample))
taxPR2 <-read.csv("taxonomy/EUK.tax.PR2.csv",row.names = 1)
euk.css <- read.csv("cleanedData/clean.EUK.CSS.csv",row.names = 1)
euk.rare <- read.csv("cleanedData/clean.EUK.rarefy.csv",row.names = 1)


#Correct col names
colnames(euk) <- c(gsub("\\.","-",gsub(".*(MD9.\\d{4})","\\1",colnames(euk[,1:88]))),colnames(euk[,89:99]))
colnames(euk.Nreps) <- c(gsub("\\.","-",gsub(".*(MD9.\\d{4})","\\1",colnames(euk.Nreps[,1:11]))),colnames(euk.Nreps[,12:22]))

#Taxonomic overview
## Some functions to help 

CountTable <- function(in.taxonomy,in.data,output="Count",some.unassigned=T){
  if(length(in.taxonomy)!=length(in.data[,1])){stop("Dataframe and corresponding taxonomy are not the same length")}
  in.taxonomy[is.na(in.taxonomy)] <- ""
  out.dat <- as.data.frame(matrix(ncol=length(in.data[1,]),nrow=length(unique(in.taxonomy))))
  rownames(out.dat) <- sort(unique(in.taxonomy))
  colnames(out.dat) <- colnames(in.data)    
  out.dat.abundance <- out.dat
  for (sample in 1:length(in.data[1,])){
    out.dat[,sample] <- table(in.taxonomy[in.data[,sample]>0])[match(sort(unique(in.taxonomy)),names(table(in.taxonomy[in.data[,sample]>0])))]
    out.dat.abundance[,sample] <- aggregate(in.data[,sample], by=list(Category=in.taxonomy), FUN=sum)[,2]
  }
  out.dat[is.na(out.dat)] <- 0
  if(some.unassigned==T){rownames(out.dat)[1] <- "Unassigned"}
  if(output=="Count"){return(out.dat)}else if(
    output=="Abundance"){return(out.dat.abundance)}
}
#Then we write a function for concatenating small abundance groups per sample
minAbundance <- function(inputtable=NA,minAbun= 0.01){
  inputtable <- rbind(inputtable,rep(0,dim(inputtable)[1]))
  rownames(inputtable)[dim(inputtable)[1]] <- "Others"
  for (row in 1:dim(inputtable)[2]){
    min <- sum(inputtable[,row])*minAbun
    others <- sum(inputtable[inputtable[,row]<min,row])
    inputtable[inputtable[,row]<min,row] <- 0
    inputtable["Others",row] <- others
    inputtable <- inputtable[rowSums(inputtable)>1,]
  }
  return(inputtable)
}

##function to make dataset binary 

make_binary <- function (df,threshold){
  df <- sapply(df, function(x) ifelse(is.numeric(x) & x < threshold, 0, 1))
  return(as.data.frame(df))
}






##Colour function 
#getPalette = colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette = colorRampPalette(brewer.pal(9, "Set3"))


#First we get rid of poor assignments - these can be changed 
taxPR2.f <- taxPR2[match(row.names(euk),taxPR2$X.1),]
taxPR2.f[is.na(taxPR2.f)] <- ""
hist(taxPR2.f$Subdivision,breaks=1000)
taxPR2.f$tax.Domain[taxPR2.f$Domain<70] <- ""
taxPR2.f$tax.Supergroup[taxPR2.f$Supergroup<40] <- ""
taxPR2.f$tax.Division[taxPR2.f$Division<40] <- ""
taxPR2.f$tax.Family[taxPR2.f$Family<40] <- ""

##Domain level taxonomy 
Domain.A <- minAbundance(CountTable(as.character(taxPR2.f$tax.Domain),euk[,1:88],output = "Abundance"),minAbun=0.01)
Domain.C <- minAbundance(CountTable(as.character(taxPR2.f$tax.Domain),euk[,1:88],output = "Count"),minAbun=0.01)
Domain.A.prop <- as.matrix(prop.table(as.matrix(Domain.A),margin = 2))
Domain.C.prop <- as.matrix(prop.table(as.matrix(Domain.C),margin = 2))
row.names(Domain.A.prop)[1] <- "Unknown"
row.names(Domain.C.prop)[1] <- "Unknown"

pdf("figures/tax.domain.pdf",width = 12,height = 9)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 6.1),xpd=TRUE)
barplot(Domain.A.prop,las=2,cex.names=0.6,col=rev(getPalette(dim(Domain.A.prop)[1])),ylab="Read Abundance")
legend(105,0.8,rev(rownames(Domain.A.prop)),fill=getPalette(dim(Domain.A.prop)[1]),cex=0.9,bty = "n",y.intersp=0.75)
barplot(Domain.C.prop,las=2,cex.names=0.6,col=rev(getPalette(dim(Domain.C.prop)[1])),ylab="ASV Counts")
legend(105,0.8,rev(rownames(Domain.C.prop)),fill=getPalette(dim(Domain.C.prop)[1]),cex=0.9,bty = "n",y.intersp=0.75)
dev.off()

##Supergroup level taxonomy 
Supergroup.A <- minAbundance(CountTable(as.character(taxPR2.f$tax.Supergroup),euk[,1:88],output = "Abundance"),minAbun=0.01)
Supergroup.C <- minAbundance(CountTable(as.character(taxPR2.f$tax.Supergroup),euk[,1:88],output = "Count"),minAbun=0.01)
Supergroup.A.prop <- as.matrix(prop.table(as.matrix(Supergroup.A),margin = 2))
Supergroup.C.prop <- as.matrix(prop.table(as.matrix(Supergroup.C),margin = 2))
row.names(Supergroup.A.prop)[1] <- "Unknown"
row.names(Supergroup.C.prop)[1] <- "Unknown"

pdf("figures/tax.supergroup.pdf",width = 12,height = 9)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 6.1),xpd=TRUE)
barplot(Supergroup.A.prop,las=2,cex.names=0.6,col=rev(getPalette(dim(Supergroup.A.prop)[1])),ylab="Read Abundance")
legend(105,0.8,rev(rownames(Supergroup.A.prop)),fill=getPalette(dim(Supergroup.A.prop)[1]),cex=0.9,bty = "n",y.intersp=0.75)
barplot(Supergroup.C.prop,las=2,cex.names=0.6,col=rev(getPalette(dim(Supergroup.C.prop)[1])),ylab="ASV Counts")
legend(105,0.8,rev(rownames(Supergroup.C.prop)),fill=getPalette(dim(Supergroup.C.prop)[1]),cex=0.9,bty = "n",y.intersp=0.75)
dev.off()


##Division level taxonomy 
Division.A <- minAbundance(CountTable(as.character(taxPR2.f$tax.Division),euk[,1:88],output = "Abundance"),minAbun=0.01)
Division.C <- minAbundance(CountTable(as.character(taxPR2.f$tax.Division),euk[,1:88],output = "Count"),minAbun=0.01)
Division.A.prop <- as.matrix(prop.table(as.matrix(Division.A),margin = 2))
Division.C.prop <- as.matrix(prop.table(as.matrix(Division.C),margin = 2))
row.names(Division.A.prop)[1] <- "Unknown"
row.names(Division.C.prop)[1] <- "Unknown"

pdf("figures/tax.Division.pdf",width = 12,height = 9)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 6.1),xpd=TRUE)
barplot(Division.A.prop,las=2,cex.names=0.6,col=rev(getPalette(dim(Division.A.prop)[1])),ylab="Read Abundance")
legend(105,0.8,rev(rownames(Division.A.prop)),fill=getPalette(dim(Division.A.prop)[1]),cex=0.9,bty = "n",y.intersp=0.75)
barplot(Division.C.prop,las=2,cex.names=0.6,col=rev(getPalette(dim(Division.C.prop)[1])),ylab="ASV Counts")
legend(105,0.8,rev(rownames(Division.C.prop)),fill=getPalette(dim(Division.C.prop)[1]),cex=0.9,bty = "n",y.intersp=0.75)
dev.off()

##Family level taxonomy 
Family.A <- minAbundance(CountTable(as.character(taxPR2.f$tax.Family),euk[,1:88],output = "Abundance"),minAbun=0.01)
Family.C <- minAbundance(CountTable(as.character(taxPR2.f$tax.Family),euk[,1:88],output = "Count"),minAbun=0.01)
Family.A.prop <- as.matrix(prop.table(as.matrix(Family.A),margin = 2))
Family.C.prop <- as.matrix(prop.table(as.matrix(Family.C),margin = 2))
row.names(Family.A.prop)[1] <- "Unknown"
row.names(Family.C.prop)[1] <- "Unknown"

pdf("figures/tax.Family.pdf",width = 12,height = 9)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 6.1),xpd=TRUE)
barplot(Family.A.prop,las=2,cex.names=0.6,col=rev(getPalette(dim(Family.A.prop)[1])),ylab="Read Abundance")
legend(108,1,rev(rownames(Family.A.prop)),fill=getPalette(dim(Family.A.prop)[1]),cex=0.5,bty = "n",y.intersp=0.75)
barplot(Family.C.prop,las=2,cex.names=0.6,col=rev(getPalette(dim(Family.C.prop)[1])),ylab="ASV Counts")
legend(108,1,rev(rownames(Family.C.prop)),fill=getPalette(dim(Family.C.prop)[1]),cex=0.4,bty = "n",y.intersp=0.75)
dev.off()

### Who da fungi!?



test <- euk.Nreps[rownames(euk.Nreps) %in% taxPR2.f$X.1[taxPR2.f$tax.Family=="Eurotiomycetes"],]
test2 <- euk[rownames(euk) %in% taxPR2.f$X.1[taxPR2.f$tax.Family=="Eurotiomycetes"],]
taxPR2.f$X.1[taxPR2.f$tax.Family=="Eurotiomycetes"]
rownames(euk.Nreps) %in% taxPR2.f$X.1[taxPR2.f$tax.Family=="Eurotiomycetes"]


test <- euk.Nreps[rownames(euk.Nreps) %in% taxPR2.f$X.1[taxPR2.f$tax.Family=="Cephaloidophoridae"],]
test2 <- euk[rownames(euk) %in% taxPR2.f$X.1[taxPR2.f$tax.Family=="Cephaloidophoridae"],]
taxPR2.f$X.1[taxPR2.f$tax.Family=="Cephaloidophoridae"]
rownames(euk.Nreps) %in% taxPR2.f$X.1[taxPR2.f$tax.Family=="Cephaloidophoridae"]



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
### our metabarcoding data will only ever amplify a fraction of total biodiversity 
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



## Now let's look at the normalised data 

richEst.css <- breakaway(round(euk.css,0))
richEstimate.css <- unlist(lapply(richEst.css,FUN = function(x){x[["estimate"]]}))
rich.mean.css <- tapply(richEstimate.css,FUN=mean,INDEX = substr(names(richEstimate.css),8,15))
names(rich.mean.css) <- gsub('\\.',"-",names(rich.mean.css))

richEst.rare <- breakaway(euk.rare[,1:72])
richEstimate.rare <- unlist(lapply(richEst.rare,FUN = function(x){x[["estimate"]]}))
rich.mean.rare <- tapply(richEstimate.rare,FUN=mean,INDEX = substr(names(richEstimate.rare),8,15))
names(rich.mean.rare) <- gsub('\\.',"-",names(rich.mean.rare))

#lets plot 

pdf("figures/richness.rarefaction.unscaled.pdf",width = 8,height = 5)
par(mar=c(5.1,4.1,1.1,1.1),mfrow=c(1,3))
plot(dates$Median[match(as.factor(substr(names(richEstimate),1,8)),dates$sampleID)],
     richEstimate,pch=16,
     ylim=c(0,4200),
     xlab="CalYrBP",
     ylab="EstimatedRichness (breakaway on non-norm)")
plot(dates$Median[match(gsub("\\.","-",substr(names(richEstimate.css),8,15)),dates$sampleID)],
     richEstimate.css,pch=16,
     ylim=c(0,4200),
     xlab="CalYrBP",
     ylab="EstimatedRichness CSS norm")
plot(dates$Median[match(gsub("\\.","-",substr(names(richEstimate.rare),8,15)),dates$sampleID)],
     richEstimate.rare,pch=16,
     ylim=c(0,4200),
     xlab="CalYrBP",
     ylab="EstimatedRichness rare norm")
dev.off()

pdf("figures/richness.rarefaction.scaled.pdf",width = 8,height = 5)
par(mar=c(5.1,4.1,1.1,1.1),mfrow=c(1,3))
plot(dates$Median[match(as.factor(substr(names(richEstimate),1,8)),dates$sampleID)],
     richEstimate,pch=16,
     xlab="CalYrBP",
     ylab="EstimatedRichness (breakaway on non-norm)")
plot(dates$Median[match(gsub("\\.","-",substr(names(richEstimate.css),8,15)),dates$sampleID)],
     richEstimate.css,pch=16,
     xlab="CalYrBP",
     ylab="EstimatedRichness CSS norm")
plot(dates$Median[match(gsub("\\.","-",substr(names(richEstimate.rare),8,15)),dates$sampleID)],
     richEstimate.rare,pch=16,
     xlab="CalYrBP",
     ylab="EstimatedRichness rare norm")
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
     xlab="Breakaway Estimate",ylab="Observed Richness")
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


#Now some taxonomic subsets of richness with 3 reps 



euk.Nreps.high.binary.3rep.PRO <- euk.Nreps.high.binary.3rep[taxPR2.f$tax.Domain=="Bacteria",]
euk.Nreps.high.binary.3rep.EUK.1 <- euk.Nreps.high.binary.3rep[taxPR2.f$tax.Domain=="Eukaryota",]
euk.Nreps.high.binary.3rep.EUK <- euk.Nreps.high.binary.3rep.EUK.1[!(rownames(euk.Nreps.high.binary.3rep.EUK.1) %in% taxPR2.f$X.1[taxPR2.f$tax.Subdivision=="Metazoa"]),]
euk.Nreps.high.binary.3rep.MET <- euk.Nreps.high.binary.3rep.EUK.1[rownames(euk.Nreps.high.binary.3rep.EUK.1) %in% taxPR2.f$X.1[taxPR2.f$tax.Subdivision=="Metazoa"],]


pdf("figures/richness.groups.pdf",width = 8,height = 5)
par(mar=c(5.1,4.1,1.1,1.1),mfrow=c(1,3))
plot(dates$Median[match(colnames(euk.Nreps.high.binary.3rep.PRO),dates$sampleID)],
     colSums(euk.Nreps.high.binary.3rep.PRO),
     pch=16,cex=1.5,ylab="ASV Richness (Prokaryotes)",xlab="CalYrBP")
plot(dates$Median[match(colnames(euk.Nreps.high.binary.3rep.EUK),dates$sampleID)],
     colSums(euk.Nreps.high.binary.3rep.EUK),
     pch=16,cex=1.5,ylab="ASV Richness (Protists)",xlab="CalYrBP")
plot(dates$Median[match(colnames(euk.Nreps.high.binary.3rep.MET),dates$sampleID)],
     colSums(euk.Nreps.high.binary.3rep.MET),
     pch=16,cex=1.5,ylab="ASV Richness (Metazoa)",xlab="CalYrBP")
dev.off()





#Beta diversity 


out <- metaMDS(vegdist(t(euk.Nreps.high.binary.3rep[1:11])),trymax = 200) 
out.j <- metaMDS(vegdist(t(euk.Nreps.high.binary.3rep[1:11]),binary = TRUE,method = "jaccard"),trymax = 200) 
pdf("figures/betadiv.BC.3reps.pdf",height = 6,width = 6)
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,ylab="",xlab="",
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1))
text(out$points[,1],out$points[,2]+0.02,labels = dates$Median[match(rownames(out$points),dates$sampleID)])
dev.off()
pdf("figures/betadiv.JC.3reps.pdf",height = 6,width = 6)
plot(out.j$points[,1],out.j$points[,2],col="darkred",cex=1.3,pch=16,ylab="",xlab="",
     ylim=c(min(out.j$points[,2])-0.1,max(out.j$points[,2])+0.1),
     xlim=c(min(out.j$points[,1])-0.1,max(out.j$points[,1])+0.1))
text(out.j$points[,1],out.j$points[,2]+0.02,labels = dates$Median[match(rownames(out.j$points),dates$sampleID)])
dev.off()



out <- metaMDS(vegdist(t(euk[1:88])))
out.j <- metaMDS(vegdist(t(euk[1:88]),binary = TRUE,method = "jaccard"),trymax = 200) 
pdf("figures/betadiv.BC.allreps.pdf",height = 8,width = 8)
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1),ylab="",xlab="")
ordihull(out,substr(colnames(euk[1:88]),1,8),draw = "polygon",col="darkgrey",lty=0)
ordispider(out,substr(colnames(euk[1:88]),1,8),lty=1,lwd=2,col="grey2")
points(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16)
for (i in 1:(length(tapply(out$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8)))-1)) {
         arrows(tapply(out$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))[[i]],
                tapply(out$points[,2],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))[[i]],
                tapply(out$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))[[i+1]],
                tapply(out$points[,2],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))[[i+1]],
                length = 0.1,lwd = 1.5,col = "red3")}
text(tapply(out$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8)),
     tapply(out$points[,2],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))+0.25,
     labels = dates$Median[match(names(tapply(out$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))),dates$sampleID)],
     col="darkblue")
dev.off()

pdf("figures/betadiv.JC.allreps.pdf",height = 8,width = 8)
plot(out.j$points[,1],out.j$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out.j$points[,2])-0.1,max(out.j$points[,2])+0.1),
     xlim=c(min(out.j$points[,1])-0.1,max(out.j$points[,1])+0.1),ylab="",xlab="")
ordihull(out.j,substr(colnames(euk[1:88]),1,8),draw = "polygon",col="darkgrey",lty=0)
ordispider(out.j,substr(colnames(euk[1:88]),1,8),lty=1,lwd=2,col="grey2")
points(out.j$points[,1],out.j$points[,2],col="darkred",cex=1.3,pch=16)
for (i in 1:(length(tapply(out$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8)))-1)) {
  arrows(tapply(out.j$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))[[i]],
         tapply(out.j$points[,2],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))[[i]],
         tapply(out.j$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))[[i+1]],
         tapply(out.j$points[,2],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))[[i+1]],
         length = 0.1,lwd = 1.5,col = "red3")}
text(tapply(out.j$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8)),
     tapply(out.j$points[,2],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))+0.25,
     labels = dates$Median[match(names(tapply(out.j$points[,1],FUN=mean,INDEX = substr(colnames(euk[1:88]),1,8))),dates$sampleID)],
     col="darkblue")
dev.off()

## Pull selected taxa

euk.selectedTaxa <- euk.Nreps[row.names(euk.Nreps) %in% c("ASV_35","ASV_621","ASV_1456","ASV_2468"),]

MTG <-read.csv("rawdata/r100.metazoa.csv",row.names =1)
MTG.genus <- MTG[MTG$tax_rank=="genus",]
MTG.selected <- MTG.genus[MTG.genus$tax_name %in% c("Gadus","Clupea","Oikopleura"),]

MTG.P <- read.csv("rawdata/r100.viridiplantae.csv",row.names = 1)
MTG.P.genus <- MTG.P[MTG.P$tax_rank=="genus",]
MTG.zostera <- MTG.P.genus[MTG.P.genus$tax_name =="Zostera",]


MTG.selected.2 <- rbind(MTG.selected,MTG.zostera)

MTG.wide <- dcast(MTG.selected.2, tax_name ~ sample2, value.var="N_reads")

euk.selectedTaxa.long <- melt(euk.selectedTaxa)

euk.selectedTaxa.long$variable

MTG.selected.2$sample2

euk.selectedTaxa.long$OTU

ASVid <- data.frame("ASV"=c("ASV_35","ASV_621","ASV_1456","ASV_2468" ),
                    "ID"=c("Zostera","Oikopleura","Gadus","Clupea"))
euk.selectedTaxa.long$ID <- ASVid$ID[match(euk.selectedTaxa.long$OTU,ASVid$ASV)]


MTB.MTG.comp <- data.frame("Sample"=c(as.character(euk.selectedTaxa.long$variable),MTG.selected.2$sample2),
                           "Value"=c(euk.selectedTaxa.long$value,cut(MTG.selected.2$N_reads,breaks = c(50,100,1000,5000,10000),labels=c(0.5,1,2.5,4))),
                           "ID"=c(euk.selectedTaxa.long$ID,MTG.selected.2$tax_name),
                           "Dataset"=c(rep("MTB",length(euk.selectedTaxa.long$variable)),
                                       rep("MTG",length(MTG.selected.2$sample2))))

MTB.com <- MTB.MTG.comp[MTB.MTG.comp$Dataset=="MTB",]
MTG.com <- MTB.MTG.comp[MTB.MTG.comp$Dataset=="MTG",]

pdf("figures/ComparisonMTBMTG.pdf",height = 4,width = 11)
par(mar=c(5.1, 7.1, 2.1, 9.1), xpd=TRUE)
plot(dates$Median[match(MTB.com$Sample,dates$sampleID)],
     as.numeric(factor(MTB.com$ID,levels=c("Zostera" ,"Gadus","Oikopleura","Clupea")))+0.1,
     pch=16,cex=MTB.com$Value/2,col="dodgerblue",
     ylim=c(0.5,4.5),
     xlab="",yaxt='n',
     ylab="")

points(dates$Median[match(MTG.com$Sample,dates$sampleID)],
       as.numeric(factor(MTG.com$ID,levels=c("Zostera" ,"Gadus","Oikopleura","Clupea")))-0.1,
       pch=16,cex=MTG.com$Value,col="darkred")

axis(2,labels=c("Zostera" ,"Gadus","Oikopleura","Clupea"),1:4,las=1)

legend(9000,4.5,col = "dodgerblue",pch=16,
       pt.cex=c(0.5,1.5,4),legend=c("  1 rep","  3 reps","  8 reps"),bty="n",y.intersp=1.5)
legend(9000,2.5,col = "darkred",pch=16,
       pt.cex=c(2,3,4),legend=c(" 100-1k reads"," 1k-5k reads"," 5k+ reads"),bty="n",y.intersp=1.5)


dev.off()

##Comparison plots

## import new MTG

## Dataset 1 all taxa @ genus
## Dataset 2 metazoa @ genus
## Dataset 3 age-dmg excluded all taxa
## Dataset 4 age-dmg excluded metazoa


#DS1
MTG.raw.DS1 <- read.csv("rawdata/Skagerrak_0_38.csv.gz")
MTG.raw.DS1 <- MTG.raw.DS1[MTG.raw.DS1$tax_rank=="genus",]
MTG.raw.DS1 <- MTG.raw.DS1[MTG.raw.DS1$N_reads>99,]
MTG.raw.DS1$sample2 <- gsub(".*(MD9-\\d{4}).*","\\1",MTG.raw.DS1$sample)
#DS2
MTG.raw.DS2 <- MTG.raw.DS1[grepl("metazoa",MTG.raw.DS1$tax_path),] 
#DS3
MTG.raw.DS3 <- read.csv('rawdata/Giulia080124/filtered_genus_all_2nc.csv',row.names = 1)
MTG.raw.DS3$sample2 <- gsub(".*(MD9-\\d{4}).*","\\1",MTG.raw.DS3$sample)
MTG.raw.DS3$Kingdom[ is.na(MTG.raw.DS3$Kingdom)] <- ""
MTG.raw.DS3 <- MTG.raw.DS3[MTG.raw.DS3$Kingdom!="Viridiplantae",]
MTG.raw.DS3 <- MTG.raw.DS3[MTG.raw.DS3$Filtering=="Ancient",]

#DS4
MTG.raw.DS4 <- read.csv('rawdata/Giulia080124/filtered_metazoans_all_2nc.csv',row.names = 1)
MTG.raw.DS4$sample2 <- gsub(".*(MD9-\\d{4}).*","\\1",MTG.raw.DS4$sample)
MTG.raw.DS4 <- MTG.raw.DS4[MTG.raw.DS4$Filtering=="Ancient",]

## Make wide data
MTG.wide.DS1 <- dcast(MTG.raw.DS1, tax_name ~ sample2, value.var="N_reads",fill = 0)
MTG.wide.DS2 <- dcast(MTG.raw.DS2, tax_name ~ sample2, value.var="N_reads",fill = 0)
MTG.wide.DS3 <- dcast(MTG.raw.DS3, tax_name ~ sample2, value.var="N_reads",fill = 0)
MTG.wide.DS4 <- dcast(MTG.raw.DS4, tax_name ~ sample2, value.var="N_reads",fill = 0)


##MTB
##DS1 all ASVs
##DS2 met ASVs

MTB.wide.DS1 <- euk.Nreps[,0:11]
MTB.wide.DS2 <- euk.Nreps[rownames(euk.Nreps) %in% taxPR2.f$X.1[taxPR2.f$tax.Subdivision=="Metazoa"],0:11]


##Alpha 


MTG.binary.DS1 <- make_binary(MTG.wide.DS1,1)[,-1]
MTG.binary.DS1 <- MTG.binary.DS1[,!grepl("NTC|BLANK",names(MTG.binary.DS1))]
MTG.binary.DS2 <- make_binary(MTG.wide.DS2,1)[,-1]
MTG.binary.DS3 <- make_binary(MTG.wide.DS3,1)[,-1]
MTG.binary.DS4 <- make_binary(MTG.wide.DS4,1)[,-1]

MTB.binary.DS1 <- make_binary(MTB.wide.DS1,1)
MTB.binary.DS2 <- make_binary(MTB.wide.DS2,1)





## Plots 

pdf("figures/figure1/richness.MTG.DS1.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(dates$Median[match(as.factor(substr(names(MTG.binary.DS1),1,8)),dates$sampleID)],
     colSums(MTG.binary.DS1),
     pch=16,
     xlab="CalYrBP",
     ylab="MTG.DS1 Genus Richness",
     xlim=c(0,8400))
dev.off()


pdf("figures/figure1/richness.MTG.DS2.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(dates$Median[match(as.factor(substr(names(MTG.binary.DS2),1,8)),dates$sampleID)],
     colSums(MTG.binary.DS2),
     pch=16,
     xlab="CalYrBP",
     ylab="MTG.DS2 Genus Richness",
     xlim=c(0,8400))
dev.off()

pdf("figures/figure1/richness.MTG.DS3.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(dates$Median[match(as.factor(substr(names(MTG.binary.DS3),1,8)),dates$sampleID)],
     colSums(MTG.binary.DS3),
     pch=16,
     xlab="CalYrBP",
     ylab="MTG.DS3 Genus Richness",
     xlim=c(0,8400))
dev.off()

pdf("figures/figure1/richness.MTG.DS4.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(dates$Median[match(as.factor(substr(names(MTG.binary.DS4),1,8)),dates$sampleID)],
     colSums(MTG.binary.DS4),
     pch=16,
     xlab="CalYrBP",
     ylab="MTG.DS4 Genus Richness",
     xlim=c(0,8400))
dev.off()

pdf("figures/figure1/richness.MTB.DS1.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(dates$Median[match(as.factor(substr(names(MTB.binary.DS1),1,8)),dates$sampleID)],
     colSums(MTB.binary.DS1),
     pch=16,
     xlab="CalYrBP",
     ylab="MTB.DS1 ASV Richness",
     xlim=c(0,8400))
dev.off()

pdf("figures/figure1/richness.MTB.DS2.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(dates$Median[match(as.factor(substr(names(MTB.binary.DS2),1,8)),dates$sampleID)],
     colSums(MTB.binary.DS2),
     pch=16,
     xlab="CalYrBP",
     ylab="MTB.DS2 ASV Richness",
     xlim=c(0,8400))
dev.off()



colSums(MTG.binary.DS1)
colSums(MTB.binary.DS1)

pdf("figures/figure1/richness.comp.MTG.DS1-MTB.DS1.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(colSums(MTG.binary.DS1),
     colSums(MTB.binary.DS1),
     pch=16,
     xlab="MTG.DS1",
     ylab="MTB.DS1")

test <- cor.test(colSums(MTG.binary.DS1),colSums(MTB.binary.DS1))
abline(lm(colSums(MTB.binary.DS1)~colSums(MTG.binary.DS1)),col="red")
legend("topright",legend=c(paste0("Corr = ",round(unname(test$estimate),3)),
                           paste0("p =",signif(test$p.value,3))),
       text.col="red",bty="n")

dev.off()


pdf("figures/figure1/richness.comp.MTG.DS1-MTB.DS2.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(colSums(MTG.binary.DS1),
     colSums(MTB.binary.DS2),
     pch=16,
     xlab="MTG.DS1",
     ylab="MTB.DS2")

test <- cor.test(colSums(MTG.binary.DS1),colSums(MTB.binary.DS2))
abline(lm(colSums(MTB.binary.DS2)~colSums(MTG.binary.DS1)),col="red")
legend("topright",legend=c(paste0("Corr = ",round(unname(test$estimate),3)),
                           paste0("p =",signif(test$p.value,3))),
       text.col="red",bty="n")

dev.off()

#here we make a little index to subset the MTB dataset as MTG loses a sample

index <- match(names(colSums(MTG.binary.DS2)),names(colSums(MTB.binary.DS1)))

pdf("figures/figure1/richness.comp.MTG.DS2-MTB.DS1.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(colSums(MTG.binary.DS2),
     colSums(MTB.binary.DS1)[index],
     pch=16,
     xlab="MTG.DS2",
     ylab="MTB.DS1")

test <- cor.test(colSums(MTG.binary.DS2),colSums(MTB.binary.DS1)[index])
abline(lm(colSums(MTB.binary.DS1)[index]~colSums(MTG.binary.DS2)),col="red")
legend("topright",legend=c(paste0("Corr = ",round(unname(test$estimate),3)),
                           paste0("p =",signif(test$p.value,3))),
       text.col="red",bty="n")

dev.off()


index <- match(names(colSums(MTG.binary.DS2)),names(colSums(MTB.binary.DS2)))

pdf("figures/figure1/richness.comp.MTG.DS2-MTB.DS2.pdf",width = 4,height = 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(colSums(MTG.binary.DS2),
     colSums(MTB.binary.DS2)[index],
     pch=16,
     xlab="MTG.DS2",
     ylab="MTB.DS2")

test <- cor.test(colSums(MTG.binary.DS2),colSums(MTB.binary.DS2)[index])
abline(lm(colSums(MTB.binary.DS2)[index]~colSums(MTG.binary.DS2)),col="red")
legend("topright",legend=c(paste0("Corr = ",round(unname(test$estimate),3)),
                           paste0("p =",signif(test$p.value,3))),
       text.col="red",bty="n")

dev.off()

#### BETA DIV 
##MTG
MTG.wide.DS1 <- dcast(MTG.raw.DS1, tax_name ~ sample2, value.var="N_reads",fill = 0)
MTG.wide.DS2 <- dcast(MTG.raw.DS2, tax_name ~ sample2, value.var="N_reads",fill = 0)
MTG.wide.DS3 <- dcast(MTG.raw.DS3, tax_name ~ sample2, value.var="N_reads",fill = 0)
MTG.wide.DS4 <- dcast(MTG.raw.DS4, tax_name ~ sample2, value.var="N_reads",fill = 0)

##MTB
MTB.wide.DS1 <- euk.Nreps[,0:11]
MTB.wide.DS2 <- euk.Nreps[rownames(euk.Nreps) %in% taxPR2.f$X.1[taxPR2.f$tax.Subdivision=="Metazoa"],0:11]


##MTG DS1

datain <- MTG.wide.DS1[,4:14]
colnames(datain) <- as.character(dates$Median[match(colnames(datain),dates$sampleID)])
out <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2))),trymax = 200)
out.j <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2)),binary = TRUE,method = "jaccard"),trymax = 200) 

pdf("figures/figure1/beta.MTG.DS1.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1),ylab="",xlab="")
for (i in 1:(length(out$points[,1])-1)) {
  arrows(out$points[i,1],
         out$points[i,2],
         out$points[i+1,1],
         out$points[i+1,2],
         length = 0.1,lwd = 1.5,col = "red3")}
points(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,)
text(out$points[,1],
     out$points[,2]+0.1,
     labels = rownames(out$points),
     col="darkblue")
dev.off()
MTG.DS1.nMDS <- out

##MTG DS2
datain <- MTG.wide.DS2[,-1]
colnames(datain) <- as.character(dates$Median[match(colnames(datain),dates$sampleID)])
out <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2))),trymax = 200)
out.j <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2)),binary = TRUE,method = "jaccard"),trymax = 200) 

pdf("figures/figure1/beta.MTG.DS2.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1),ylab="",xlab="")
for (i in 1:(length(out$points[,1])-1)) {
  arrows(out$points[i,1],
         out$points[i,2],
         out$points[i+1,1],
         out$points[i+1,2],
         length = 0.1,lwd = 1.5,col = "red3")}
points(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,)
text(out$points[,1],
     out$points[,2]+0.1,
     labels = rownames(out$points),
     col="darkblue")
dev.off()
MTG.DS2.nMDS <- out

##MTG DS3
datain <- MTG.wide.DS3[,-1]
colnames(datain) <- as.character(dates$Median[match(colnames(datain),dates$sampleID)])
out <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2))),trymax = 200)
out.j <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2)),binary = TRUE,method = "jaccard"),trymax = 200) 

pdf("figures/figure1/beta.MTG.DS3.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1),ylab="",xlab="")
for (i in 1:(length(out$points[,1])-1)) {
  arrows(out$points[i,1],
         out$points[i,2],
         out$points[i+1,1],
         out$points[i+1,2],
         length = 0.1,lwd = 1.5,col = "red3")}
points(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,)
text(out$points[,1],
     out$points[,2]+0.1,
     labels = rownames(out$points),
     col="darkblue")
dev.off()
MTG.DS3.nMDS <- out

##DS4
datain <- MTG.wide.DS4[,-1]
colnames(datain) <- as.character(dates$Median[match(colnames(datain),dates$sampleID)])
out <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2))),trymax = 200)
out.j <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2)),binary = TRUE,method = "jaccard"),trymax = 200) 

pdf("figures/figure1/beta.MTG.DS4.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1),ylab="",xlab="")
for (i in 1:(length(out$points[,1])-1)) {
  arrows(out$points[i,1],
         out$points[i,2],
         out$points[i+1,1],
         out$points[i+1,2],
         length = 0.1,lwd = 1.5,col = "red3")}
points(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,)
text(out$points[,1],
     out$points[,2]+0.1,
     labels = rownames(out$points),
     col="darkblue")
dev.off()
MTG.DS4.nMDS <- out



##MTB DS1
datain <- MTB.wide.DS1
colnames(datain) <- as.character(dates$Median[match(colnames(datain),dates$sampleID)])
out <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2))),trymax = 200)
out.j <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2)),binary = TRUE,method = "jaccard"),trymax = 200) 


pdf("figures/figure1/beta.MTB.DS1.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1),ylab="",xlab="")
for (i in 1:(length(out$points[,1])-1)) {
  arrows(out$points[i,1],
         out$points[i,2],
         out$points[i+1,1],
         out$points[i+1,2],
         length = 0.1,lwd = 1.5,col = "red3")}
points(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,)
text(out$points[,1],
     out$points[,2]+0.1,
     labels = rownames(out$points),
     col="darkblue")
dev.off()
MTB.DS1.nMDS <- out



##MTB DS2
datain <- MTB.wide.DS2
colnames(datain) <- as.character(dates$Median[match(colnames(datain),dates$sampleID)])
out <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2))),trymax = 200)
out.j <- metaMDS(vegdist(t(prop.table(as.matrix(datain),2)),binary = TRUE,method = "jaccard"),trymax = 200) 


pdf("figures/figure1/beta.MTB.DS2.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,
     ylim=c(min(out$points[,2])-0.1,max(out$points[,2])+0.1),
     xlim=c(min(out$points[,1])-0.1,max(out$points[,1])+0.1),ylab="",xlab="")
for (i in 1:(length(out$points[,1])-1)) {
  arrows(out$points[i,1],
         out$points[i,2],
         out$points[i+1,1],
         out$points[i+1,2],
         length = 0.1,lwd = 1.5,col = "red3")}
points(out$points[,1],out$points[,2],col="darkred",cex=1.3,pch=16,)
text(out$points[,1],
     out$points[,2]+0.1,
     labels = rownames(out$points),
     col="darkblue")
dev.off()
MTB.DS2.nMDS <- out

##Procrustes 

#MTB.DS1 vs MTG.DS1
proc <- procrustes(MTB.DS1.nMDS$points,MTG.DS1.nMDS$points,symmetric = TRUE,scale = TRUE)
summary(proc)
proc.t <- protest(MTB.DS1.nMDS$points,MTG.DS1.nMDS$points, scores = "sites", permutations = 9999)

pdf("figures/figure1/beta.MTB.DS1-MTG.DS1.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(proc,type = "text",cex = 0.4)
legend("topleft",legend=c(paste0("SS = ",round(proc.t$ss,3)),
                          paste0("p = ",round(proc.t$signif,3))),text.col = "red",bty="n")
dev.off()

#MTB.DS1 vs MTG.DS2
proc <- procrustes(MTB.DS1.nMDS$points[2:11,],MTG.DS2.nMDS$points,symmetric = TRUE,scale = TRUE)
summary(proc)
proc.t <- protest(MTB.DS1.nMDS$points[2:11,],MTG.DS2.nMDS$points, scores = "sites", permutations = 9999)

pdf("figures/figure1/beta.MTB.DS1-MTG.DS2.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(proc,type = "text",cex = 0.4)
legend("topleft",legend=c(paste0("SS = ",round(proc.t$ss,3)),
                          paste0("p = ",round(proc.t$signif,3))),text.col = "red",bty="n")
dev.off()

#MTB.DS2 vs MTG.DS1
proc <- procrustes(MTB.DS2.nMDS$points,MTG.DS1.nMDS$points,symmetric = TRUE,scale = TRUE)
summary(proc)
proc.t <- protest(MTB.DS2.nMDS$points,MTG.DS1.nMDS$points, scores = "sites", permutations = 9999)

pdf("figures/figure1/beta.MTB.DS2-MTG.DS1.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(proc,type = "text",cex = 0.4)
legend("topleft",legend=c(paste0("SS = ",round(proc.t$ss,3)),
                          paste0("p = ",round(proc.t$signif,3))),text.col = "red",bty="n")
dev.off()

#MTB.DS2 vs MTG.DS2
proc <- procrustes(MTB.DS2.nMDS$points[2:11,],MTG.DS2.nMDS$points,symmetric = TRUE,scale = TRUE)
summary(proc)
proc.t <- protest(MTB.DS2.nMDS$points[2:11,],MTG.DS2.nMDS$points, scores = "sites", permutations = 9999)

pdf("figures/figure1/beta.MTB.DS2-MTG.DS2.pdf",height = 5,width = 5)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(proc,type = "text",cex = 0.4)
legend("topleft",legend=c(paste0("SS = ",round(proc.t$ss,3)),
                          paste0("p = ",round(proc.t$signif,3))),text.col = "red",bty="n")
dev.off()

### new plot comparing selected taxa

## Pull selected taxa

euk.selectedTaxa <- euk.Nreps[row.names(euk.Nreps) %in% c("ASV_35","ASV_621","ASV_1456","ASV_2468"),]

MTG.selected <- MTG.raw.DS1[MTG.raw.DS1$tax_name %in% c("Gadus","Clupea","Oikopleura","Zostera"),]

MTG..selected.wide <- dcast(MTG.selected.2, tax_name ~ sample2, value.var="N_reads",fill = 0)

euk.selectedTaxa.long <- melt(euk.selectedTaxa)

euk.selectedTaxa.long$variable

MTG.selected.2$sample2

euk.selectedTaxa.long$OTU

ASVid <- data.frame("ASV"=c("ASV_35","ASV_621","ASV_1456","ASV_2468" ),
                    "ID"=c("Zostera","Oikopleura","Gadus","Clupea"))
euk.selectedTaxa.long$ID <- ASVid$ID[match(euk.selectedTaxa.long$OTU,ASVid$ASV)]


MTB.MTG.comp <- data.frame("Sample"=c(as.character(euk.selectedTaxa.long$variable),MTG.selected.2$sample2),
                           "Value"=c(euk.selectedTaxa.long$value,as.numeric(as.character(cut(MTG.selected.2$N_reads,breaks = c(50,100,1000,5000,10000),labels=c(0.5,1,2.5,4))))),
                           "ID"=c(euk.selectedTaxa.long$ID,MTG.selected.2$tax_name),
                           "Dataset"=c(rep("MTB",length(euk.selectedTaxa.long$variable)),
                                       rep("MTG",length(MTG.selected.2$sample2))))

MTB.com <- MTB.MTG.comp[MTB.MTG.comp$Dataset=="MTB",]
MTG.com <- MTB.MTG.comp[MTB.MTG.comp$Dataset=="MTG",]

#ancient data only 
MTG.raw.DS3. <- read.csv('rawdata/Giulia080124/filtered_genus_all_2nc.csv',row.names = 1)
MTG.raw.DS3$sample2 <- gsub(".*(MD9-\\d{4}).*","\\1",MTG.raw.DS3$sample)
MTG.raw.DS3$Kingdom[ is.na(MTG.raw.DS3$Kingdom)] <- ""
MTG.raw.DS3 <- MTG.raw.DS3[MTG.raw.DS3$Filtering=="Ancient",]
MTG.raw.DS3.s <- MTG.raw.DS3[MTG.raw.DS3$tax_name %in% c("Gadus","Clupea","Oikopleura","Zostera"),]


MTG.com.a <- dcast(MTG.raw.DS3.s, tax_name ~ sample2, value.var="N_reads",fill = 0)
MTG.com.a.plot <- melt(MTG.com.a)
MTG.com.a.plot$values <- cut(MTG.com.a.plot$value,breaks = c(50,100,1000,5000,10000),labels=c(0.5,1,2.5,4))



pdf("figures/figure1/ComparisonMTBMTG.pdf",height = 4,width = 11)
par(mar=c(5.1, 7.1, 2.1, 9.1), xpd=TRUE)
plot(dates$Median[match(MTB.com$Sample,dates$sampleID)],
     as.numeric(factor(MTB.com$ID,levels=c("Zostera" ,"Gadus","Oikopleura","Clupea")))+0.1,
     pch=16,cex=MTB.com$Value/2,col="dodgerblue",
     ylim=c(0.5,4.5),
     xlab="",yaxt='n',
     ylab="")

points(dates$Median[match(MTG.com$Sample,dates$sampleID)],
       as.numeric(factor(MTG.com$ID,levels=c("Zostera" ,"Gadus","Oikopleura","Clupea")))-0.1,
       pch=16,cex=MTG.com$Value,col="pink")
points(dates$Median[match(MTG.com.a.plot$variable,dates$sampleID)],
       as.numeric(factor(MTG.com.a.plot$tax_name,levels=c("Zostera" ,"Gadus","Oikopleura","Clupea")))-0.1,
       pch=16,cex=as.numeric(as.character(MTG.com.a.plot$values)),col="darkred")

axis(2,labels=c("Zostera" ,"Gadus","Oikopleura","Clupea"),1:4,las=1)

legend(9000,4.5,col = "dodgerblue",pch=16,
       pt.cex=c(0.5,1.5,4),legend=c("  1 rep","  3 reps","  8 reps"),bty="n",y.intersp=1.5)
legend(9000,2.5,col = "darkred",pch=16,
       pt.cex=c(2,3,4),legend=c(" 100-1k reads"," 1k-5k reads"," 5k+ reads"),bty="n",y.intersp=1.5)

dev.off()




### Extra stuff
#length age relationship

euk$ASV_len <- nchar(euk$unname.rawSeqs..match.row.names.expSamples...names.rawSeqs...) 

test <- melt(euk, measure.vars=1:88, variable.name="Sample", value.name="nReads")
test2 <- test[test$nReads>5,]

pdf("figures/lengthAgeAllReps.pdf",width=7,height=5.5)
plot(jitter(dates$Median[match(substr(test2$Sample,1,8),dates$sampleID)],amount = 200),
     jitter(test2$ASV_len,amount = 0.3),pch=16,cex=0.3,
     xlab="CalYrBP",
     ylab="ASV Length")
dev.off()

euk.Nreps$ASV_len <- nchar(euk.Nreps$unname.rawSeqs..match.row.names.expSamplesNreps...names.rawSeqs...) 

test3 <- melt(euk.Nreps, measure.vars=1:11, variable.name="Sample", value.name="nReps")
test4 <- test3[test3$nReps>0,]
test5 <- test4[test4$phylum=="Chordata",]

pdf("figures/lengthAge.nReps.pdf",width=11,height=7)
plot(jitter(dates$Median[match(test4$Sample,dates$sampleID)],amount = 150),
     jitter(test4$ASV_len,amount = 0.15),pch=16,cex=0.3,
     xlab="CalYrBP",
     ylab="ASV Length")
dev.off()


pdf("figures/lengthAge.nReps.cols.pdf",width=11,height=7)
par(mar=c(5.1, 4.1, 2.1, 4.1), xpd=TRUE)
palette(c("grey","blue","darkgreen","darkgrey","red"))
plot(jitter(dates$Median[match(test4$Sample,dates$sampleID)],amount = 150),
     jitter(test4$ASV_len,amount = 0.15),pch=16,cex=0.3,
     xlab="CalYrBP",
     ylab="ASV Length",
     col=as.factor(test4$superkingdom))
points(jitter(dates$Median[match(test5$Sample,dates$sampleID)],amount = 150),
       jitter(test5$ASV_len,amount = 0.15),pch=16,cex=0.5,col="purple")
legend(9000,130,legend = c("Eukaryotes","Bacteria","Archaea","Chordata","No Assign"),pch=15,col=c("red","darkgreen","blue","purple","grey"),cex=0.5, bty='n')
dev.off()



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

ordihull(out,substr(names(euk[,1:88]),1,8))







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













