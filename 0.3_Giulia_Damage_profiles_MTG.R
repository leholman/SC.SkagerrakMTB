######Damage profiles MTG ######

#packages
library(ggridges)
library(viridis)
library(tidyverse) 
require(gplots)
library(scales)
library(reshape2)
library(vegan)
library(rioja)
library(readxl)
library(ggplot2)
library(lvplot)
library(dplyr)
library(gghighlight)
library(ggpubr)
library(plotly)
library(gridExtra)


# import data (metaDMG output)
dt <- read_csv("rawdata/Skagerrak_0_38.csv.gz")
metadataAge <- read_csv("MTG.metadataAge.csv")
dt$Age <- metadataAge$Median[match(gsub("^[^-]+-[^-]+-([^_]+)_.*$", "\\1", dt$sample), metadataAge$Sample)]
dt$Age[is.na(dt$Age)] <- "Control" #replacing NAs in blanks with 'Control'

##adding age v2
dt$Median_age <-metadataAge$Median[match(gsub("^[^-]+-[^-]+-([^_]+)_.*$", "\\1", dt$sample), metadataAge$Sample)]
dt$Median_age[is.na(dt$Median_age)] <- "Control" #replacing NAs in blanks with 'Control'



#functions
# mode function 2 returns the mean in case there’s less than 2 values, NAs, and uses the density function
estimate_mode <- function(x) {
  if (length(x) == 2) {
    mean(x)
  } else if (length(x) == 1){
    x
  }else if (length(x) == 0){
    NA
  }else {
    d <- density(x)
    d$x[which.max(d$y)]
  }
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]}

#########Building the damage model on Viridiplantae#########
#Maximum amount of damage 
DamMin = 0.00
#MAP_significance
MS0 = 0
# Minimum reads for parsing taxa
MinRead = 10
# Minimum mean read length
MinLength = 35

#subsetting the table
dt1 <- dt %>% filter(MAP_damage > DamMin, N_reads >= MinRead, mean_L > MinLength, MAP_significance  > MS0,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("MD9", sample))

# calculating stats
dt4 <- dt1 %>%
  group_by(sample) %>%
  summarise(mode = tapply(MAP_damage,sample, getmode))

dt4a <- dt1 %>%
  group_by(sample) %>%
  summarise(A_mode = tapply(MAP_damage,sample, estimate_mode))

dt5 <- dt1 %>%
  group_by(sample) %>%
  summarise(min = tapply(MAP_damage,sample, min))

dt6 <- dt1 %>%
  group_by(sample) %>%
  summarise(max = tapply(MAP_damage,sample, max))

dt7 <- dt1 %>%
  group_by(sample) %>%
  summarise(median = tapply(MAP_damage,sample, median))

dt8 <- dt1 %>%
  group_by(sample) %>%
  summarise(sd = tapply(MAP_damage,sample, sd))

# binding them to table, and adding ages 
dt9 <- cbind(dt4, dt4a[,2], dt5[,2], dt6[,2], dt7[,2], dt8[,2])

#adding age
##
dt9$Age <-metadataAge$Median[match(gsub("^[^-]+-[^-]+-([^_]+)_.*$", "\\1", dt9$sample), metadataAge$Sample)]

#preliminary plot age vs. max damage
dt9 <-dt9[order(dt9$Age),]
plot(dt9$Age, dt9$max)

#checking the classes of the rows and changing these to be numeric
sapply(dt9, class)
dt9[, 2:8] <- lapply(dt9[, 2:8], as.numeric)

# plotting the calculated stats
pdf(file="figures/DNAdamagePlantsReads10.pdf")
plot(dt9$median, dt9$Age, type ="b", ylim = rev(range(c(0,9000))), xlim=c(0,0.3), xlab="DNA damage", ylab="Age in Years BP")
points(dt9$min,dt9$Age, col="red", pch="*")
lines(dt9$min,dt9$Age, col="red",lty=2)
points(dt9$max,dt9$Age, col="blue", pch="+")
lines(dt9$max,dt9$Age, col="blue",lty=2)
points(dt9$mode,dt9$Age, col="grey", pch=1)
lines(dt9$mode,dt9$Age, col="grey",lty=2)
points(dt9$A_mode,dt9$Age, col="green", pch=1)
lines(dt9$A_mode,dt9$Age, col="green",lty=2)
legend(0.2,0,legend=c("Min","Median","Max", "Mode", "A_mode"), col=c("red","black","blue", "grey", "green"),
       pch=c("o","*","+"),lty=c(1,2,3), ncol=1)
dev.off()

######### plotting minimum thresholds
#calculating stats
df4 <-dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 10) %>%
  summarise('10' = tapply(MAP_damage,sample, min))
df4$cat <-"Minimum"

df4a <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 50) %>%
  summarise('50' = tapply(MAP_damage,sample, min))
df4a$cat <-"Minimum"

df5 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 75) %>%
  summarise('75' = tapply(MAP_damage,sample, min))
df5$cat <-"Minimum"

df6 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 100) %>%
  summarise('100' = tapply(MAP_damage,sample, min))
df6$cat <-"Minimum"

df7 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 125) %>%
  summarise('125' = tapply(MAP_damage,sample, min))
df7$cat <-"Minimum"

df8 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 150) %>%
  summarise('150' = tapply(MAP_damage,sample, min))
df8$cat <-"Minimum"

df9 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 300) %>%
  summarise('300' = tapply(MAP_damage,sample, min))
df9$cat <-"Minimum"

df10 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 350) %>%
  summarise('350' = tapply(MAP_damage,sample, min))
df10$cat <-"Minimum"

df11 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 375) %>%
  summarise('375' = tapply(MAP_damage,sample, min))
df11$cat <-"Minimum"


df12 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 500) %>%
  summarise('500' = tapply(MAP_damage,sample, min))
df12$cat <-"Minimum"

df13 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 1000) %>%
  summarise('1000' = tapply(MAP_damage,sample, min))
df13$cat <-"Minimum"

df14 <-dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 10) %>%
  summarise('10' = tapply(MAP_damage,sample, max))
df14$cat <-"Maximum"

df15 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 50) %>%
  summarise('50' = tapply(MAP_damage,sample, max))
df15$cat <-"Maximum"

df16 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 75) %>%
  summarise('75' = tapply(MAP_damage,sample, max))
df16$cat <-"Maximum"

df17 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 100) %>%
  summarise('100' = tapply(MAP_damage,sample, max))
df17$cat <-"Maximum"

df18 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 125) %>%
  summarise('125' = tapply(MAP_damage,sample, max))
df18$cat <-"Maximum"

df19 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 150) %>%
  summarise('150' = tapply(MAP_damage,sample, max))
df19$cat <-"Maximum"

df20 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 300) %>%
  summarise('300' = tapply(MAP_damage,sample, max))
df20$cat <-"Maximum"

df21 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 350) %>%
  summarise('350' = tapply(MAP_damage,sample, max))
df21$cat <-"Maximum"

df22 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 375) %>%
  summarise('375' = tapply(MAP_damage,sample, max))
df22$cat <-"Maximum"

df23 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 500) %>%
  summarise('500' = tapply(MAP_damage,sample, max))
df23$cat <-"Maximum"

df24 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 1000) %>%
  summarise('1000' = tapply(MAP_damage,sample, max))
df24$cat <-"Maximum"

df25 <-dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 10) %>%
  summarise('10' = tapply(MAP_damage,sample, median))
df25$cat <-"Median"

df26 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 50) %>%
  summarise('50' = tapply(MAP_damage,sample, median))
df26$cat <-"Median"

df27 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 75) %>%
  summarise('75' = tapply(MAP_damage,sample, median))
df27$cat <-"Median"

df28 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 100) %>%
  summarise('100' = tapply(MAP_damage,sample, median))
df28$cat <-"Median"

df29 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 125) %>%
  summarise('125' = tapply(MAP_damage,sample, median))
df29$cat <-"Median"

df30 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 150) %>%
  summarise('150' = tapply(MAP_damage,sample, median))
df30$cat <-"Median"

df31 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 300) %>%
  summarise('300' = tapply(MAP_damage,sample, median))
df31$cat <-"Median"

df32 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 350) %>%
  summarise('350' = tapply(MAP_damage,sample, median))
df32$cat <-"Median"

df33 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 375) %>%
  summarise('375' = tapply(MAP_damage,sample, median))
df33$cat <-"Median"

df34 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 500) %>%
  summarise('500' = tapply(MAP_damage,sample, median))
df34$cat <-"Median"

df35 <- dt1 %>%
  group_by(sample) %>%
  filter(N_reads >= 1000) %>%
  summarise('1000' = tapply(MAP_damage,sample, median))
df35$cat <-"Median"

X1 <- melt(df4)
X2 <- melt(df4a)
X3 <- melt(df5)
X4 <- melt(df6)
X5 <- melt(df7)
X6 <- melt(df8)
X7 <- melt(df9)
X8 <- melt(df10)
X9 <- melt(df11)
X10 <- melt(df12)
X11 <- melt(df13)
X12 <- melt(df14)
X13 <- melt(df15)
X14 <- melt(df16)
X15 <- melt(df17)
X16 <- melt(df18)
X17 <- melt(df19)
X18 <- melt(df20)
X19 <- melt(df21)
X20 <- melt(df22)
X21 <- melt(df23)
X22 <- melt(df24)
X23 <- melt(df25)
X24 <- melt(df26)
X25 <- melt(df27)
X26 <- melt(df28)
X27 <- melt(df29)
X28 <- melt(df30)
X29 <- melt(df31)
X30 <- melt(df32)
X31 <- melt(df33)
X32 <- melt(df34)
X33 <- melt(df35)
mydata <- rbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, X16, 
                X17, X18, X19, X20, X21, X22, X23, X24, X25, X26, X27, X28, X29, X30, X31, X32, X33)  

#renaming columns
colnames(mydata) <- c("sample", "category", "N_reads", "DNAdamage")

#adding age
mydata$new <-metadataAge$Median[match(gsub("^[^-]+-[^-]+-([^_]+)_.*$", "\\1", mydata$sample), metadataAge$Sample)]
names(mydata)[names(mydata) == 'new'] <- 'age'

pdf(file="figures/MTG.DNAdamageMultiReads.pdf")
ggplot(data=mydata) +
  geom_line(aes(x=as.numeric(age), y=DNAdamage, color=N_reads)) + scale_color_brewer(palette="Spectral") +
  facet_wrap(~category) +
  xlab("Years BP") +
  ylab("DNA damage") +
  labs(color = "Minimum number \n of reads")
dev.off()

## Highlight taxa with more than 500 reads 
mfd <- mydata %>% filter(N_reads == 500) 

## Plot plant taxa and add the min, max and median, size = N_reads. Save as pdf
pdf(file = "figures/MTG.DNAdamageJitterPlot_N_reads.pdf", width = 8, height = 4)
ggplot() +
  geom_jitter(data = dt1, aes(x=as.numeric(Median_age), y=MAP_damage, size = N_reads), alpha =0.5) +
  gghighlight(N_reads > 500) +
  geom_line(data = mfd, aes(x=age, y=DNAdamage, color=category), size = 1) + 
  scale_color_brewer(palette="Spectral") +
  geom_point(data = mfd, aes(x=age, y=DNAdamage, color=category), size = 2) +
  xlab("Years BP") +
  ylab("DNA damage") +
  labs(color = "Values for Plants with \n>500 reads", size = "Number of reads")
dev.off()

## Plot plant taxa and add the min, max and median, size = MAP_significance Save as pdf
pdf(file = "figures/MTG.DNAdamage_JitterPlot_significance.pdf", width = 8, height = 4)
ggplot(data = dt1) +
  geom_jitter(aes(x=as.numeric(Median_age), y=MAP_damage, size = MAP_significance), alpha =0.5) +
  gghighlight(N_reads > 500) +
  geom_line(data = mfd, aes(x=age, y=DNAdamage, color=category), size = 1) + 
  scale_color_brewer(palette="Spectral") +
  geom_point(data = mfd, aes(x=age, y=DNAdamage, color=category), size = 2) +
  xlab("Years BP")+
  ylab("DNA damage") +
  labs(color = "Values for Plants with \n>500 reads", size = "Significance \nfor Taxa with >500 reads")
dev.off()

########Filtering animal taxa (Meazoa) >100 reads based on the values of the plant damage model built above.
#Map_significance
MS2 = 2
# Minimum reads for parsing taxa
MinRead1 = 100
# Minimum mean read length
MinLength1 = 35

# subsetting the table

tb1 <- dt %>% filter(N_reads >= 100, mean_L > MinLength, MAP_significance  > MS2,  grepl("",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", sample))
dt_plant <- dt %>% filter(N_reads >= 500, mean_L > MinLength, MAP_significance  > MS2,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", sample),tax_name!="Zostera")

unique(dt_plant$tax_name) 

##Density plot for plants > 500 reads ####to be checked #######
dtX <- dt1 %>% filter(N_reads >= 500)  

pdf("figures/MTG.PLant_DamageDensityPlot.pdf", height = 6,width = 4) 
dt_plant %>%
  mutate(Median_age = fct_relevel(Median_age,
                                  "8440.5" , "7413", "6499.5",  "5637",  "4800", "2639.5",  "2016.5" , "1586", "997", "465"))%>%
  ggplot(aes(x = MAP_damage, y = Median_age, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01) +
  theme_bw() +
  theme(legend.position = "none") + xlab("MAP_damage") + ylab("Age in Years BP") 
dev.off()

# create filtering table 
## calculate stats and combine data 
pb4 <-dt_plant %>%
  group_by(sample) %>%
  summarise(mode = tapply(MAP_damage,sample, getmode))

pb4a <- dt_plant %>%
  group_by(sample) %>%
  summarise(A_mode = tapply(MAP_damage,sample, estimate_mode))

pb5 <- dt_plant %>%
  group_by(sample) %>%
  summarise(min = tapply(MAP_damage,sample, min))

pb6 <- dt_plant %>%
  group_by(sample) %>%
  summarise(max = tapply(MAP_damage,sample, max))

pb7 <- dt_plant %>%
  group_by(sample) %>%
  summarise(median = tapply(MAP_damage,sample, median))

pb8 <- dt_plant %>%
  group_by(sample) %>%
  summarise(sd = tapply(MAP_damage,sample, sd))

# binding them to table, and adding ages 
pb9 <- cbind(pb4, pb4a[,2], pb5[,2], pb6[,2], pb7[,2], pb8[,2])

#add age
pb9$Age <-metadataAge$Median[match(gsub("^[^-]+-[^-]+-([^_]+)_.*$", "\\1",pb9$sample), metadataAge$Sample)]


#checking the classes of the rows and changing these to be numeric
sapply(pb9, class)
x <- ncol(pb9)
pb9[, 2:x] <- lapply(pb9[, 2:x], as.numeric)
pb9.w <- melt(pb9, id.vars = c("sample", "Age"),measure.vars = c("max","median","min"))

# Join the two data frames based on the sample ID
joined_data <- left_join(tb1, pb9, by = "sample")
joined_data$Age <-metadataAge$Median[match(gsub("^[^-]+-[^-]+-([^_]+)_.*$", "\\1",joined_data$sample), metadataAge$Sample)]


# Filter the joined data frame based on the minimum and maximum values for all genuses
filtered_data_genus <- joined_data %>%
  filter(MAP_damage >= min & MAP_damage <= max)
filtered_data_genus

#add columns
filtered_data_genus <- filtered_data_genus %>%
  mutate(Col_fact = "Fitting")
filtered_data_genus <- filtered_data_genus %>%
  mutate(Filtering = "Ancient")

#dataset discarded taxa min threshold
filtered_data_genus_min <- joined_data %>%
  filter(MAP_damage < min)
filtered_data_genus_min

#add columns then
filtered_data_genus_min <- filtered_data_genus_min %>%
  mutate(Col_fact = "Outlier_min")
filtered_data_genus_min <- filtered_data_genus_min %>%
  mutate(Filtering = "Not_ancient")

#dataset discarded taxa max threshold
filtered_data_genus_max <- joined_data %>%
  filter(MAP_damage > max)
filtered_data_genus_max

#add columns
filtered_data_genus_max <- filtered_data_genus_max %>%
  mutate(Col_fact = "Outlier_max")
filtered_data_genus_max <- filtered_data_genus_max %>%
  mutate(Filtering = "Ancient")

#making one table all genuses
filtered_data_genus_all <- rbind(filtered_data_genus, filtered_data_genus_max, filtered_data_genus_min)
filtered_data_genus_all <- filtered_data_genus_all[filtered_data_genus_all$tax_name!="Loxodonta",]

# Filter the joined data frame based on the minimum and maximum values for only Metazoa
filtered_data_metazoan <- filtered_data_genus_all %>% filter(grepl("Metazoa",tax_path), grepl("\\bgenus\\b", tax_rank))
filtered_data_metazoan$Age <- metadataAge$Median[match(gsub("^[^-]+-[^-]+-([^_]+)_.*$", "\\1", filtered_data_metazoan$sample), metadataAge$Sample)]


###save csv
write.csv(filtered_data_metazoan, file = "cleanedData/filtered_metazoan_all.csv")
write.csv(filtered_data_genus_all, file = "cleanedData/filtered_genus_all.csv")

#####ORIGINAL####
#plot metazoan reads and identify outliers based on the damage model
pdf(file = "figures/OLDMTG.DNAdamageJitterPlot_finalAnimals100readsWithOutliers.pdf", width = 9, height = 6)
ggplot() +
  geom_jitter(data = subset(joined_data, grepl("Metazoa", tax_path)), aes(x=Age, y=MAP_damage, size = N_reads), alpha =0.5) +
  geom_line(data = mfd, aes(x=age, y=DNAdamage, color=category), size = 1) + 
  scale_color_brewer(palette="Spectral") +
  geom_point(data = mfd, aes(x=age, y=DNAdamage, color=category), size = 2) +
  geom_text(data = subset(joined_data, grepl("Metazoa", tax_path) & MAP_damage < min ), aes(x =Age, y = MAP_damage, label = tax_name), size = 3, color = "black") +
  geom_text(data = subset(joined_data, grepl("Metazoa", tax_path) & MAP_damage > max ), aes(x =Age, y = MAP_damage, label = tax_name), size = 3, color = "black") +
  scale_size_continuous(trans="sqrt") +
  xlab("Age in Years BP")
dev.off()

#####LUKE's version ####
#plot metazoan reads and identify outliers based on the damage model


pdf(file = "figures/MTG.DNAdamageJitterPlot_finalAnimals100readsWithOutliers.pdf", width = 7, height = 5)
ggplot() +
  geom_jitter(data = subset(joined_data, grepl("Metazoa", tax_path)), aes(x=Age, y=MAP_damage, size = N_reads), alpha =0.3) +
  geom_line(data = pb9.w, aes(x=Age, y=value, color=variable), size = 1) + 
  scale_color_brewer(palette="Spectral") +
  geom_point(data = pb9.w, aes(x=Age, y=value, color=variable), size = 1) +
  geom_text(data = subset(joined_data, grepl("Metazoa", tax_path) & MAP_damage < min ), aes(x =Age, y = MAP_damage, label = tax_name), size = 3, color = "black") +
  geom_text(data = subset(joined_data, grepl("Metazoa", tax_path) & MAP_damage > max ), aes(x =Age, y = MAP_damage, label = tax_name), size = 3, color = "black") +
  #scale_size_continuous(trans="sqrt") +
  xlab("Cal Years BP")
dev.off()

