######## Edit distance plots #########
###packages
library(tidyverse) 
library(reshape2)
library(ggplot2)
library(gridExtra)
###########
#make count tables
count_table1 <- table(Eleph2650$ED)
print(count_table1)
data_matrix1 <- as.matrix(count_table1)

count_table2 <- table(Mamm2650$ED)
print(count_table2)
data_matrix2 <- as.matrix(count_table2)

#make proportion tables
b1 <- prop.table(data.matrix(data_matrix1), margin=2)
colSums(prop.table(b1, margin=2)) # should give 100 for each column

b2 <- prop.table(data.matrix(data_matrix2), margin=2) 
colSums(prop.table(b2, margin=2)) # should give 100 for each column

#transposing to long format table for ggplot
long_format1 <- melt(b1)
long_format1 <- long_format1[, -2]
colnames(long_format1) <- c("ED", "Read_prop")

long_format2 <- melt(b2)
long_format2 <- long_format2[, -2]
colnames(long_format2) <- c("ED", "Read_prop")

long_format1_Ele <- cbind(long_format1, Taxa = "Elephant")
long_format2_Mamm <- cbind(long_format2, Taxa = "Mammoth")

#join data Mammoth-Elephant reads
BindTable_E_M <- bind_rows(long_format1_Ele, long_format2_Mamm)

#make combined plots
pdf(file = "Read_prop_ED_all.pdf", width = 4, height = 6)
ggplot(BindTable_E_M, aes(x=ED , y=Read_prop)) +
  geom_path(aes(col=Taxa))+
  geom_point(aes(col=Taxa))+
  scale_color_manual(values = c("red", "blue"))+
  labs(x = "Mismatches", y = "Read proportions") +
  theme_minimal() + theme(aspect.ratio=10/2)
dev.off()

#or make single plots
p1 <- ggplot(long_format1, aes(x = ED, y = Read_prop)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +  # Add points (optional)
  labs(title = "Elephant reads", x = "Mismatches", y = "Read proportions") +
  theme_minimal()
p1

p2 <- ggplot(long_format2, aes(x = ED, y = Read_prop)) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "red", size = 3) +  # Add points (optional)
  labs(title = "Mammoth reads", x = "Mismatches", y = "Read proportions") +
  theme_minimal()
p2

grid.arrange(p1, p2,
              ncol = 2, nrow = 1)
