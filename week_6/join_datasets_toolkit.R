#James Diao
#July 23, 2016
#Kohane Lab

library(dplyr)
library(tidyr)
library(ggbiplot)

setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
HST.folder.6 <- function(file.name) { paste(getwd(), file.name, sep = "/HST-2016/week_6/") }
merged <- read.csv(file = HST.folder.6("merged_output.csv"), stringsAsFactors = F)
### Sanity checks
#mean(merged$Start==merged$Stop) #Confirmed start = stop
#mean(merged$Assertion=="pathogenic")
#sort(levels(merged$Chromosome)) #All present
merged$Stop <- NULL
merged$Assertion <- NULL
merged <- rename(merged, Position = Start)
str(merged)
length(unique(merged$Disease)) == 2858 #what...
to.rm <- (merged$Ref=="-") | (merged$Alt=="-")
mean(to.rm) # == 0.05531253
merged <- merged[-to.rm,]
merged <- unite(merged, col = "ID", Chromosome, Position, Ref, Alt, sep = "-", remove = FALSE)

merged$ID %>% length() # == 10286
merged$ID %>% unique %>% length() # == 8441



