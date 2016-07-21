#James Diao
#July 11, 2016
#Kohane Lab

library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
load(file = paste(getwd(), "ACMG.files", sep = "/HST-2016/week_5/"))
populations <- c("ASW","CEU","YRI")

map <- read.table(file = paste(getwd(), "phase3map.txt",sep = "/1000G/"), header = T)
submap <- as.data.frame(sapply(populations, function(pop) {
  map$pop == pop
}))
submap <- cbind(submap, "Total" = TRUE)
row.names(submap) <- map$sample

population.data <- sapply(ACMG.1000g, function(data) {
  sapply(1:ncol(submap), function(i) {
    (subset(data, select = c(rep(FALSE,10), submap[,i])) %>% colSums > 0) %>% mean
  })
}) %>% t %>% tbl_df
colnames(population.data) <- c(populations,"Total")

write.table(population.data, sep = "\t", quote = FALSE,
  paste(getwd(), "Population_Variant_Status.txt",sep="/HST-2016/week_5/"))







