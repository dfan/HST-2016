#James Diao
#July 12, 2016
#Kohane Lab

library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
load(file = paste(getwd(), "ACMG.files", sep = "/HST-2016/week_5/"))

sapply(1:num.genes, function(i) {
  output <- ACMG.1000g[[i]][,-c(1:10)] %>% colSums
}) %>% rowSums %>% hist(xlab = "Variants per Patient")



