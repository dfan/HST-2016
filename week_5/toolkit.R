#James Diao
#July 12, 2016
#Kohane Lab

library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
load(file = paste(getwd(), "ACMG.files", sep = "/HST-2016/week_5/"))

sapply(ACMG.1000g, nrow) -> varnum
plot(varnum,xlab = "Number of Variants", main = NULL)
plot(gene.tx.length, xlab = "Gene Length (tx region)", main = NULL)
cor(varnum,gene.tx.length)

sapply(1:num.genes, function(i) {
  output <- ACMG.1000g[[i]][,-c(1:10)] %>% colSums
}) %>% rowSums -> var.sites
hist(var.sites, xlab = "(0|1)  (1|1)  (1|0)  Sites per Patient", main = NULL, breaks = 30)

sapply(1:num.genes, function(i) {
  output <- ACMG.1000g[[i]][,-c(1:10)] %>% colSums
  nrow(ACMG.1000g[[i]]) - output
}) %>% rowSums -> null.sites
hist(null.sites, xlab = "0|0 Sites per Patient", main = NULL, breaks = 30)





#Size
#gene.tx.length
#Variant Number
#