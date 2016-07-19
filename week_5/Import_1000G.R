#James Diao
#July 11, 2016
#Kohane Lab

library(dplyr)
library(tidyr)
library(scrapeR)
library(RMySQL)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
load(file = paste(getwd(), "ACMG.files", sep = "/HST-2016/week_5/"))

#########################################################
###  Function for importing 1000 genomes / ExAC data  ###
#########################################################

import.file.exac <- function(gene) {
  setwd("/Users/jamesdiao/Documents/Kohane_Lab/ExAC")
  name <- paste(gene,"exac.csv", sep = "_")
  output <- read.csv(paste(getwd(),name,sep="/"), stringsAsFactors = FALSE)
  output[nchar(paste(output$Alternate,output$Reference))==3,]
}

import.file.1000g <- function(gene, cut) {
  setwd("/Users/jamesdiao/Documents/Kohane_Lab/1000G")
  name <- paste(gene,"genotypes.vcf", sep = "_")
  output <- read.table(paste(getwd(),name,sep="/"), stringsAsFactors = FALSE)
  #Add header
  if(cut) output <- output[,1:length(header)]
  names(output)[1:length(header)] <- header
  #Remove all single alt indels
  output <- output[nchar(output$REF)==1,] #deletions
  output <- output[grepl(",",output$ALT) | nchar(output$ALT)==1,] #insertions
  #Add AF Column
  af <- sapply(output$INFO, function(x) strsplit(x,";")[[1]][2] )
  af <- sapply(af, function(x) substring(x,4,nchar(x)))
  names(af) <- NULL
  output <- cbind("AF"=af, output)
  output$AF <- as.character(output$AF)
  #Separate paired and unpaired
  alt <- strsplit(output$ALT,",")
  paired <- sapply(alt,length)!=1
  if (sum(paired)!=0) {
    sub.out <- output[paired,]
    output <- output[!paired,]
    #View all rows with multiple alts
    #sub.out[,1:6]
    #Remove all paired indels
    master <- lapply(lapply(alt[paired],nchar), function(x) x==1)
    new.af <- lapply(sub.out$AF, function(x) strsplit(x,",")[[1]])
    new.alt <- lapply(sub.out$ALT, function(x) strsplit(x,",")[[1]])
    sub.out$AF <- lapply(1:nrow(sub.out), function(x) new.af[[x]][master[[x]]])
    sub.out$ALT <- lapply(1:nrow(sub.out), function(x) new.alt[[x]][master[[x]]])
    #Add back the rows that are now single missense
    add.back <- sub.out[sapply(sub.out$AF, length)==1,]
    add.back$AF <- unlist(add.back$AF)
    add.back$ALT <- unlist(add.back$ALT)
    output <- rbind(output,add.back)
    #Extract (or actually, just sum together) double missenses
    sub.out <- sub.out[sapply(sub.out$AF, length) >1,]
    if (nrow(sub.out)!=0) {
      #how many cycles
      max.alts <- max(sapply(sub.out$ALT,length))
      for (i in 1:max.alts) {
        #Each temp is used holds the ith REF-ALT match
        temp <- sub.out
        temp$ALT <- sapply(1:nrow(sub.out), function(x) sub.out$ALT[[x]][i])
        temp$AF <- sapply(1:nrow(sub.out), function(x) sub.out$AF[[x]][i])
        output <- rbind(output,temp[!is.na(temp$AF),])
      }
    }
  }
  if (!cut)
    output[,11:ncol(output)] <- apply(output[,11:ncol(output)], 2, function(y) grepl("1",y))
  output
}


#########################
###  Import all data  ###
#########################

# Import ExAC data for all ACMG
ACMG.exac <- lapply(gene.list, import.file.exac)
#Take relevant subset (position, ref, alt, allele.freq)
ACMG.exac.id <- lapply(ACMG.exac, function(x) {
  columns <- select(x, Chrom, Position, Reference, Alternate, Allele.Frequency) %>%
    arrange(Chrom, Position, Reference, Alternate) %>% unite("ID", c(Chrom,Position, Reference, Alternate), sep= "-")
})

# Import 1000G data for all ACMG
ACMG.1000g <- lapply(gene.list, import.file.1000g, cut = F)
#Take relevant subset (position, ref, alt, allele.freq)
ACMG.1000g.id <- lapply(ACMG.1000g, function(x) {
  columns <- select(x, CHROM, POS, REF, ALT, AF) %>%
    arrange(CHROM, POS, REF, ALT) %>% unite("ID", c(CHROM, POS, REF, ALT), sep= "-")
})
names(ACMG.exac) = names(ACMG.exac.id) = names(ACMG.1000g) = names(ACMG.1000g.id) = gene.list

setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
#save.image(file = paste(getwd(), "ACMG.files", sep = "/HST-2016/week_5/"))









