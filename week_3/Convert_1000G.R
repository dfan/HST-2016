#James Diao
#6/27/2016

library(scrapeR)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/")
LMM.page <- scrape(url="http://personalizedmedicine.partners.org/Laboratory-For-Molecular-Medicine/Tests/Cardiomyopathy/HCM-Panel.aspx", headers=FALSE, parse=TRUE)
HCM.panel <- levels(readHTMLTable(LMM.page[[1]])[[1]]$Gene)
num.genes <- length(HCM.panel)
header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

import.file.exac <- function(gene) {
  name <- paste(gene,"exac.csv", sep = "_")
  output <- read.csv(paste(getwd(),name,sep="/ExAC_HCM/"), stringsAsFactors = FALSE)
  output[nchar(paste(output$Alternate,output$Reference))==3,]
}

import.file.1000g <- function(gene, cut) {
  name <- paste(gene,"genotypes.vcf", sep = "_")
  output <- read.table(paste(getwd(),name,sep="/1000G_HCM/"), stringsAsFactors = FALSE)
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

data.1000g.full <- lapply(HCM.panel, function(x) import.file.1000g(x, FALSE))
data.exac.full <- lapply(HCM.panel, import.file.exac) #0.5 seconds, 3 Mb
names(data.1000g.full) <- HCM.panel
names(data.exac.full) <- HCM.panel

#Take relevant subset (position, ref, alt, allele.freq)
data.1000g <- lapply(data.1000g.full, function(x) {
  columns <- sapply(c("POS","REF","ALT","AF"), function(y) which(colnames(x) %in% y))
  x[order(x$POS),columns]
})
data.exac <- lapply(data.exac.full, function(x) {
  columns <- sapply(c("Position","Reference","Alternate","Allele.Frequency"), function(y) which(colnames(x) %in% y))
  x[order(x$Position),columns]
})
#Combine position-ref-alt for left-join
data.1000g.id <- lapply(data.1000g, function(x) cbind(apply(x[,-4],1, function(y) paste(y,collapse = "-")),x[,4]) )
data.exac.id <- lapply(data.exac, function(x) cbind(apply(x[,-4],1, function(y) paste(y,collapse = "-")),x[,4]) )

#Left-join
data.total.id <- lapply(1:num.genes, function(i) {
  shared <- intersect(data.1000g.id[[i]][,1], data.exac.id[[i]][,1])
  output <- data.frame(shared,
                       as.numeric(data.1000g.id[[i]][ data.1000g.id[[i]][,1] %in% shared, 2]),
                       as.numeric(data.exac.id[[i]][ data.exac.id[[i]][,1] %in% shared, 2]))
  colnames(output) <- c("ID","1000G","Exac")
  output
})
names(data.total.id) <- HCM.panel

#gold_std thresholds

gold_std <- c(penetrance = 0.8, prevalence = 0.002, allelic.het = 0.0033)
af_threshold <- as.numeric(gold_std['prevalence'] * gold_std['allelic.het'] / gold_std['penetrance'])

#Logical vector for which variants passed the cutoff
var.logical.1 <- lapply(1:num.genes, function(i) {
  data.1000g.id[[i]][,1] %in%
    data.exac.id[[i]][(as.numeric(data.exac.id[[i]][,2]) < af_threshold),1]
})

var.logical.2 <- lapply(1:num.genes, function(i) {
  max_pop <- apply((select(data.exac.full[[i]], contains("Allele.Count."))) /
                     (select(data.exac.full[[i]], contains("Allele.Number."))),1,max)
  data.1000g.id[[i]][,1] %in%
    data.exac.id[[i]][(max_pop < af_threshold),1]
})

var.logical <- var.logical.1
var.pass <- data.frame(ID = unlist(lapply(1:num.genes, function(i) data.1000g.id[[i]][,1])),
                       Pass = unlist(var.logical))
glimpse(var.pass)
mean(var.pass$Pass)

empty.row <- rep(FALSE,ncol(data.1000g.full[[1]])-10)
patient.status <- sapply(1:num.genes, function(i) {
  sub <- data.1000g.full[[i]][var.logical[[i]],11:ncol(data.1000g.full[[i]])] #Subset only variants in allele_pass
  if (nrow(sub)==0) empty.row
  sub %>% colSums>0
})
std_ps <- rowSums(patient.status)>0
names(std_ps) <- NULL
mean(std_ps)


mean(var.pass$Pass)
std_var <- var.pass
rm(var.pass)

pop_max_af <- unlist(lapply(data.exac.full, function(dt) {
  apply((select(dt, contains("Allele.Count."))) /
          (select(dt, contains("Allele.Number."))),1,max)
}))




#save.image("~/Documents/Kohane_Lab/week_4/ExAC_1000G_HCM.RData")







