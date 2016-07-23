# James Diao 
# June 20, 2016

library(RMySQL)
library(shiny)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/Week_1/")
all_cons <- dbListConnections(MySQL())
for (con in all_cons)
  dbDisconnect(con)
con <- dbConnect(MySQL(), user = 'root', 
                 password = 'root', dbname = 'kohane_lab', host = 'localhost',
                 unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
query <- "select Position, Allele_Count,Allele_Number,Allele_Frequency,
                Allele_Count_African, Allele_Number_African,
                Allele_Count_East_Asian, Allele_Number_East_Asian,
                Allele_Count_European_Non_Finnish, Allele_Number_European_Non_Finnish,
                Allele_Count_Finnish, Allele_Number_Finnish,
                Allele_Count_Latino, Allele_Number_Latino,
                Allele_Count_Other, Allele_Number_Other,
                Allele_Count_South_Asian, Allele_Number_South_Asian from"
input.mybpc3 <- dbGetQuery(con, paste(query,"MYBPC3"))
input.myh7 <- dbGetQuery(con, paste(query,"MYH7"))

#nts.gene = number of nts in that gene
#Assumption: nts = range in chromosomal positions
nts.mybpc3 <- max(input.mybpc3$Position)-min(input.mybpc3$Position)+1
nts.myh7 <- max(input.myh7$Position)-min(input.myh7$Position)+1

# Takes ancestry (African) and gene dataset input (input.MYBPC3) and outputs "Position, Count, Number, Freq" 
# This also adds all counts from duplicate positions (totals appear to be the same across duplicates).
process <- function(ancestry, gene.input) {
  #Collect and name ancestry-specific columns
  out <- data.frame(gene.input["Position"], 
                    gene.input[paste("Allele_Count_", ancestry, sep = "")], 
                    gene.input[paste("Allele_Number_", ancestry, sep = "")])
  colnames(out) <- c("Position", "Count", "Number")
  #Collect nonunique items
  pos_table <- table(out$Position)
  nonunique <- names(pos_table)[pos_table>1]
  #create add_on = all duplicates that have been combined into a single data frame
  add_on <- as.data.frame(t(sapply(nonunique, function(x) 
            c(sum(out$Count[out$Position==x]), 
              mean(out$Number[out$Position==x])))))
  add_on <- cbind(rownames(add_on), add_on)
  colnames(add_on) <- colnames(out)
  #remove all duplicates
  out <- out[!(out$Position %in% nonunique),]
  #add combined duplicates back in
  out <- rbind(out, add_on)
  #remove all variants that have 0 counts
  out <- out[out$Count!=0,]
  rownames(out) <- NULL
  #calculate and append Freq = count/number
  out$Freq <- out$Count/out$Number
  out
}

#Extract data for MYBPC3 and MYH7
ancestries <- c("African", "East_Asian", "European_Non_Finnish", "Finnish","Latino","Other","South_Asian")
data.mybpc3 <- lapply(ancestries, function(x) process(x, input.mybpc3))
data.myh7 <- lapply(ancestries, function(x) process(x, input.myh7))
names(data.mybpc3) = names(data.myh7) = c("Afr", "E.Asian", "NFE", "Finnish","Latino","Other","S.Asian")
#Plot the number of unique variants
counts <- rbind(sapply(data.mybpc3, nrow),sapply(data.myh7, nrow))
rownames(counts) <- c("MYBPC3","MYH7")
bp <- barplot(counts, ylim = c(0,max(apply(counts,2,sum))*1.2), main = "Unique Variant Counts by Ethnicity",
              ylab = "Number of Unique Variants", col = c("Black","Red"), beside = T)
text(bp, counts, counts, pos = 3)
legend("topright",c("MYBPC3","MYH7"), pch = 19, col = c("Black","Red"))

#Looks like there's on avg 4-5 counts/chromosome
counts_per_chr <- sum(input.mybpc3$Allele_Count)/mean(input.mybpc3$Allele_Number)

#Assume independence
#Multiple all pairs of sequence probabilities, multiply by 1/nts
#Sum these up
ntd_one_per <- function(freq, nts) {
  freq <- freq/sum(freq)
  2*sum(apply(combn(freq,2),2,prod))*2/nts
}

#Calculate nucleotide diversity for mybpc3
ntd.mybpc3 <- sapply(data.mybpc3, function(x) ntd_one_per(x$Freq,nts.mybpc3))
ntd.myh7 <- sapply(data.myh7, function(x) ntd_one_per(x$Freq,nts.myh7))
ntd.1 <- rbind(ntd.mybpc3, ntd.myh7)
bp <- barplot(ntd.1, ylim = c(0,max(ntd.1)*1.3), main = "Unique Variant Counts by Ethnicity",
              ylab = "Number of Unique Variants", col = c("Black","Red"), beside = T)
legend("topright",c("MYBPC3","MYH7"), pch = 19, col = c("Black","Red"))
ntd.1

#Takes around 1.5 seconds
system.time(sapply(data.mybpc3, function(x) ntd_one_per(x$Freq,nts.mybpc3)))

#For testing purposes
freq <- data.mybpc3[[4]]$Freq
nts <- nts.mybpc3
#sum(apply(combn(freq,4),2,prod))*4/nts + sum(apply(combn(freq,3),2,prod))*3/nts + sum(apply(combn(freq,2),2,prod))*2/nts

#m <- 1
#nt_diversity <- -(1-2*m/nts)/nts*log(sum(apply(combn(freq,4)^2,2,prod)))/nts
#nt_diversity

pi.avg <- function(nts, m) {
  sapply(0:m, function(i) (2*m-i)/nts)
}
pi.weights <- function(n, m) {
  sapply(0:m, function(i) 
    choose(n, 2*m-i)*
      choose(2*m-i, i)
  )
}
#from n, pick m, with m possible overlaps
pi.est <- function (nts, n, m) {
  sum(pi.avg(nts,m)*pi.weights(n,m))/sum(pi.weights(n,m))
}

### Calculating nucleotide diversity - Part 2
ntd_m_per <- function(freq, nts, m) {
  freq <- freq/sum(freq)
  2*pi.est(nts,length(freq),m)*sum(apply(combn(freq,2*m),2,prod))
}
ntd.mybpc3 <- sapply(data.mybpc3, function(x) ntd_m_per(x$Freq, nts.mybpc3, 1))
ntd.myh7 <- sapply(data.myh7, function(x) ntd_m_per(x$Freq, nts.myh7, 1))
ntd.m.2 <- rbind(ntd.mybpc3, ntd.myh7)
#Correlation for mybpc3
cor(ntd.1[1,], ntd.m.2[1,])
#Correlation for myh7
cor(ntd.1[2,], ntd.m.2[2,])

# For m = 2 and 77  variants (Fin), we need to sum up choose(77 , 4) values = 1.5 million
# For m = 2 and 340 variants (Afr), we need to sum up choose(340, 4) values = 500 million. 
# For m = 2 and 850 variants (NFE), we need to sum up choose(850, 4) values = 20,000 million.

