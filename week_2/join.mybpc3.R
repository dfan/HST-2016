#James Diao
#6/23/2016

setwd("/Users/jamesdiao/Documents/Kohane_Lab/Week_2/")
mybpc3 <- read.table(paste(getwd(),"MYBPC3_genotypes.txt",sep = "/"))
clinvar.mybpc3 <- read.csv(paste(getwd(),"clinvar_result.csv",sep = "/"), stringsAsFactors = F)
location <- sort(clinvar.mybpc3$Location)
#collect ONLY missense
location <- location[location!=""]
location <- unique(sapply(strsplit(location,' - '), function(x) x[1]))
#clinvar.mybpc3 <- clinvar.mybpc3[clinvar.mybpc3$Location %in% location, ]

output <- sapply(sapply(as.character(mybpc3$V8), function(x) strsplit(x,";")), function(x) x[2])
output <- paste(gsub("AF=", "sum(", output),")",sep = "")
output <- sapply(output, function(x) eval(parse(text=x)))
names(output) <- mybpc3$V2
plot(output)
output

#ntd.mybpc3
genotype <- mybpc3[,10:ncol(mybpc3)]
colnames(genotype) <- 1:ncol(genotype)
rownames(genotype) <- mybpc3$V2
set.1 <- apply(genotype,2,function(x)
  sapply(x, function(y) 
    as.numeric(strsplit(as.character(y),"[|]")[[1]][1])
  ))
set.2 <- apply(genotype,2,function(x)
  sapply(x, function(y) 
    as.numeric(strsplit(as.character(y),"[|]")[[1]][2])
  ))
set <- cbind(set.1, set.2)

ntd_pair <- function(pair, nts) {
  sum(xor(pair[1,], pair[2,]))/nts
}

ntd.simul <- function(cohort, pairs, nts) {
  sapply(1:pairs, function(x) ntd_pair(
    cohort[sample(nrow(cohort),2),], nts
  ))
}

hist(colSums(set))
ntd.mybpc3 <- ntd.simul(set, 10000, max(mybpc3$V2)-min(mybpc3$V2)+1)
plot(ntd.mybpc3)
mean(ntd.mybpc3)




