#James Diao
#6/27/2016

load("~/Documents/Kohane_Lab/week_3/ExAC_HCM.RData")
library(dplyr)
library(tidyr)
exac_1000g_means <- sapply(data.total.id, function(x) apply(x[,2:3],2,mean))
exac_1000g_cor <- sapply(data.total.id, function(x) cor(x[,2], x[,3]))
colnames(exac_1000g_means) <- HCM.panel
names(exac_1000g_cor) <- HCM.panel

#1000G and ExAC
plot(sapply(data.total.id, function(dti) {
  perc.diff <- (dti[,2]-dti[,3])/dti[,2]
  mean(perc.diff[!is.infinite(perc.diff)])
  }), type = 'h', ylab = "(1000G-ExAC)/(1000G)")


#Correlations
bp <- barplot(exac_1000g_cor[!is.na(exac_1000g_cor)], main = "Correlations", ylim = c(0,1.1), ylab = "Correlation (r)")
exac_1000g_cor

#Means
bp <- barplot(exac_1000g_means, main = "Means ",ylab = "Means", col = c("Black","Red"), las = 2, beside = T)
legend("topright",c("1000G","Exac"), pch = 19, col = c("Black","Red"))
exac_1000g_means

#Variant Counts
counts <- rbind(sapply(data.1000g.id, nrow), 
                sapply(data.exac.id, nrow),
                sapply(data.total.id, nrow))
rownames(counts) <- c("1000G", "ExAC", "Both")
bp <- barplot(counts, ylim = c(0,1000), main = "Variant Counts         ", las = 2,
              ylab = "Unique Variants", col = c("Black","Red","Blue"), beside = T)
legend("topright",c("1000G","Exac","Both"), pch = 19, col = c("Black","Red","Blue"))
counts
apply(counts,1,median)

### Part 2

penetrance <- seq(0,1,0.05)
prevalence <- 0.002
allelic.het <- 0.001
af_threshold <- prevalence*allelic.het/penetrance
names(af_threshold) <- penetrance
af_threshold
plot(af_threshold, ylab = "Allele Frequency Threshold", xlab = "Penetrance Cutoff")

gene <- "MYBPC3"
var <- 1
gene.exac.id <- data.exac.id[[gene]]
gene.1000g.full <- data.1000g.full[[gene]]
gene.logical.1000g <- unlist(data.1000g.full[[gene]][var,11:ncol(data.1000g.full[[gene]])])
mean(gene.logical.1000g)

#Logical vector for which genes passed the cutoff
genes.pass <- lapply(af_threshold, function(af) {
  ap <- gene.exac.id[(as.numeric(gene.exac.id[,2]) < af),1]
  apply(gene.1000g.full[,c(3,5,6)],1, function(y) paste(y,collapse = "-")) %in% ap
})
mean.genes.pass <- sapply(genes.pass, mean)
mean.genes.pass
plot(names(mean.genes.pass),mean.genes.pass, xlab = "Penetrance Cutoff", ylab = "Proportion 1000G Genes Passed")

#Logical vector for which patients are diseased
patient.status <- lapply(genes.pass, function(gp) 
  gene.1000g.full[gp,11:ncol(gene.1000g.full)] %>% rbind(FALSE) %>% colSums>0 #Subset only variants in allele_pass
)

#Proportion of HCM-positive
sapply(patient.status, mean)
plot(penetrance, .Last.value, ylab = "Proportion HCM-Positive", ylim = c(0,1))

# True positives P(CT | AT) = P(AT & CT) / P(AT) : 
# Each row of *input* is a variant, each column is a patient. 
# Take each variant vector AND patient status vector to find number of true positives
# Take sum(patient status vector) to find all actual positives

# False positives P(CT | AF) = P(CT & AF) / P(AF) : 
# Each row of *input* is a variant, each column is a patient. 
# Take each variant vector AND !(patient status vector) to find number of false positives
# Take sum(!patient status vector) to find all actual falsehoods

roc <- sapply(patient.status, function(ps) 
  rbind(sum(gene.logical.1000g & ps) / sum(ps), #tpr, #NaN means everyone is healthy (disease = 0)
        sum(gene.logical.1000g & !ps) / sum(!ps)) #fpr, #NaN means everyone is diseased (healthy = 0)
)
rownames(roc) <- c("True Positive Rate", "False Positive Rate")
roc <- roc[, !is.na(colSums(roc))]
plot(names(roc[1,]), roc[1,], xlab = "Penetrance", ylab = "True Positive Rate", type = 'l')
plot(names(roc[2,]), roc[2,], xlab = "Penetrance", ylab = "False Positive Rate", type = 'l')
plot(roc[2,], roc[1,], xlab = "False Positive Rate", ylab = "True Positive Rate")
roc


#Contains IDs of the alleles in which genes passed the cutoff
#allele_pass <- lapply(af_threshold, function(af) gene.exac.id[(as.numeric(gene.exac.id[,2]) < af),1])
#Proportion of passes from exac
#sapply(allele_pass, length) / nrow(gene.exac.id)
#Logical vector for which genes passed the cutoff
#genes.pass <- lapply(allele_pass, function(ap) apply(gene.1000g.full[,c(3,5,6)],1, function(y) paste(y,collapse = "-")) %in% ap)
#proportion of passes from 1000G
#sapply(genes.pass, mean)



#FIRST IDENTIFY TRUE SICK
#gs_penetrance <- 0.8 #Add 0
#gs_prevalence <- 0.002
#gs_allelic.het <- 0.004
#gs_af_threshold <- prevalence*allelic.het/penetrance
#gs_af_threshold
#Contains IDs of the alleles in which genes passed the cutoff
#allele_pass <- gene.exac.id[(as.numeric(gene.exac.id[,2]) < gs_af_threshold),1]
#Proportion of passes from exac
#length(allele_pass)/ nrow(gene.exac.id)
#Logical vector for which genes passed the cutoff
#genes.pass <- apply(gene.1000g.full[,c(3,5,6)],1, function(y) paste(y,collapse = "-")) %in% allele_pass
#proportion of passes from 1000G
#mean(genes.pass)
#Logical vector for which patients are diseased
#output <- gene.1000g.full[genes.pass,] #Subset only variants in allele_pass
#patient.status <- colSums(rbind(FALSE, apply(output[,11:ncol(output)], 2, function(y) grepl("1",y)))) > 0
#Proportion of HCM-positive
#mean(patient.status)
