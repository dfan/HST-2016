#James Diao
#6/27/2016

load("~/Documents/Kohane_Lab/week_3/ExAC_HCM.RData")

exac_1000g_means <- sapply(data.total.id, function(x) apply(x[,2:3],2,mean))
exac_1000g_cor <- sapply(data.total.id, function(x) cor(x[,2], x[,3]))
colnames(exac_1000g_means) <- HCM.panel
names(exac_1000g_cor) <- HCM.panel

#Correlations
bp <- barplot(exac_1000g_cor[!is.na(exac_1000g_cor)], las = 2, ylim = c(0,1.1), ylab = "Correlation")
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
bp <- barplot(counts, ylim = c(0,1000), las = 2,
              ylab = "Unique Variants", col = c("Black","Red","Blue"), beside = T)
legend("topright",c("1000G","Exac","Both"), pch = 19, col = c("Black","Red","Blue"))
counts
apply(counts,1,median)

### Part 2
penetrance <- seq(0,1,0.05)
prevalence <- 0.002
allelic.het <- 0.001
af_threshold <- prevalence*allelic.het/penetrance
af_threshold
plot(penetrance, af_threshold)

#Output lowest allele freq in each gene
min_af <- sapply(data.exac, function(x) {
  nums <- as.numeric(x[,4])
  min(nums[nums!=0])
  })
barplot(min_af, main = "Lowest Allele Frequency in Each Gene", las = 2)
abline(h=af_threshold)

#Output how many variants in exac pass the cutoff?
keep_prop <- sapply(data.exac, function(x) sum(x[,4] < af_threshold)) / sapply(data.exac, nrow)
barplot(keep_prop, ylim = c(0,0.5), las = 2)
keep_prop

#Contains IDs of the alleles in which genes passed the cutoff
allele_pass <- lapply(data.exac.id, function(x) x[(as.numeric(x[,2]) < af_threshold),1])

#Logical vector for which genes passed the cutoff
genes.pass <- lapply(1:num.genes, function(i) 
  apply(data.1000g.full[[i]][,c(3,5,6)],1, function(y) paste(y,collapse = "-")) %in% allele_pass[[i]]
)
mean.genes.pass <- sapply(genes.pass, mean)
names(mean.genes.pass) <- HCM.panel
barplot(mean.genes.pass, las = 2, main = "Mean Proportion Variants Retained As Pathogenic")

#Logical vector for which patients are diseased
patient.status <- lapply(1:num.genes, function(i) {
  output <- data.1000g.full[[i]][genes.pass[[i]],] #Subset only variants in allele_pass
  colSums(rbind(FALSE, apply(output[,11:ncol(output)], 2, function(y) grepl("1",y)))) > 0
})
diseased <- sapply(patient.status, mean)
names(diseased) <- HCM.panel
barplot(diseased, las = 2, main = "Proportion Patients with HCM")

#Convert genotype data to logical
logical.1000g <- lapply(data.1000g.full, function(input) {
  apply(input[,11:ncol(input)], 2, function(y) grepl("1",y))
})

# True positives P(CT | AT) = P(AT & CT) / P(AT) : 
# Each row of full is a variant, each column is a patient. 
# Collect a vector of patient status <-- this is important
# Take each variant vector AND patient status vector to find number of true positives
# Take sum(patient status vector) to find all actual positives
tpr <- lapply(1:num.genes, function(i) {
  apply(logical.1000g[[i]], 1, function(y) 
      sum(y & patient.status[[i]])
    ) / sum(patient.status[[i]])
})
sapply(tpr, mean) #NaN means everyone is healthy (disease = 0)

# False positives P(CT | AF) = P(CT & AF) / P(AF) : 
# Each row of full is a variant, each column is a patient. 
# Take each variant vector AND !(patient status vector) to find number of false positives
# Take sum(!patient status vector) to find all actual falsehoods
fpr <- lapply(1:num.genes, function(i) {
  apply(logical.1000g[[i]], 1, function(y) 
    sum(y & !patient.status[[i]])
  ) / sum(!patient.status[[i]])
})
sapply(fpr, mean) #NaN means everyone is diseased (healthy = 0)

#Next step: convert to 1D version

#Taking only 1000G variants (rows) that pass the cutoff
#data.1000g.filt <- lapply(1:num.genes, function(i) {
#  output <- data.1000g.full[[i]][genes.pass[[i]],] #Subset only variants in allele_pass
#  if (nrow(output)==0) return(NULL)
#  else rbind(FALSE,apply(output[,11:ncol(output)], 2, function(y) grepl("1",y)))
#    #Extra "FALSE" row added for matrix properties
#})

#diseased <- lapply(data.1000g.filt, function(x) {
#  if (length(x)==0) return(NULL)
#  as.matrix(x[,colSums(x) > 0]) #Subset only patients with variant
#})

#Average variants/patient
#sapply(diseased, function(x) mean(colMeans(x)))

#names(disease_frac) <- HCM.panel
#disease_frac <- sapply(diseased, ncol) / sapply(data.1000g.filt, ncol)
#barplot(disease_frac, las = 2)

#variant_frac <- sapply(data.1000g.filt, nrow) / sapply(data.1000g, nrow)
#variant_frac 
#barplot(disease_frac, las=2)




