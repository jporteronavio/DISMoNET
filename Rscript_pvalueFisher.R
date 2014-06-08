# R script to calculate matrix p-value (positive dataset, negative dataset, positive rest of genome, negative rest of  genome)

#check multtest is installed

library(multtest)

#setwd("/path/R")
#source("/path/pvalueFisher.R")
#data<-read.table("/path/tmp", header=TRUE, sep="\t")
#a<-matrix(c(data$DatasetPositive[i], data$DatasetNegative[i], data$GenomePositive[i], data$GenomeNegative[i]), 2)


ftest<-NULL
data<-read.table("./data/R/contingencytable.tmp", sep="\t")
for(i in 1:nrow(data)){
a<-matrix(c(data$V1[i], data$V2[i], data$V3[i], data$V4[i]), 2)

tmp1<-fisher.test(a, alternative = "greater")$p.value
ftest<-c(ftest,tmp1)
rm(tmp1)}

procs <- c("Bonferroni", "BH")

if (nrow(data)>1) {
	pvaladjF <- mt.rawp2adjp(ftest, procs)

	pvaladjFOK<-cbind(pvaladjF$index, pvaladjF$adjp)

	write.table(pvaladjFOK, "./data/R/resultsR")
} else {
	
	write.table(ftest, "./data/R/resultsR")

}
