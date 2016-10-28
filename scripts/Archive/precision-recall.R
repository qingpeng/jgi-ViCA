setwd("~/Documents/rna_virus/poster2")
library(ROCR)
library(viridis)
#data with 2608 RdRp containing contigs from metatranscriptomes and 10,000 viral contigs
data<-read.csv("all.txt", header=F,sep=",")

hist(qlogis(data[data$V2==0,1]), col=rgb(0.392,.549,0.564,0.7), xlim=c(-20,10), breaks=50,ylim=c(0,1250), main="Logit of RdRP viral contigs and bacterial contigs", ylab="Non-viral Frequency",xlab="Logit")
par(new=T)
hist(qlogis(data[(data$V2==1),1]), col=rgb(0.827,0.492,0.156,0.7),xlim=c(-20,10), breaks=50,ylim=c(0,250), main="", xlab="", ylab="",yaxt="n")
axis(4)

rdrp.pd <- prediction(data$V1,data$V2)
rdrp.pe <- performance(rdrp.pd, measure="prec", x.measure="rec")
plot(rdrp.pe, colorize=T,colorize.palette=viridis(256,1), lwd=4, main="Classification of RdRP contigs")


metat <- read.csv("gamma_meta_predictionAndLabels.txt", header=F,sep=",")
par(mar=c(5.1,4.1,4.1,2.1))
hist(metat[metat$V1<0.5,1], col=rgb(0.392,.549,0.564,1), xlim=c(0,1), breaks=50,ylim=c(0,150000), main="Scores of metatranscriptome contigs", ylab="Frequency",xlab="Score")
par(new=T)
hist(metat[metat$V1>0.5,1], col=rgb(0.827,0.492,0.156,1),xlim=c(0,1), breaks=50,ylim=c(0,2000), main="", xlab="", ylab="",yaxt="n")
axis(4)
abline(v=0.5, lty=2)


head(perf)
logisticreg<-read.csv("pr2_data.csv")
svm.pd <- prediction(svm$X2,svm$X1)
logisticreg.pd <- prediction(logisticreg$X2,logisticreg$X1)
svm.pe <- performance(svm.pd, measure="prec", x.measure="rec")
logisticreg.pe <- performance(logisticreg.pd, measure="prec", x.measure="rec")
plot(svm.pe, downsampling=0.1, colorize=T, lwd=2, main="SVM classifier")
plot(logisticreg.pe, downsampling=0.1, colorize=T, lwt=2, main="Logistic classifier")

#Visualize the split in samples
pdf('SVMeval.pdf',width=8.5,height=11)
par(mfrow=c(2,1))
plot(svm.pe, downsampling=0.1, colorize=T, colorize.palette=viridis(256,1), lwd=4, main="SVM classifier")
hist(svm[svm$X1=='nonviral',2], col=rgb(1,0,0,0.7), xlim=c(-30,30), breaks=120,ylim=c(0,100000), main="Scores of viral and nonviral classifications", ylab="Nonviral Frequency",xlab="Score")
par(new=T)
hist(svm[(svm$X1=='viral'),2], col=rgb(0,1,1,0.7),xlim=c(-30,30), breaks=20,ylim=c(0,1000), main="", xlab="", ylab="",yaxt="n")
axis(4)
dev.off()

#Visualize the split in samples
pdf('LReval.pdf',width=8.5,height=11)
par(mfrow=c(2,1))
plot(logisticreg.pe, downsampling=0.1, colorize=T,colorize.palette=viridis(256,1), lwd=4, main="Logistic classifier")
hist(logisticreg[logisticreg$X1=='nonviral',2], col=rgb(1,0,0,0.7), xlim=c(-30,30), breaks=120,ylim=c(0,100000), main="Scores of viral and nonviral classifications", ylab="Nonviral Frequency",xlab="Score")
par(new=T)
hist(logisticreg[(logisticreg$X1=='viral'),2], col=rgb(0,1,1,0.7),xlim=c(-30,30), breaks=20,ylim=c(0,1000), main="", xlab="", ylab="",yaxt="n")
axis(4)
dev.off()



hist(data[data$V1=='non-virus',2], col=rgb(1,0,0,0.7), xlim=c(-30,30), breaks=120,ylim=c(0,100000), main="Scores of viral and nonviral classifications", ylab="Nonviral Frequency",xlab="Score")
par(new=T)
hist(data[(data$V1=='virus'),2], col=rgb(0,1,1,0.7),xlim=c(-30,30), breaks=60,ylim=c(0,10000), main="", xlab="", ylab="",yaxt="n")
axis(4)
