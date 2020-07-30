#   Logistic regression without covariates
setwd("D:/VirtualBox_Ubuntu/Ubuntu_LTS_shared folders/STAT530/hapmap1")
nopc <- read.table("nopc.assoc.logistic",header=TRUE)
head(nopc[order(nopc$P),])

# convert pca file to txt file
evec <- read.table("qcd.pca.evec");
fid <- sapply(as.character(evec$V1),function(x){ strsplit(x,":")[[1]][1]; });
iid <- sapply(as.character(evec$V1),function(x){ strsplit(x,":")[[1]][2]; });
out <- data.frame(FID=fid,IID=iid,
                  PC1=evec$V2,
                  PC2=evec$V3,
                  PC3=evec$V4);
write.table(out,file="pcs.txt",row.names=FALSE,quote=FALSE);

# Analysis of PC1
pc1 <- read.table("pc1.assoc.logistic",header=TRUE)
pc1 <- pc1[pc1$TEST=="ADD",]
head(pc1[order(pc1$P),])

# Analysis of PC2
pc2 <- read.table("pc2.assoc.logistic",header=TRUE)
pc2 <- pc2[pc2$TEST=="ADD",]
head(pc2[order(pc2$P),])

# Analysis of PC3
pc3 <- read.table("pc3.assoc.logistic",header=TRUE)
pc3 <- pc3[pc3$TEST=="ADD",]
head(pc3[order(pc3$P),])

#  Plotting PC1
setwd("D:/VirtualBox_Ubuntu/Ubuntu_LTS_shared folders/STAT530/hapmap1")
pc1 <- read.table("pc1.assoc.logistic",header=TRUE)
pc1 <- pc1[pc1$TEST=="ADD",]
library(qqman)
png(file="Manhattan.png",width=500,heigh=300)
manhattan(pc1[,c("SNP","CHR","BP","P")])
dev.off()
