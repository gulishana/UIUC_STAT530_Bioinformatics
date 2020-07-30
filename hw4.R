## Problem 2
setwd("D:/???Courses???/STAT530/Homework/HW4_eQTL")# for homework
setwd("C:/Users/Gulishana/Desktop/HW4_eQTL") # for Notebook PDF

## Part 1: markers report
library(MatrixEQTL);
useModel = modelLINEAR; 
SNP_file_name = "hw4_markers.txt";
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Part 2: transcripts report of each brain region (10 regions in total)
region <- c("CRBL","FCTX","HIPP","MEDU","OCTX","PUTM","SNIG","TCTX","THAL","WHMT");
for(r in region){
  
  expression_file_name = paste("expr_",r,".txt",sep="");
  output_file_name = paste("res_",r,".txt",sep="");
  pvOutputThreshold = 1e-5;
  errorCovariance = numeric();
  gene = SlicedData$new();
  gene$fileDelimiter = " ";       # the SPACE character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1+291704;   # one row of column labels and all exons
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  cvrt = SlicedData$new();        # no covariates
  ## Run the analysis
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = TRUE);
  
}

## Part 3: find the most significant eQTLs in each brain region (10 regions in total)
res <- read.table("res_CRBL.txt",header=TRUE);
print(res[which.min(res$p.value),]);
## Then repeat the same code for the other 9 regions

## OR perform a for loop for all 10 regions together
region <- c("CRBL","FCTX","HIPP","MEDU","OCTX","PUTM","SNIG","TCTX","THAL","WHMT");
for(r in region){
  res <- read.table(paste("res_",r,".txt",sep=""), header=TRUE);
  print(cbind(r,res[which.min(res$p.value),]));
}

  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Problem 3
## Part 1: perform an unpooled FDR adjustment of p-values at 1% FDR
n <- 26493*34373;
alpha <- 0.01;
region <- c("CRBL","FCTX","HIPP","MEDU","OCTX","PUTM","SNIG","TCTX","THAL","WHMT");
for(r in region){
  res <- read.table(paste("res_",r,".txt",sep=""), header=TRUE);
  cat(r,":",sum(p.adjust(res$p.value,method="fdr",n=n)<=alpha),"\n");
}


## Part 2: pool p-values first and perform a global FDR adjustment of p-values at 1% FDR
region <- c("CRBL","FCTX","HIPP","MEDU","OCTX","PUTM","SNIG","TCTX","THAL","WHMT");
p <- rep(list(NA),length(region));
for(i in 1:length(region)){
  res <- read.table(paste("res_",region[i],".txt",sep=""), header=TRUE);
  p[[i]] <- res$p.value;
}
p <- unlist(p);
alpha <- 0.01;
sum(p.adjust(p,method="fdr",n=26493*34373*10)<=alpha);



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Problem 4
## Identify the significant eQTLs in each brain region
n <- 26493*34373;
alpha <- 0.01;
pairs <- matrix(NA,nrow=0,ncol=2);
region <- c("CRBL","FCTX","HIPP","MEDU","OCTX","PUTM","SNIG","TCTX","THAL","WHMT");
for(r in region){
  res <- read.table(paste("res_",r,".txt",sep=""), header=TRUE);
  keep <- p.adjust(res$p.value,method="fdr",n=n)<=alpha;
  pairs <- rbind(pairs,res[keep,1:2]);
}

pairs <- pairs[!duplicated(pairs),];
save(pairs,file = "pairs.Rdata");
nrow(pairs);



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Problem 5
## Hierarchical clustering and plotting heatmap and dendrograms
## Part 1: calculate the Z-scores and save the data
load("pairs.Rdata")

library(data.table);
markers <- fread("hw4_markers.txt");
setkey(markers,"id");
zscore <- matrix(NA,nrow=684,ncol=0);
region <- c("CRBL","FCTX","HIPP","MEDU","OCTX","PUTM","SNIG","TCTX","THAL","WHMT");
for(r in region){
  expr <- fread(paste("expr_",r,".txt",sep=""));
  setkey(expr,"ExprID");
  zscore <- cbind(zscore,apply(pairs,1,function(x){
    g <- as.numeric(expr[x[2],-1,with=F]);
    s <- as.numeric(markers[x[1],-1,with=F]);
    p <- summary(lm(g~s))$coef[2,4];
    return(-qnorm(p/2));
  }));
}

colnames(zscore) <- region;
save(zscore,file = "zscore.Rdata");


## heatmapping
load("zscore.Rdata")
hc.eqtl <- hclust(as.dist((1-cor(t(zscore)))/2),method="complete");

# 50 colors between blue and red
colors <- bluered(50)
library(gplots);
heatmap.2(zscore,Rowv=as.dendrogram(hc.eqtl),
          labRow=FALSE,
          col=colors,key=TRUE,density.info="none",trace="none",cexCol=1);



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Problem 6
## DAVID analysis

load("pairs.Rdata")
load("zscore.Rdata")
hc.eqtl <- hclust(as.dist((1-cor(t(zscore)))/2),method="complete");
clust <- cutree(hc.eqtl,k=3);
table(clust);

for(i in 1:3){
  trans <- sapply(as.character(pairs[clust==i,"gene"]),function(x){
    substr(x,2,nchar(x));
    });
  write.table(unique(trans),file=paste("cluster_",i,".txt",sep=""),
            row.names=FALSE,col.names=FALSE,quote=FALSE);
}

library(data.table);
cluster_1 <- fread("cluster_1.txt"); nrow(cluster_1);
cluster_2 <- fread("cluster_2.txt"); nrow(cluster_2);
cluster_3 <- fread("cluster_3.txt"); nrow(cluster_3);

