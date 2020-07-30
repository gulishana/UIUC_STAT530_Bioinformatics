# install package edgeR
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR") #run two times to install

## Problem 2

## Part(a)
setwd("C:/Users/Gulishana/Desktop/HW6_dm_counts")

library(edgeR);
targets=readTargets(file="dm_targets.txt");
raw_counts=readDGE(targets$file,comment.char="_",header=F);
dim(raw_counts$counts)

## filter out genes with low expression
keep=rowSums(cpm(raw_counts)>1)>=3;
counts=raw_counts[keep,keep.lib.sizes=F];

## recalculate the library size
names(counts);
dim(counts$counts);

## calculate TMM normalization factors
counts=calcNormFactors(counts);

## fit a two-way ANOVA
design=model.matrix(~sex+mating+sex*mating,data=targets);
print(design)

## use tagwise dispersion
counts=estimateGLMCommonDisp(counts,design);
counts=estimateGLMTrendedDisp(counts,design);
counts=estimateGLMTagwiseDisp(counts,design);

## fit the GLM regression
fit=glmFit(counts,design,dispersion=counts$tagwise.dispersion);
head(fit$coefficients);
genes=c("FBgn0000018","FBgn0000042","FBgn0025712")
fit$coefficients[genes,]


## Part(b)
## make contrast to test Ho: beta0+beta3=0
colnames(design)=c("intercept","male","polygamy","male_polygamy")
contrast=makeContrasts(male+male_polygamy,levels=make.names(colnames(design)));
test=glmLRT(fit,contrast = contrast);
de_poly=decideTestsDGE(test,p.value=0.1)
table(de_poly)


## Part(c)
## One way: use the male- and female-biased genes in polygamous flies from Part(b)
male_biased=which(de_poly==1)
female_biased=which(de_poly==-1)
unbiased=which(de_poly==0)

## make contrast to test Ho: -(beta2+beta3)=0
contrast=makeContrasts(-(polygamy+male_polygamy),levels=make.names(colnames(design)));

test_male_biased=glmLRT(fit[male_biased,],contrast = contrast);
male=decideTestsDGE(test_male_biased,p.value=0.1)
table(male);

test_female_biased=glmLRT(fit[female_biased,],contrast = contrast);
female=decideTestsDGE(test_female_biased,p.value=0.1)
table(female)

test_unbiased=glmLRT(fit[unbiased,],contrast = contrast);
un=decideTestsDGE(test_unbiased,p.value=0.1)
table(un)

## test the difference significance between samples
t.test(male,un,alternative = "less")
t.test(female,un,alternative = "greater")
t.test(un,mu=0,alternative = "two.sided")


## Another way:
## First to identify the male- and female-biased genes in monogamous flies by FDR 0.1
bias_test=glmLRT(fit,coef=2);

de_mono=decideTestsDGE(bias_test,p.value=0.1); 
table(de_mono);
male_biased=which(de_mono==1)
female_biased=which(de_mono==-1)
unbiased=which(de_mono==0)

## make contrast to test Ho: -(beta2+beta3)=0
contrast=makeContrasts(-(polygamy+male_polygamy),levels=make.names(colnames(design)));

test_male_biased=glmLRT(fit[male_biased,],contrast = contrast);
male=decideTestsDGE(test_male_biased,p.value=0.1)
table(male);

test_female_biased=glmLRT(fit[female_biased,],contrast = contrast);
female=decideTestsDGE(test_female_biased,p.value=0.1)
table(female)

test_unbiased=glmLRT(fit[unbiased,],contrast = contrast);
un=decideTestsDGE(test_unbiased,p.value=0.1)
table(un)

## test the difference significance between samples
t.test(male,un,alternative = "less")
t.test(female,un,alternative = "greater")



## Problem 3
n=1000
sims=200
mu_rep=rep(0,length=n)
mu_seq=seq(0,5,length=n)
MLE_vs_JS=matrix(NA,nrow=sims,ncol=4)
colnames(MLE_vs_JS)=c("MLE_rep","JS_rep","MLE_seq","JS_seq")

for(i in 1:sims){
  X_rep=rnorm(n,mean=mu_rep,sd=1);
  X_seq=rnorm(n,mean=mu_seq,sd=1);
  MLE_rep=X_rep
  MLE_seq=X_seq
  JS_rep =X_rep*( 1 - (n-2)/(t(X_rep)%*%X_rep) );
  JS_seq =X_seq*( 1 - (n-2)/(t(X_seq)%*%X_seq) );
  
  MLE_vs_JS[i,1]=sum( (MLE_rep-mu_rep)^2 )/n;
  MLE_vs_JS[i,2]=sum( ( JS_rep-mu_rep)^2 )/n;
  MLE_vs_JS[i,3]=sum( (MLE_seq-mu_seq)^2 )/n;
  MLE_vs_JS[i,4]=sum( ( JS_seq-mu_seq)^2 )/n;
}

colSums(MLE_vs_JS)/sims