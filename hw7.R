## Problem 1

## Preprocess the multiple myeloma data
setwd("D:/???Courses???/STAT530/Homework/HW7")
setwd("C:/Users/Gulishana/Desktop/HW")

rm(list=ls());
library(data.table);

pData <- fread("GSE24080_series_matrix.txt",
               skip = 29, nrows = 37, header = FALSE);
expr <- fread("GSE24080_series_matrix.txt",
              skip = 67, nrow = 54675, header = TRUE, fill = TRUE);

## training data: train = 1
train <- sapply(pData[10,-1,with=F],function(x){
  as.numeric(strsplit(x," ")[[1]][2]=="Training");
});

## remove "fake" subjects
train[1:4] = 0;

## testing data: test = 1 (categories: Training, Validation, MAQC_Remove)
test <- sapply(pData[10,-1,with=F],function(x){
  as.numeric(strsplit(x," ")[[1]][2]=="Validation");
});

## outcome
Y <- sapply(pData[15,-1,with=F],function(x){
  as.numeric(strsplit(x,": ")[[1]][3]);
});

## age
age <- sapply(pData[11,-1,with=F],function(x){
  as.numeric(strsplit(x,": ")[[1]][2]);
});

## sex
sex <- sapply(pData[12,-1,with=F],function(x){
  as.numeric(strsplit(x,": ")[[1]][2]=="female");
});

## training data
X.train <- cbind(as.matrix(t(expr[,-1,with=F]))[train==1,],
                 age[train==1],sex[train==1]);
colnames(X.train) <- c(unlist(expr[,1,with=F]),"age","sex");
Y.train <- Y[train==1];

## testing data
X.test<- cbind(as.matrix(t(expr[,-1,with=F]))[test==1,],
               age[test==1],sex[test==1]);
colnames(X.test) <- c(unlist(expr[,1,with=F]),"age","sex");
Y.test <- Y[test==1];




## Part (a)
library(glmnet)
nrow(X.train)
nrow(X.test)


## Training: find the lambda that gives the smallest 5-fold Cross-Validation misclassification error
CV5=cv.glmnet(X.train,Y.train,family="binomial",alpha=1,type.measure="class",nfolds=5) 
names(CV5)
lambda=CV5$lambda.min;  lambda
plot(CV5)


## Training: build a logistic regression model on the training data
fit_lasso = glmnet(X.train,Y.train,family="binomial",alpha=1,lambda=lambda)
summary(fit_lasso)
head(fit_lasso$beta)
which(fit_lasso$beta!=0)
length(which(fit_lasso$beta!=0))


## Test: use cutoff c=0.5 to report misclassification rate of the test data
pred.Y.test=predict.glmnet(fit_lasso,X.test,type="response")

misclass_rate= mean( ifelse( Y.test!=(pred.Y.test>=0.5),1,0 ) ); misclass_rate





## Part (b)

## Calculate the ratio of "0" and "1" in the true response of the test data, respectively
length(which(Y.test==1))/length(Y.test)

length(which(Y.test==0))/length(Y.test)





## Part (c)

## Preprocess data
ncol(X.train)
ncol(X.test)
X.train=X.train[,1:2000]
X.test=X.test[,1:2000]
data_training=data.frame(Y=Y.train,X=X.train)
data_test=data.frame(Y=Y.test,X=X.test)


## Training: build a random forest model on the training data
library(ranger)
fit_rf=ranger(Y~.,data=data_training,classification=TRUE,write.forest = TRUE)
summary(fit_rf)


## Test: use cutoff c=0.5 to report misclassification rate of the test data
pred.Y.test=predict(fit_rf,data=data_test,type="response")
pred.Y.test=pred.Y.test$predictions

misclass_rate= mean( ifelse( Y.test!=(pred.Y.test>=0.5),1,0 ) ); misclass_rate
