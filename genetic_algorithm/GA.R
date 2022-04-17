#utf-8 encoding
setwd("D:\\胶质瘤\\上传的代码\\3Sscore\\genetic_algorithm")###please change this path to the "genetic_algorithm" file
library(MASS)
library(genalg)
library(timeROC)
library(survivalsvm)
library(survminer)
library(survival)
library(survivalROC)
library("gbm")
library(randomForestSRC)
library(ipflasso)
library(party)
data=read.table(".\\693\\693riskfile.csv",sep = ",",row.names = 1,header = T)
###
gatest=read.table(".\\TCGA\\TCGAriskfile.csv",sep = ",",row.names = 1,header = T)
gaevaluate <- function(indices) {
  result = 1
  if (sum(indices) > 0) {
    newdata=cbind(data[,1:2],data[,(which(indices[]==1)+2)])
    newgatest=cbind(gatest[,1:2],gatest[,(which(indices[]==1)+2)])
    colnames(newgatest)=colnames(newdata)
    
    rsfmodel <- rfsrc(Surv(futime, fustat) ~ ., newgatest, ntree = 50, 
                      nodedepth = 10,importance = TRUE)
    
    pred = predict(object = rsfmodel, 
                   newdata = newgatest)$predicted
    predicresult=Hmisc::rcorr.cens(-pred, Surv(newgatest$futime, newgatest$fustat))
    
    result=1-predicresult[1]
    
  }
  result
}
monitor <- function(obj) {
  minEval = min(obj$evaluations);
  plot(obj, type="hist");
}
woppa <- rbga.bin(size=47, mutationChance=0.05, zeroToOneRatio=10,
                  evalFunc=gaevaluate, verbose=TRUE, monitorFunc=monitor,popSize=50, iters=50)
plot(woppa)
