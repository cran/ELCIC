## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!require("devtools")) {
#    install.packages("devtools")
#  }
#  devtools::install_github("chencxxy28/ELCIC")

## -----------------------------------------------------------------------------
library(ELCIC)
library(MASS)
library(mvtnorm)

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  set.seed(28582)
#  inter.total=20
#  count.matrix<-matrix(0,nrow=4,ncol=7)
#  candidate.sets<-list(c(1,2),c(1,3),c(1,4),c(1,2,3),c(1,2,4),c(1,3,4),c(1,2,3,4))
#  rownames(count.matrix)<-c("ELCIC","AIC","BIC","GIC")
#  colnames(count.matrix)<-candidate.sets
#  samplesize<-100
#  ov=2
#  
#  for(iter in 1:inter.total)
#  {
#      #generate data
#  simulated.data<-glm.generator(beta=c(0.5,0.5,0.5,0),samplesize=samplesize,rho=0.5,dist="poisson")
#  y<-simulated.data[["y"]]
#  x<-simulated.data[["x"]]
#  
#  criterion.all<-ELCIC.glm(x=x,y=y,candidate.sets=candidate.sets,name.var=NULL,dist="poisson")
#  
#  #print(iter)
#  index.used<-cbind(1:4,apply(criterion.all,1,which.min))
#  count.matrix[index.used]<-count.matrix[index.used]+1
#  }
#  count.matrix/inter.total

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  set.seed(28582)
#  inter.total=20
#  count.matrix<-matrix(0,nrow=4,ncol=7)
#  candidate.sets<-list(c(1,2),c(1,3),c(1,4),c(1,2,3),c(1,2,4),c(1,3,4),c(1,2,3,4))
#  rownames(count.matrix)<-c("ELCIC","AIC","BIC","GIC")
#  colnames(count.matrix)<-candidate.sets
#  samplesize<-100
#  ov=2
#  
#  for(iter in 1:inter.total)
#  {
#      #generate data
#  simulated.data<-glm.generator(beta=c(0.5,0.5,0.5,0),samplesize=samplesize,rho=0.5,dist="NB",ov=2)
#  y<-simulated.data[["y"]]
#  x<-simulated.data[["x"]]
#  
#  criterion.all<-ELCIC.glm(x=x,y=y,candidate.sets=candidate.sets,name.var=NULL,dist="poisson")
#  
#  #print(iter)
#  index.used<-cbind(1:4,apply(criterion.all,1,which.min))
#  count.matrix[index.used]<-count.matrix[index.used]+1
#  }
#  count.matrix/inter.total

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  library(PoisNor)
#  library(bindata)
#  library(geepack)

## ----message=FALSE,eval=FALSE-------------------------------------------------
#  set.seed(28582)
#  inter.total=20
#  samplesize<-300
#  time=3
#  dist="poisson"
#  candidate.sets<-list(c(1,2),c(1,3),c(1,4),c(1,2,3),c(1,2,4),c(1,3,4),c(1,2,3,4))
#  candidate.cor.sets<-c("independence","exchangeable","ar1")
#  
#  count.matrix.elcic<-count.matrix.qic<-matrix(0,nrow=length(candidate.cor.sets),ncol=length(candidate.sets))
#  rownames(count.matrix.elcic)<-rownames(count.matrix.qic)<-candidate.cor.sets
#  colnames(count.matrix.elcic)<-colnames(count.matrix.qic)<-candidate.sets
#  
#  for(iter in 1:inter.total)
#  {
#      #generate data
#  data.corpos<-gee.generator(beta=c(-1,1,0.5,0),samplesize=samplesize,time=time,num.time.dep=2,num.time.indep=1,rho=0.4,x.rho=0.2,dist="poisson",cor.str="exchangeable",x.cor.str="exchangeable")
#  
#  y<-data.corpos$y
#  x<-data.corpos$x
#  id<-data.corpos$id
#  r<-rep(1,samplesize*time)
#  
#  criterion.elcic<-ELCIC.gee(x=x,y=y,r=r,id=id,time=time,candidate.sets=candidate.sets,name.var=NULL,dist="poisson",candidate.cor.sets=candidate.cor.sets)
#  
#  criterion.qic<-QICc.gee(x=x,y=y,id=id,dist=dist,candidate.sets=candidate.sets, name.var=NULL, candidate.cor.sets=candidate.cor.sets)
#  
#  index.used.elcic<-which(criterion.elcic==min(criterion.elcic),arr.ind = TRUE)
#  count.matrix.elcic[index.used.elcic]<-count.matrix.elcic[index.used.elcic]+1
#  
#  index.used.qic.row<-which(candidate.cor.sets==rownames(criterion.qic))
#  index.used.qic.col<-which.min(criterion.qic)
#  count.matrix.qic[index.used.qic.row,index.used.qic.col]<-count.matrix.qic[index.used.qic.row,index.used.qic.col]+1
#  
#  print(iter)
#  }
#  count.matrix.elcic/inter.total
#  count.matrix.qic/inter.total

