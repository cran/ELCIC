data("geesimdata")

x<-geesimdata$x

test_that("output error: y is not a vector",{expect_error(ELCIC.gee.single(x=geesimdata$x,y=as.data.frame(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="poisson",corstr="exchangeable"),"y should be in a vector format")})

test_that("output error: x is not a matrix",{expect_error(ELCIC.gee.single(x=as.data.frame(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="poisson",corstr="exchangeable"),"x should be in a matrix format")})

test_that("output error: dist is not undefined",{expect_error(ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="gamma",corstr="exchangeable"),"Invalid type of dist. It should be one of gaussian,binomial,poisson")})


test_that("output error: non-unique index",{expect_error(ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=c(1,1,3),name.var=NULL,dist="poisson",corstr="exchangeable"),"Invalid candidate model provided")})

test_that("output error: undefined candidate model",{expect_error(ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2","x3","x4"),dist="poisson",corstr="exchangeable"),"Invalid candidate model provided")})

test_that("output error: non-unique variable name",{expect_error(ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x1"),dist="poisson",corstr="exchangeable"),"Invalid candidate model provided")})

test_that("output error: invalid correlation structure for outcomes in gee without missing",{expect_error(ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar2"),"Invalid type of correlation structure for outcomes. It should be one of ar1,exchangeable,independence")})

output1<-ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar1",joint=TRUE)
output2<-ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar1",joint=TRUE)
test_that("output equal: equal values given both index and name, when joint=true",{expect_equal(output1,output2)})

output1<-ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar1",joint=FALSE)
output2<-ELCIC.gee.single(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar1",joint=FALSE)
test_that("output equal: equal values given both index and name, when joint=false",{expect_equal(output1,output2)})

output1<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=TRUE)
output2<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=TRUE)
test_that("output equal: equal values given both index and name, when joint=true",{expect_equal(output1,output2)})


output1<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),dist="poisson",candidate.cor.sets="ar1",joint=TRUE)
output2<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),dist="poisson",candidate.cor.sets="ar1",joint=TRUE)
test_that("output equal: equal values given both index and name, when joint=true",{expect_equal(output1,output2)})


output1<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=FALSE)
output2<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=FALSE)
test_that("output equal: equal values given both index and name, when joint=false",{expect_equal(output1,output2)})


output1<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),dist="poisson",candidate.cor.sets="ar1",joint=FALSE)
output2<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),dist="poisson",candidate.cor.sets="ar1",joint=FALSE)
test_that("output equal: equal values given both index and name, when joint=false",{expect_equal(output1,output2)})


output1<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=FALSE)
output2<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets=c("ar1","exchangeable"),joint=FALSE)
test_that("output equal: equal values given more than 1 correlation structures, when joint=false",{expect_equal(output1,output2)})


output1<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="independence",joint=FALSE)
output2<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets=c("independence","exchangeable"),joint=FALSE)
test_that("output equal: equal values given more than 1 correlation structures, when joint=false",{expect_equal(output1,output2)})

r<-rep(1,nrow(x))
time_index<-rep(seq_len(3),times=nrow(x)/3)
r[time_index==3]<-rbinom(nrow(x)/3,1,0.5)
x.missing<-geesimdata$x
y.missing<-geesimdata$y
x.missing[r==0,-1]<-NA
y.missing[r==0]<-NA
output1<-ELCIC.gee(x=x.missing,y=(geesimdata$y),r=r,id=geesimdata$id,time=3,candidate.sets=list(c(1,2)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="exchangeable",joint=TRUE)
output2<-ELCIC.gee(x=x.missing,y=y.missing,r=r,id=geesimdata$id,time=3,candidate.sets=list(c(1,2)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets=c("exchangeable"),joint=TRUE)
test_that("output equal: equal values given more than 1 correlation structures, when joint=TRUE",{expect_equal(output1,output2)})

output1<-ELCIC.gee(x=x.missing,y=(geesimdata$y),r=r,id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="exchangeable",joint=FALSE)
output2<-ELCIC.gee(x=x.missing,y=y.missing,r=r,id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets=c("exchangeable"),joint=FALSE)
test_that("output equal: equal values given more than 1 correlation structures, when joint=FALSE",{expect_equal(output1,output2)})


#generate other dist
set.seed(515413)
samplesize<-200
time<-3
geesimdata<-gee.generator(beta=c(-1,1,0.5,0),samplesize=samplesize,time=time,num.time.dep=2,num.time.indep=1,rho=0.4,x.rho=0.2,dist="gaussian",cor.str="exchangeable",x.cor.str="exchangeable")
rownames(geesimdata$x)<-seq_len(nrow(geesimdata$x))
colnames(geesimdata$x)<-c("intercept","x1","x2","x3")

output1<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(geesimdata$x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="gaussian",candidate.cor.sets="ar1",joint=FALSE)
output2<-ELCIC.gee(x=(geesimdata$x),y=(geesimdata$y),r=rep(1,nrow(geesimdata$x)),id=geesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="gaussian",candidate.cor.sets=c("ar1","exchangeable"),joint=FALSE)
test_that("output equal: equal values given more than 1 correlation structures, when joint=false",{expect_equal(output1,output2)})

