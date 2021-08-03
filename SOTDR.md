Scalar on time-by-distribution regression and its application to
modelling the associations between daily-living physical activity and
cognitive functions in Alzheimer’s Disease
================
Rahul Ghosal
04/21/2021











% Transpose

% Inverse

% Trace

% epigraph

% vector

% vector element

% n-th vector

% vector

% vector

% vector element

% matrix

% matrix element

% matrix

% matrix

% n-th matrix

# Introduction

This document presents an illustration of the SOTDR and SOTDR-L method
proposed in Ghosal et al. (2021) using continuously collected physical
activity (PA) data on AD and CNC participants from KU-ADC. First, we
illustrate modelling cognitive status (CNC vs AD) using subject-specific
average PA (Section 3.1, Model 1 of the paper).

``` r
#####Load the data and preprocess for modelling############
#preprocessing for extracting average PA
### minute-csv folder contains acceloremetry data for each subject #############
path=getwd()
path=paste(path,"/minute-csv",sep="")
setwd(path)
filelist<-list.files(path)
myfiles <-lapply(filelist, read.csv)
subjdf<-matrix(0,nrow = 92,ncol=6)
subjdf<-as.data.frame(subjdf)
#actmat<-matrix(0,nrow = 92,ncol=144)
for(i in 1:92)
{
  mydata<-myfiles[[i]]
  mydata$ADStatus<-ifelse(mydata$ADStatus=="Yes",1,0)
  mydata$Sex<-ifelse(mydata$Sex=="Male",1,0)
  #Date, Time , Vector.Magnitude
  subdata<-mydata[,c("Date","Time","Vector.Magnitude")]
  library(chron)
  subdata$Time<-60 * 24 * as.numeric(times(subdata$Time))
  meanact<-mean(subdata$Vector.Magnitude)
  subjdf[i,]<-c(mydata$id[1],mydata$Age[1],mydata$ADStatus[1],mydata$Sex[1],mydata$YearsOfEducation[1],meanact)
  }
names(subjdf)<-c("id","age","adstatus","Sex","education","meanact")
#summary(subjdf)
#removing NA value#################################################
subjdf[2,]$education<-c(NA)
subjdf$education[2]<-mean(subjdf$education,na.rm=TRUE)
indad<-which(subjdf$adstatus==1)
#####ploting Average PA by cognitive status (CNC vs AD)
library(vioplot)
vioplot(subjdf$meanact[-indad],subjdf$meanact[indad],ylab="Average physical activity",names=c("CNC","AD"),col=c("blue","red"),cex.lab=1.4)
```

![](SOTDR_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
#### Illustrating regression using subject-specific average PA
glm1<-glm(adstatus~age+Sex+education+meanact,data=subjdf,family = "binomial")
summary(glm1)
```

    ## 
    ## Call:
    ## glm(formula = adstatus ~ age + Sex + education + meanact, family = "binomial", 
    ##     data = subjdf)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.6641  -0.8532  -0.3705   0.7287   2.1532  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  7.608403   3.566866   2.133 0.032918 *  
    ## age         -0.051176   0.038077  -1.344 0.178948    
    ## Sex          2.133907   0.554386   3.849 0.000119 ***
    ## education   -0.224314   0.091423  -2.454 0.014144 *  
    ## meanact     -0.005541   0.002050  -2.703 0.006873 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 125.400  on 91  degrees of freedom
    ## Residual deviance:  92.314  on 87  degrees of freedom
    ## AIC: 102.31
    ## 
    ## Number of Fisher Scoring iterations: 5

Higher subject-specific average PA is significantly associated
(alpha=0.05) with a lower odds of AD.

Next, we ilustrate a functional data analysis approach in the paper
using subject-specific temporal curves (Section 3.2, Model 2 of the
paper). We display the average diurnal PA curves (smoothed) for CNC and
AD group.

``` r
#####Load the data and preprocess for temporal modelling############
####Code for extracting diurnal mean activity profile from minute level data####
rm(list=ls())
path=getwd()
path=paste(path,"/minute-csv",sep="")
setwd(path)
filelist<-list.files(path)
myfiles <-lapply(filelist, read.csv)
subjdf<-matrix(0,nrow = 92,ncol=5)
subjdf<-as.data.frame(subjdf)
actmat<-matrix(0,nrow = 92,ncol=144)
for(i in 1:92)
{
  mydata<-myfiles[[i]]
  mydata$ADStatus<-ifelse(mydata$ADStatus=="Yes",1,0)
  mydata$Sex<-ifelse(mydata$Sex=="Male",1,0)
  ##Date, Time , Vector.Magnitude
  subdata<-mydata[,c("Date","Time","Vector.Magnitude")]
  library(chron)
  subdata$Time<-60 * 24 * as.numeric(times(subdata$Time))
  n<-nlevels(subdata$Date)
  lev<-levels(subdata$Date)
  tbtick<-seq(0,1440,by=10)
  actvec<-c()
  l=length(tbtick)-1
  for (j in 1:l)
  {
    tempdataj<-subdata[subdata$Time>=tbtick[j]&subdata$Time <tbtick[j+1],]
    actvec[j]<-mean(tempdataj$Vector.Magnitude)
  }
 subjdf[i,]<-c(mydata$id[1],mydata$Age[1],mydata$ADStatus[1],mydata$Sex[1],mydata$YearsOfEducation[1])
actmat[i,]<-actvec
}
names(subjdf)<-c("id","age","adstatus","Sex","education")
#removing NA value#################################################
subjdf[2,]$education<-c(NA)
subjdf$education[2]<-mean(subjdf$education,na.rm=TRUE)
##################
subjdf$actmat<-actmat
################plot average smoothed diurnal activity profiles#########
tbtick<-seq(0,1440,by=10)
binmid<-tbtick+5
binmid<-binmid[-145]
indAD<-which(subjdf$adstatus==1)
actad<-subjdf$actmat[indAD,]
actcontrol<-subjdf$actmat[-indAD,]
cnavg<-colMeans(actcontrol)
adavg<-colMeans(actad)
l1<-loess.smooth(binmid, cnavg,span = 0.25)
l2<-loess.smooth(binmid, adavg,span=0.25)
plot(l1$x/60,l1$y,type="l",col="blue",xaxt="n",xlab="Time of the Day",ylab="Average activity",cex.lab=1.2)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE)
lines(l2$x/60,l2$y,type="l",col="red",lty=2)  
legend('topleft',c("CNC ","AD") , 
       lty=c(1,2), col=c("blue", "red"), bty='n', cex=.95)
```

![](SOTDR_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#SOFR for modelling association of diurnal average PA
#curve with cognitive status
library(refund)
fit.lf0 <- pfr(adstatus~ age+Sex+education+lf(actmat,argvals=(binmid/60), k=12, bs="ps",m=2),data = subjdf,family='binomial')
summary(fit.lf0)
```

    ## 
    ## Family: binomial 
    ## Link function: logit 
    ## 
    ## Formula:
    ## adstatus ~ age + Sex + education + s(x = actmat.tmat, by = L.actmat, 
    ##     k = 12, bs = "ps", m = 2)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  6.54844    3.61481   1.812 0.070054 .  
    ## age         -0.04002    0.03919  -1.021 0.307151    
    ## Sex          2.11057    0.55287   3.817 0.000135 ***
    ## education   -0.21277    0.09004  -2.363 0.018131 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                         edf Ref.df Chi.sq p-value  
    ## s(actmat.tmat):L.actmat   2      2  7.554  0.0229 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.272   Deviance explained = 26.7%
    ## -REML = 63.213  Scale est. = 1         n = 92

``` r
aa<-coef(fit.lf0,n=144)

###### plot diurnal effect of PA profile on log odds of AD ##########
par(mar=c(5,5,4,2))
plot(binmid/60,aa$value,xaxt="n",xlab="Time of the day",ylab=expression(paste(beta(t))),main="Estimated temporal effect of diurnal PA on log odds of AD",type="l",ylim=c(-.0020,0.0008),cex.lab=1.2)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE,cex = 1.2)
abline(h=0)
upper<-as.numeric(aa$value)+1.96*aa$se
lower<-as.numeric(aa$value)-1.96*aa$se
lines(aa$acargvals.argvals,upper,lty=2)
lines(aa$acargvals.argvals,lower,lty=2)
abline(v=seq(4,20,by=4),lty=3,col="gray20")
```

![](SOTDR_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

Higher PA during morning hours (∼ 10 am-3 p.m) is significantly
associated (*α* = 0.05) with a lower odds of AD. Next, we perform
distributional data analysis using subject-specific quantile functions
of PA (Section 3.3, Model 3 of the paper). We derive the
subject-specific quantile functions and plot the avreage quantile
functions of PA for CNC and AD.

``` r
#####Load the data and preprocess for distributional modelling############
rm(list=ls())
path=getwd()
path=paste(path,"/minute-csv",sep="")
setwd(path)
filelist<-list.files(path)
myfiles <-lapply(filelist, read.csv)
subjdf<-matrix(0,nrow = 92,ncol=5)
subjdf<-as.data.frame(subjdf)
actmat<-matrix(0,nrow = 92,ncol=101)
for(i in 1:92)
{
  mydata<-myfiles[[i]]
  mydata$ADStatus<-ifelse(mydata$ADStatus=="Yes",1,0)
  mydata$Sex<-ifelse(mydata$Sex=="Male",1,0)
  #Date, Time , Vector.Magnitude
  subdata<-mydata[,c("Date","Time","Vector.Magnitude")]
  library(chron)
  subdata$Time<-60 * 24 * as.numeric(times(subdata$Time))
  n<-nlevels(subdata$Date)
  lev<-levels(subdata$Date)
  
  actvec<-quantile(subdata$Vector.Magnitude,probs=c(seq(0,1,l=101)),na.rm=TRUE)
  
  #cos_coeff = ActCosinor(x = actvec, window = 10) #window of 10 minutes
  subjdf[i,]<-c(mydata$id[1],mydata$Age[1],mydata$ADStatus[1],mydata$Sex[1],mydata$YearsOfEducation[1])
  actmat[i,]<-actvec
  }
names(subjdf)<-c("id","age","adstatus","Sex","education")
#summary(subjdf)
subjdf[2,]$education<-c(NA)
subjdf$education[2]<-mean(subjdf$education,na.rm=TRUE)
##################
subjdf$actmat<-actmat
################plot average activity quantile function ###########
p<-seq(0,1,l=101)
indAD<-which(subjdf$adstatus==1)
actad<-subjdf$actmat[indAD,]
actcontrol<-subjdf$actmat[-indAD,]
cnavg<-colMeans(actcontrol)
adavg<-colMeans(actad)
l1<-loess.smooth(p, cnavg,span = 0.5)
l2<-loess.smooth(p, adavg,span=0.5)
par(mar =c(5.1, 4.1, 4.1, 2.1) )
plot(p,cnavg,type="l",col="blue",xlab="Quantile level p",ylab="Average activity quantile function Q(p)",ylim=c(0,3800),cex.lab=1.2)
lines(p,adavg,type="l",col="red",lty=2)  
legend('topleft',c("CNC ","AD") , 
       lty=c(1,2), col=c("blue","red"), bty='n', cex=.95)
```

![](SOTDR_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#distributional modelling using subject-specific quantile functions of PA
library(refund)
p<-seq(0,1,l=101)
fit.lf <- pfr(adstatus~ age+Sex+education+lf(actmat,argvals=p, k=12, bs="ps",m=2),data = subjdf,family="binomial")
summary(fit.lf)
```

    ## 
    ## Family: binomial 
    ## Link function: logit 
    ## 
    ## Formula:
    ## adstatus ~ age + Sex + education + s(x = actmat.tmat, by = L.actmat, 
    ##     k = 12, bs = "ps", m = 2)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 10.58763    4.13933   2.558   0.0105 *  
    ## age         -0.07215    0.04278  -1.687   0.0917 .  
    ## Sex          2.52749    0.62385   4.051 5.09e-05 ***
    ## education   -0.16671    0.09244  -1.803   0.0713 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                           edf Ref.df Chi.sq p-value  
    ## s(actmat.tmat):L.actmat 2.914  3.221  10.63  0.0172 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =   0.33   Deviance explained = 32.6%
    ## -REML =  53.04  Scale est. = 1         n = 92

``` r
aa<-coef(fit.lf,n=101)
par(mar=c(5,5,4,2))
plot(p[71:101],aa$value[71:101],xlab=" Quantile level p",ylab=expression(paste(beta(p))),main="Distributional effect on log odds of AD",ylim=c(-0.08,0.08),type="l",cex.lab=1.2)
abline(h=0)
upper<-as.numeric(aa$value[71:101])+(1.96*aa$se[71:101])
lower<-as.numeric(aa$value[71:101])-(1.96*aa$se[71:101])
lines(aa$acargvals.argvals[71:101],upper,lty=2)
lines(aa$acargvals.argvals[71:101],lower,lty=2)
abline(v=seq(0.75,0.95,by=0.05),lty=3,col="gray20")
```

![](SOTDR_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

We observe that higher maximal quatile levels (*p* ∈ (0.90, 1)) of PA
are significantly associated with lower odds of AD.

# SOTDR

We illustrate the SOTDR method (Section 4.1, Model 4 of the paper) via
time-by-distribution data objects of PA for modelling cognitive status.

``` r
#####Load the data and preprocess for time-by-distribution modelling############
rm(list=ls())
path=getwd()
path=paste(path,"/minute-csv",sep="")
setwd(path)
filelist<-list.files(path)
myfiles <-lapply(filelist, read.csv)
subjdf<-matrix(0,nrow = 92,ncol=5)
subjdf<-as.data.frame(subjdf)
actmatlist<-list()
for(i in 1:92)
{
  mydata<-myfiles[[i]]
  mydata$ADStatus<-ifelse(mydata$ADStatus=="Yes",1,0)
  mydata$Sex<-ifelse(mydata$Sex=="Male",1,0)
  #Date, Time , Vector.Magnitude
  subdata<-mydata[,c("Date","Time","Vector.Magnitude")]
  library(chron)
  subdata$Time<-60 * 24 * as.numeric(times(subdata$Time))
  n<-nlevels(subdata$Date)
  lev<-levels(subdata$Date)
  tbtick<-seq(0,1440,by=10)
  actmat<-matrix(0,nrow=144,ncol=101)
  l=length(tbtick)-1
  p<-seq(0,1,l=101)
  for (j in 1:l)
  {for(k in 1:length(p)){
    
    tempdataj<-subdata[subdata$Time>=tbtick[j]&subdata$Time <tbtick[j+1],]
    actmat[j,k]<-quantile(tempdataj$Vector.Magnitude,probs = p[k],na.rm = TRUE)
  }
  }
  subjdf[i,]<-c(mydata$id[1],mydata$Age[1],mydata$ADStatus[1],mydata$Sex[1],mydata$YearsOfEducation[1])
  actmatlist[[i]]<-actmat
  }

names(subjdf)<-c("id","age","adstatus","Sex","education")
summary(subjdf)
summary(subjdf)
subjdf[2,]$education<-c(NA)
subjdf$education[2]<-mean(subjdf$education,na.rm=TRUE)
###SAVING THE EXTRACTED DATA
save(subjdf,file="subjdf.Rdata")
save(actmatlist,file="twodimact.Rdata")
```

``` r
################Load pre-extracted data#####################
load("twodimact.Rdata")
load("subjdf.Rdata")
subjdf$education[2]<-mean(subjdf$education,na.rm=TRUE)
tbtick<-seq(0,1440,by=10)
binmid<-tbtick+5
binmid<-binmid[-145]
p<-seq(0,1,l=101)
####################plot average TD objects and their difference###############
indAD<-which(subjdf$adstatus==1)
actadlist<-actmatlist[indAD]
actcontrollist<-actmatlist[-indAD]
cnavg<-Reduce("+", actcontrollist) / length(actcontrollist)
adavg<-Reduce("+", actadlist) / length(actadlist)
library(fields)
par(mfrow=c(2,2))
image.plot(binmid/60,p,adavg,xaxt="n",xlab="Time of the Day", ylab = "p",main = "Average Q(t,p) surface AD ",cex.axis=1.25,
           cex.lab=1.25,cex.main=1.25,
           legend.shrink=0.75,legend.line=-1.5,zlim=c(0,3800),xlim=c(0,24))
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE,cex = 1.2)


image.plot(binmid/60,p,cnavg,xaxt="n",xlab="Time of the Day", ylab = "p",main = "Average Q(t,p) surface CNC",cex.axis=1.25,
           cex.lab=1.25,cex.main=1.25,
           legend.shrink=0.75,legend.line=-1.5,zlim=c(0,3800))
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE,cex = 1.2)



diffactmat<-cnavg-adavg
image.plot(binmid/60,p,diffactmat,xaxt="n",xlab="Time of the Day", ylab = "p",main = "difference surface Q(t,p) between CNC and AD",cex.axis=1.25,
           cex.lab=1.25,cex.main=1.25,
           legend.shrink=0.75,legend.line=-1.5)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE,cex = 1.2)
```

![](SOTDR_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

The model fitting prcedure is illustrated below.

``` r
#start modelling two dim SOFR
##create the basis in two dimesions
n<-nrow(subjdf)
T<-binmid/60
library(fda)
knotsT = seq(min(T),max(T),l=10)
knotsP = seq(min(p),max(p),l=10)
norder = 4
# this implies the number of basis functions
nbasisT = length(knotsT) + norder - 2
nbasisP = length(knotsP) + norder - 2
dayrngT = c(min(T),max(T))
dayrngP = c(0,1)
bbasisT = create.bspline.basis(dayrngT,nbasisT,norder,knotsT)
bbasisP = create.bspline.basis(dayrngP,nbasisP,norder,knotsP)
deltap<-p[2]-p[1]
deltat<-T[2]-T[1]

Wlist<-list()
for (i in 1:n)
{Wlist[[i]]<-matrix(0,nrow=nbasisT,ncol = nbasisP)
for(k in 1:nbasisT )
{for(l in 1: nbasisP)
{meat<-actmatlist[[i]]
tbas<-eval.basis(T,bbasisT)[,k]
pbas<-  eval.basis(p,bbasisP)[,l]
val<-(deltap*deltat)*as.numeric(t(tbas)%*%meat%*%pbas)
Wlist[[i]][k,l]<-val
}
}
}
Wiklmat<-matrix(0,nrow =n,ncol = nbasisP*nbasisT)
for (i in 1:n)
{Wiklmat[i,]<-as.vector(t(Wlist[[i]]))
}

newdf<-subjdf[,c(2,3,4,5)]
newdf<-cbind(newdf,Wiklmat)
names(newdf)[5:148]<-paste("var",1:144)

####step 1 LASSO
x <- model.matrix(adstatus~., newdf)[,-1]
# Convert the outcome (class) to a numerical variable
y <- newdf$adstatus
library(glmnet)
set.seed(123)
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial",nfolds = n)
# Fit the final model on the training data
model <- glmnet(x, y, alpha = 1, family = "binomial",lambda = cv.lasso$lambda.min)
gamma<-model$beta
indsel<-which(gamma!=0) 
varsel<-names(newdf[,-2])[indsel]
#varsel
#### Step 2: fit glm with selected components
pfrfit2<-glm(adstatus~age+Sex+education+`var 60`+`var 72`+`var 96`
             ,family = "binomial",data=newdf)
summary(pfrfit2) 
```

    ## 
    ## Call:
    ## glm(formula = adstatus ~ age + Sex + education + `var 60` + `var 72` + 
    ##     `var 96`, family = "binomial", data = newdf)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.1211  -0.7471  -0.3092   0.7377   1.8839  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 12.3680373  4.5914463   2.694  0.00707 ** 
    ## age         -0.0893098  0.0473063  -1.888  0.05904 .  
    ## Sex          2.6373712  0.6757296   3.903  9.5e-05 ***
    ## education   -0.1743943  0.0954234  -1.828  0.06761 .  
    ## `var 60`    -0.0108673  0.0059868  -1.815  0.06949 .  
    ## `var 72`     0.0005886  0.0073404   0.080  0.93609    
    ## `var 96`    -0.0137289  0.0072087  -1.904  0.05685 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 125.400  on 91  degrees of freedom
    ## Residual deviance:  82.866  on 85  degrees of freedom
    ## AIC: 96.866
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
pfrfit21<-glm(adstatus~age+Sex+education
             ,family = "binomial",data=newdf)
summary(pfrfit21) 
```

    ## 
    ## Call:
    ## glm(formula = adstatus ~ age + Sex + education, family = "binomial", 
    ##     data = newdf)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.8351  -0.8531  -0.4825   0.9068   2.1392  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  4.13738    2.96804   1.394 0.163325    
    ## age         -0.02373    0.03363  -0.706 0.480376    
    ## Sex          1.96063    0.50909   3.851 0.000118 ***
    ## education   -0.22573    0.08557  -2.638 0.008341 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 125.40  on 91  degrees of freedom
    ## Residual deviance: 100.99  on 88  degrees of freedom
    ## AIC: 108.99
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
anova(pfrfit2,pfrfit21,test = "LRT")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: adstatus ~ age + Sex + education + `var 60` + `var 72` + `var 96`
    ## Model 2: adstatus ~ age + Sex + education
    ##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)    
    ## 1        85     82.866                         
    ## 2        88    100.988 -3  -18.122 0.000415 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
###plot estimated beta(t,p)##
gamma<-model$beta
betacoef<-gamma[-c(1:3)]
betacoef[60]<-pfrfit2$coefficients[5]
betacoef[72]<-pfrfit2$coefficients[6]
betacoef[96]<-pfrfit2$coefficients[7]

beta.tp<-function(t,p)
{tbasvec<-eval.basis(t,bbasisT)
pbasvec<-  eval.basis(p,bbasisP)
alpha_klmat<-matrix(betacoef,nrow = 12,ncol = 12,byrow = TRUE)
val<-as.numeric((tbasvec)%*%alpha_klmat%*%t(pbasvec))
val
}
tseq<-T
pseq<-p
betamat<-beta.tp(tseq,pseq)
coefs<- matrix(betamat,nrow = 144,ncol = 101)
coefmat<-matrix(0,nrow = 144,ncol = 101)
for(i in 1:144){
  for(j in 1:101)
  {
    coefmat[i,j]<-beta.tp(tseq[i],pseq[j])
  }
}

colorTable<- designer.colors(39, c( "blue","white", "red") )
brks<- seq( -0.01,0.009,l=40) 
image.plot(binmid/60,p,coefmat,breaks=brks, col=colorTable,xaxt="n",xlab="Time of the Day", ylab = "Quantile level p",main="Time-by-distribution effect on log odds of AD",cex.lab=1.2)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE,cex = 1.2)
abline(v=seq(4,20,by=4),lty=3,col="gray20")
```

![](SOTDR_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Increased maximal capacity of PA during the morning (∼ 7 a.m- 10 a.m)
and in the afternoon (∼ 3 p.m- 5 p.m) is found to be associated with
lower odds of AD.

# SOTDR-L

Next, we illustrate the SOTDR-L approach for modelling cognitive status.
We extract the first 4 diurnal time-varying L-moments and plot them.

``` r
#####Load the data and preprocess for extracting########### 
##########time-varying  L-moments############
rm(list=ls())
path=getwd()
path=paste(path,"/minute-csv",sep="")
setwd(path)
filelist<-list.files(path)
myfiles <-lapply(filelist, read.csv)
subjdf<-matrix(0,nrow = 92,ncol=5)
subjdf<-as.data.frame(subjdf)
actmat1<-matrix(0,nrow = 92,ncol=144)
actmat2<-matrix(0,nrow = 92,ncol=144)
actmat3<-matrix(0,nrow = 92,ncol=144)
actmat4<-matrix(0,nrow = 92,ncol=144)
for(i in 1:92)
{
  mydata<-myfiles[[i]]
  mydata$ADStatus<-ifelse(mydata$ADStatus=="Yes",1,0)
  mydata$Sex<-ifelse(mydata$Sex=="Male",1,0)
  #mydata$Vector.Magnitude<-log(1+mydata$Vector.Magnitude)
  #Date, Time , Vector.Magnitude
  subdata<-mydata[,c("Date","Time","Vector.Magnitude")]
  library(chron)
  subdata$Time<-60 * 24 * as.numeric(times(subdata$Time))
  n<-nlevels(subdata$Date)
  lev<-levels(subdata$Date)
  tbtick<-seq(0,1440,by=10)
  actvec1<-c()
  actvec2<-c()
  actvec3<-c()
  actvec4<-c()
  l=length(tbtick)-1
  library(lmom)
  for (j in 1:l)
  {
  tempdataj<-subdata[subdata$Time>=tbtick[j]&subdata$Time <tbtick[j+1],]
  #calculate first 4 L-moments
    tempmu<-samlmu(tempdataj$Vector.Magnitude,nmom=4,ratios=FALSE) 
    actvec1[j]<-tempmu[1]
    actvec2[j]<-tempmu[2]
    actvec3[j]<-tempmu[3]
    actvec4[j]<-tempmu[4]
    }
subjdf[i,]<-c(mydata$id[1],mydata$Age[1],mydata$ADStatus[1],mydata$Sex[1]
              ,mydata$YearsOfEducation[1])
  actmat1[i,]<-actvec1
  actmat2[i,]<-actvec2
  actmat3[i,]<-actvec3
  actmat4[i,]<-actvec4
}
names(subjdf)<-c("id","age","adstatus","Sex","education")
#summary(subjdf)
subjdf[2,]$education<-c(NA)
##################
subjdf$actmat1<-actmat1
subjdf$actmat2<-actmat2
subjdf$actmat3<-actmat3
subjdf$actmat4<-actmat4
#####################Plotting diurnal L-moments####################
subjdf$education[2]<-mean(subjdf$education,na.rm=TRUE)
tbtick<-seq(0,1440,by=10)
binmid<-tbtick+5
binmid<-binmid[-145]
indAD<-which(subjdf$adstatus==1)
par(mfrow=c(2,2))
par(mar = c(4.5, 4.5, 2, 2)) 
act1ad<-subjdf$actmat1[indAD,]
act1control<-subjdf$actmat1[-indAD,]
cnavg1<-colMeans(act1control)
adavg1<-colMeans(act1ad)
l1<-loess.smooth(binmid, cnavg1,span = 0.25)
l2<-loess.smooth(binmid, adavg1,span = 0.25)
plot(l1$x/60,l1$y,type="l",col="blue",xaxt="n",xlab="Time of the day",ylab="L1(t)",cex.lab=1.1)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE)
lines(l2$x/60,l2$y,type="l",col="red",lty=2)  
legend('topleft',c("CNC ","AD") , 
       lty=c(1,2), col=c("blue", "red"), bty='n', cex=.95)

##l2
act2ad<-subjdf$actmat2[indAD,]
act2control<-subjdf$actmat2[-indAD,]
cnavg2<-colMeans(act2control)
adavg2<-colMeans(act2ad)
l1<-loess.smooth(binmid, cnavg2,span = 0.25)
l2<-loess.smooth(binmid, adavg2,span = 0.25)
plot(l1$x/60,l1$y,type="l",col="blue",xaxt="n",xlab="Time of the day",ylab="L2(t)",cex.lab=1.1)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE)
lines(l2$x/60,l2$y,type="l",col="red",lty=2)  
legend('topleft',c("CNC ","AD") , 
       lty=c(1,2), col=c("blue", "red"), bty='n', cex=.95)
##l3
act3ad<-subjdf$actmat3[indAD,]
act3control<-subjdf$actmat3[-indAD,]
cnavg3<-colMeans(act3control)
adavg3<-colMeans(act3ad)
l1<-loess.smooth(binmid, cnavg3,span = 0.25)
l2<-loess.smooth(binmid, adavg3,span = 0.25)
plot(l1$x/60,l1$y,type="l",col="blue",xaxt="n",xlab="Time of the day",ylab="L3(t)",cex.lab=1.1)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE)
lines(l2$x/60,l2$y,type="l",col="red",lty=2)  
legend('topleft',c("CNC ","AD") , 
       lty=c(1,2), col=c("blue", "red"), bty='n', cex=.95)
##l4
act4ad<-subjdf$actmat4[indAD,]
act4control<-subjdf$actmat4[-indAD,]
cnavg4<-colMeans(act4control)
adavg4<-colMeans(act4ad)
l1<-loess.smooth(binmid, cnavg4,span = 0.25)
l2<-loess.smooth(binmid, adavg4,span = 0.25)
plot(l1$x/60,l1$y,type="l",col="blue",xaxt="n",xlab="Time of the day",ylab="L4(t)",ylim=c(10,60),cex.lab=1.1)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE)
lines(l2$x/60,l2$y,type="l",col="red",lty=2)  
legend('topleft',c("CNC ","AD") , 
       lty=c(1,2), col=c("blue", "red"), bty='n', cex=.95)
```

![](SOTDR_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

We fit the SOTDR-L model as illustrated in Section 4.2 of the paper.

``` r
#Modelling cognitive status

library(refund)
#use variable selection with gel
library(refund)
df1<-subjdf$actmat1
df2<-subjdf$actmat2
df3<-subjdf$actmat3
df4<-subjdf$actmat4
w.sm1 = fpca.sc(df1,pve=0.95, var=TRUE,argvals =binmid/60,nbasis = 10 )
w.sm2 = fpca.sc(df2,pve=0.95, var=TRUE,argvals =binmid/60,nbasis = 10 )
w.sm3 = fpca.sc(df3,pve=0.95, var=TRUE,argvals =binmid/60,nbasis = 10 )
w.sm4 = fpca.sc(df4,pve=0.95, var=TRUE,argvals =binmid/60,nbasis = 10 )

scmat1<-w.sm1$scores
scmat2<-w.sm2$scores
scmat3<-w.sm3$scores
scmat4<-w.sm4$scores
newdf<-subjdf[,c(2,3,4,5)]
trainnew<-cbind(newdf,scmat1,scmat2,scmat3,scmat4)
names(trainnew)[5:29]<-paste("sc",1:25,sep = "")

group<-c(0,0,0,rep(1,ncol(scmat1)),rep(2,ncol(scmat2)),rep(3,ncol(scmat3)),rep(4,ncol(scmat4)))
Groupvar<-as.factor(group)

# Dumy code categorical predictor variables
x <- model.matrix(adstatus~., trainnew)[,-1]
# Convert the outcome (class) to a numerical variable
y <- trainnew$adstatus
library(grpreg)
fit<- grpreg(x,y,Groupvar,penalty="gel",family="binomial",nlambda = 10000)
cvfit<-select(fit,crit="EBIC",df="active")
gamma<-cvfit$beta
indsel<-which(gamma[-1]!=0)
Groupsel<-group[indsel]
grind<-unique(Groupsel) #3
grind #3rd order Lmoment selected
```

    ## [1] 0 3

``` r
fit.lf2 <- pfr(adstatus~ age+Sex+education+lf(actmat3,argvals=binmid/60, k=7, bs="ps",m=2),data = subjdf,family="binomial")
aa<-coef(fit.lf2,n=144)

######plot diurnal effect of L3(t)#################
par(mar=c(5,5,4,2))
plot(binmid/60,aa$value,xaxt="n",xlab="Time of the day (t)",ylab=expression(paste(beta(t))),main="Estimated temporal effect of
L3(t) on log odds of AD",type="l",ylim=c(-.006,0.001),cex.lab=1.2)
xtick<-seq(0, 24, by=4)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 1, xpd = TRUE,cex = 1.2)
abline(h=0)
upper<-as.numeric(aa$value)+1.96*aa$se
lower<-as.numeric(aa$value)-1.96*aa$se
lines(aa$acargvals3.argvals,upper,lty=2)
lines(aa$acargvals3.argvals,lower,lty=2)
```

![](SOTDR_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

We observe that an increase in the value of third order L-moment of
physical activity, during the window (8 a.m- 6 p.m) is associated with a
lower odds of AD.

<br><br><br>
