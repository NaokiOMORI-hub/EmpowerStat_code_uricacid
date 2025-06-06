## Purpose: This script runs a smooth-curve-fitting approach based on a generalized additive model using Empower Stat statistical analysis software.
## Author: Naoki Omori                                             
## naoki.shimane.medical@gmail.com                                      
## Department of Neurology, Shimane University
## Izumo-shi, Shimane, Japan

# Smooth-curve-fitting model：Caudate ncleus head_Left
Sys.setlocale(category = 'LC_ALL', locale = 'English_United States.1252'); 
.libPaths(file.path(R.home(),'library')); 
library(doBy); 
options(timeout=600); 
library(plotrix); 
library(stringi); 
library(stringr); 
library(survival); 
library(rms); 
library(nnet); 
library(car); 
library(mgcv); 
pdfwd<-6; pdfht<-6; 
load('C:/EmpowerXYS/A123/Analysis/data/UAstandardfibnegative.Rdata'); 
if (length(which(ls()=='EmpowerStatsR'))==0) EmpowerStatsR<-get(ls()[1]); 
names(EmpowerStatsR)<-toupper(names(EmpowerStatsR)); 
originalVNAME<-names(EmpowerStatsR); 
ofname<-'data_1_tbl'; 
vname<-c(NA,'AGE','SEX','SEX.0','SEX.1','CAUDATEHEAD.LEFT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','ETAT.CLIBRE.0','ETAT.CLIBRE.1','ETAT.CLIBRE.2','ETAT.CLIBRE.3','HT','HT.0','HT.1','DM','DM.0','DM.1','DL','DL.0','DL.1')[-1]; 
vlabel<-c(NA,'AGE','SEX','  0','  1','CAUDATEHEAD.LEFT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','  0','  1','  2','  3','HT','  0','  1','DM','  0','  1','DL','  0','  1')[-1]; 
varused4this <- c('AGE','SEX','CAUDATEHEAD.LEFT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','HT','DM','DL'); 
pkgs<-c('mgcv','gdata'); 
for (g in pkgs) {  
if (!(g %in% rownames(installed.packages()))) install.packages(g,repos='https://cloud.r-project.org'); 
}
library(mgcv); 
library(gdata); 
WD <- EmpowerStatsR; rm(EmpowerStatsR); gc(); 
title<-'Spline smoothing plot'; 
wd.subset=''; 
weights.var<- NA; 
yvname<-c('CAUDATEHEAD.LEFT_ROI'); 
ydist<-c('gaussian'); 
ylink<-c('identity'); 
ylv<-c(0); 
par1<-1; 
xvname<-c('AGE','SEX','BMI','EGFR','HT','DM','DL','ETAT.CLIBRE'); 
sxf<-c(0,0,0,0,0,0,0,0); 
xlv<-c(0,2,0,0,2,2,2,4); 
svname<-c('UA'); 
sdf<-c(0); 
slv<-c(0); 
par3<-1; 
timevar<- NA; 
vname.start<- NA; 
vname.stop<- NA; 
subjvname<- NA; 
colvname<- NA; 
chk<- 1; 
dec<-4; 
##R package## mgcv gdata ##R package##;
vec2shift<-function(vnew,vorg,f, opt) {
  if (is.na(f[1])) {
    mean1<-mean(vorg)
    if (opt=="logit") mean1<-log(mean1/(1-mean1))
    if (opt=="log") mean1<-log(mean1)
    vnew<-vnew+(mean1-mean(vnew))
  } else {
    mean1<-tapply(vorg,factor(f),mean)
    if (opt=="logit") mean1<-log(mean1/(1-mean1))
    if (opt=="log") mean1<-log(mean1)
    mean2<-tapply(vnew,factor(f),mean);  meand<-mean1-mean2;  lvf<-levels(factor(f))
    for (z in (1:length(lvf))) {vnew[factor(f)==lvf[z]]<-vnew[factor(f)==lvf[z]]+meand[z]; }
  }
  return(vnew)
}
getNumber<-function(str, n) {
  str<-substr(str,2,nchar(str)-1)
  for (i in (1:nchar(str))) {if (substr(str,i,i)==",") {p=i; break}; }
  ifelse(n==1,return(substr(str,1,p-1)),return(substr(str,p+1,nchar(str))))
}
legLocate<-function(x,y) {
  x[is.infinite(y)]<-NA
  y[is.infinite(y)]<-NA
  xmin<-min(x,na.rm=TRUE); xmax<-max(x,na.rm=TRUE)
  ymin<-min(y,na.rm=TRUE); ymax<-max(y,na.rm=TRUE)
  yoff<-(ymax-ymin); tmp<-table(cut(x,3),cut(y,4))
  tmp.r=which.min(tmp[,4]);tmp.c=4
  if (tmp[2,1]==0) {tmp.r=2;tmp.c=1}
  if (tmp[1,1]==0) {tmp.r=1;tmp.c=1}
  if (tmp[3,1]==0) {tmp.r=3;tmp.c=1}
  if (tmp[2,4]==0) {tmp.r=2;tmp.c=4}
  if (tmp[1,4]==0) {tmp.r=1;tmp.c=4}
  if (tmp[3,4]==0) {tmp.r=3;tmp.c=4}
  pos.y<-colnames(tmp)[tmp.c];  pos.x<-rownames(tmp)[tmp.r];  pct<-0.15
  if (tmp.c==4) {
     if (min(tmp[,4])>0) {pct<-0.3}
     ymax<-ymax+yoff*pct; legy<-ymax;ymin<-ymin-yoff*0.1
  } 
  if (tmp.c==1) {
     if (min(tmp[,1])>0) {pct<-0.3}
     legy<-as.numeric(getNumber(pos.y,2));ymin<-ymin-yoff*pct;ymax=ymax+yoff*0.1
  } 
  legx<-as.numeric(getNumber(pos.x,1))
  return(cbind(xmin,xmax,ymin,ymax,legx,legy))
}
mat2htmltable<-function(mat) {
  t1<- apply(mat,1,function(z) paste(z,collapse="</td><td>"))
  t2<- paste("<tr><td>",t1,"</td></tr>")
  return(paste(t2,collapse=" "))
}
setgam<-function(fml,yi) {
  if (ydist[yi]=="") ydist[yi]<-"gaussian"
  if (ydist[yi]=="exact") ydist[yi]<-"binomial"
  if (ydist[yi]=="breslow") ydist[yi]<-"binomial"
  if (ydist[yi]=="gaussian") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=gaussian(link="identity"))
  if (ydist[yi]=="binomial") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=binomial(link="logit"))
  if (ydist[yi]=="poisson") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=poisson(link="log"))
  if (ydist[yi]=="gamma") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=Gamma(link="inverse"))
  if (ydist[yi]=="negbin") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=negbin(c(1,10), link="log"))
  return(mdl)
}
gam2htmltable<-function(mdl) {
  gs<-summary(mdl)
  np<-length(gs$p.coeff)
  coe<-gs$p.table
  if (gs$family[[2]]=="log" | gs$family[[2]]=="logit") {
    cnames<-c(colnames(coe),"exp(est)","95%CI low","95%CI upp")
    coe<- cbind(coe, exp(coe[,1]), exp(coe[,1]-1.96*coe[,2]), exp(coe[,1]+1.96*coe[,2]))
  }
  if (gs$family[[2]]=="identity") {
    cnames<-c(colnames(coe),"95%CI low","95%CI upp")
    coe<- cbind(coe, coe[,1]-1.96*coe[,2], coe[,1]+1.96*coe[,2])
  }
  oo1<-cbind(c("",rownames(coe)),rbind(cnames,round(coe,dec)))
  oo<-c("</br>Linear terms effect<table border=3>",mat2htmltable(oo1),"</table>")
  if (!is.null(gs$pTerms.table)) {
    xsq<-gs$pTerms.table
    oo2<-cbind(c("",rownames(xsq)),rbind(colnames(xsq),round(xsq,dec)))
    oo<-c(oo, "</br>Chi-square tests for linear terms<table border=3>",mat2htmltable(oo2),"</table>")
  }
  if (!is.null(gs$s.table)) {
   stb<-gs$s.table
   oo3<-cbind(c("",rownames(stb)),rbind(colnames(stb),round(stb,dec)))
   oo<-c(oo, "</br>Approximate significance of smooth terms<table border=3>",mat2htmltable(oo3),"</table>")
  }
  p0<-c("N:", gs$n)
  p1<-c("Adj. r-square:", round(gs$r.sq,4))
  p2<-c("Deviance explained:", round(gs$dev.expl,4))
  p3<-c("UBRE score (sp.criterion):", round(gs$sp.criterion,4))
  p4<-c("Scale estimate:", gs$scale)
  p5<-c("family:", gs$family[[1]])
  p6<-c("link function:", gs$family[[2]])
  oo4<-rbind(p0,p1,p2,p3,p4,p5,p6)
  oo<-c(oo, "</br>Model statistics<table border=3>",mat2htmltable(oo4),"</table>")
  return(oo)
}
gam2pngs<-function(mdl,yi,xi) {
  pred<-predict.gam(mdl,type="terms",se.fit=TRUE)
  mfit<-NA; sfit<-NA; tmp.cname<-NA; kk0<-NA; mfit0<-NA;  ww0<-NULL
  if (xi==0) {kb=1; ke=ns;} else {kb=xi; ke=xi;}
  for (k in (kb:ke)) {
    if (slv[k]==0) {
      mfit<-cbind(mfit,apply(cbind(0,pred$fit[,sxterms[k,]]),1,sum));
      sfit<-cbind(sfit,apply(cbind(0,pred$se.fit[,sxterms[k,]]),1,sum));
      tmp.cname<-c(tmp.cname,svname[k])
      kk0<-c(kk0,k)
    }
  }
  tmp.cname<-tmp.cname[-1];  kk0<-kk0[-1]
  mfit<-matrix(mfit[,-1],ncol=length(tmp.cname)); colnames(mfit)<-paste(tmp.cname,".fit",sep="");
  sfit<-matrix(sfit[,-1],ncol=length(tmp.cname)); colnames(sfit)<-paste(tmp.cname,".se",sep="");
  if (!is.na(colvname)) {tmpfac<-wd[,colvname];} else {tmpfac<-NA;} 
  if (mdl$family[2]=="logit") {
    mfit0<-mfit;  mfit0.low<-mfit0-1.96*sfit;  mfit0.upp<-mfit0+1.96*sfit;
    colnames(mfit0.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit0.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit0)<-paste(tmp.cname,".fit",sep="");
    ww0<-cbind(wd,mfit0); if (is.na(colvname)) ww0<-cbind(ww0,mfit0.low,mfit0.upp)
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac,"logit"))
    mfit.low<-mfit-1.96*sfit; mfit.low<-matrix(exp(mfit.low)/(1+exp(mfit.low)),ncol=length(tmp.cname))
    mfit.upp<-mfit+1.96*sfit; mfit.upp<-matrix(exp(mfit.upp)/(1+exp(mfit.upp)),ncol=length(tmp.cname))
    mfit<-matrix(exp(mfit)/(1+exp(mfit)),ncol=length(tmp.cname))
    colnames(mfit.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit)<-paste(tmp.cname,".fit",sep="");
    ww<-cbind(wd,mfit); if (is.na(colvname)) ww<-cbind(ww,mfit.low,mfit.upp)
  } else if (mdl$family[2]=="log") {
    mfit0<-mfit;  mfit0.low<-mfit0-1.96*sfit;  mfit0.upp<-mfit0+1.96*sfit;
    colnames(mfit0.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit0.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit0)<-paste(tmp.cname,".fit",sep="");
    ww0<-cbind(wd,mfit0); if (is.na(colvname)) ww0<-cbind(ww0,mfit0.low,mfit0.upp)
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac,"log"))
    mfit.low<-mfit-1.96*sfit; mfit.low<-matrix(exp(mfit.low),ncol=length(tmp.cname))
    mfit.upp<-mfit+1.96*sfit; mfit.upp<-matrix(exp(mfit.upp),ncol=length(tmp.cname))
    mfit<-matrix(exp(mfit),ncol=length(tmp.cname))
    colnames(mfit.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit)<-paste(tmp.cname,".fit",sep="");
    ww<-cbind(wd,mfit); if (is.na(colvname)) ww<-cbind(ww,mfit.low,mfit.upp)
  } else if (mdl$family[2]=="identity") {
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac," "))
    ww<-cbind(wd,mfit,sfit)
  } else {
    ww<-cbind(wd,mfit,sfit)
  }
  if (xi!=0) {xf<-paste(ofname,yvname[yi],svname[xi],"gam.xls",sep="_");} else {xf<-paste(ofname,yvname[yi], "gam.xls", sep="_");}
  write.table(ww,file=xf,row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE,quote=FALSE)
  if (!is.null(ww0)) {
    if (xi!=0) {xf<-paste(ofname,yvname[yi],svname[xi],"gam0.xls",sep="_");} else {xf<-paste(ofname,yvname[yi], "gam0.xls", sep="_");}
    write.table(ww0,file=xf,row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE,quote=FALSE)
  }
  px<-c(20,1:9); gg<-"";
  for (k in kk0) {
     cname1<-paste(svname[k],".fit",sep="");  y.tmp<-mfit[,cname1]; y0.tmp<-NULL
     if (mdl$family[2]=="logit" | mdl$family[2]=="log") {
       cname2<-paste(svname[k],".low",sep=""); y.low<-mfit.low[,cname2]
       cname3<-paste(svname[k],".upp",sep=""); y.upp<-mfit.upp[,cname3]
	   y0.tmp<-mfit0[,cname1]; y0.low<-mfit0.low[,cname2]; y0.upp<-mfit0.upp[,cname3]
     } else {
       cname2<-paste(svname[k],".se",sep=""); se.tmp<-sfit[,cname2]; 
       y.low<-y.tmp-1.96*se.tmp; y.upp<-y.tmp+1.96*se.tmp
     }
     x.tmp<-wd[,svname[k]]; 
     pngf<-paste(ofname,yvname[yi],svname[k],"smooth.png",sep="_")
     pdff<-paste(ofname,yvname[yi],svname[k],"smooth.pdf",sep="_")
     pngf0<-paste(ofname,yvname[yi],svname[k],"smooth1.png",sep="_")
     pdff0<-paste(ofname,yvname[yi],svname[k],"smooth1.pdf",sep="_")
     pngf2<-paste(ofname,yvname[yi],svname[k],"smooth2.png",sep="_")
     pdff2<-paste(ofname,yvname[yi],svname[k],"smooth2.pdf",sep="_")
     if (is.na(colvname)) {
       if (chk==1) {tmp.col<-c("red","blue");} else {tmp.col<-rep("black",2);}
       xy<-legLocate(c(x.tmp,x.tmp),c(y.low,y.upp))
       png(pngf,width=720,height=560)
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       dev.off()

       pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       dev.off()

       png(pngf0,width=720,height=560)
       tmp.ord<-order(x.tmp);  x.tmp0<-x.tmp[tmp.ord]; 
       y.tmp0<-y.tmp[tmp.ord]; y.low0<-y.low[tmp.ord]; y.upp0<-y.upp[tmp.ord]
       plot(y.tmp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       rug(x.tmp0)
       dev.off()

       pdf(pdff0,width=pdfwd, height=pdfht, family="GB1");
       tmp.ord<-order(x.tmp);  x.tmp0<-x.tmp[tmp.ord]; 
       y.tmp0<-y.tmp[tmp.ord]; y.low0<-y.low[tmp.ord]; y.upp0<-y.upp[tmp.ord]
       plot(y.tmp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       rug(x.tmp0)
       dev.off()
	   
       if (!is.null(y0.tmp)) {
	     xy0<-legLocate(c(x.tmp,x.tmp),c(y0.low,y0.upp))
		 if (mdl$family[2]=="logit") {yb0<-paste("Log OR of", yb.utf8[yi]);} else {yb0<-paste("Log RR of", yb.utf8[yi]);}
	     png(pngf2,width=720,height=560)
         tmp.ord<-order(x.tmp);    x.tmp0<-x.tmp[tmp.ord]; 
         y0.tmp0<-y0.tmp[tmp.ord]; y0.low0<-y0.low[tmp.ord]; y0.upp0<-y0.upp[tmp.ord]
         plot(y0.tmp0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
         par(new=TRUE); 
         plot(y0.low0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
         par(new=TRUE); 
         plot(y0.upp0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb0, xlab=sb[k])
         rug(x.tmp0)
         dev.off()
	   }

       rm(tmp.ord,y.tmp0,x.tmp0,y.low0,y.upp0)
       xy<-legLocate(x.tmp,wd[,1])
       pngf1<-paste(ofname,yvname[yi],svname[k],"scatter.png",sep="_")
       pdff1<-paste(ofname,yvname[yi],svname[k],"scatter.pdf",sep="_")

       png(pngf1,width=720,height=560)
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1], type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(wd[,1]~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),type="p",pch=1,cex=0.5, ylab=yb.utf8[yi], xlab=sb.utf8[k])    
       dev.off()

       pdf(pdff1,width=pdfwd, height=pdfht, family="GB1");
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1], type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(wd[,1]~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),type="p",pch=1,cex=0.5, ylab=yb.utf8[yi], xlab=sb.utf8[k])    
       dev.off()

     } else {
       if (chk==1) {tmp.col<-rainbow(ncg);tmp.col1<-c("red","blue")} else {tmp.col<-rep("black",ncg);tmp.col1<-c("black","black")}
       for (b in (1:ncg)) {
         y00<-y.tmp[wd[,colvname]==colv.lv[b]]; x00<-x.tmp[wd[,colvname]==colv.lv[b]]
         xy1<-legLocate(x00,y00)
         pngf1<-paste(ofname,yvname[yi],svname[k],colvname,colv.lv[b],"smooth.png",sep="_")
         pdff1<-paste(ofname,yvname[yi],svname[k],colvname,colv.lv[b],"smooth.pdf",sep="_")

         png(pngf1,width=720,height=560)
         plot(y00~x00,ylim=c(xy1[3],xy1[4]),xlim=c(xy1[1],xy1[2]),col=tmp.col1[1],type="p", pch=20, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         dev.off()

         pdf(pdff1,width=pdfwd, height=pdfht, family="GB1");
         plot(y00~x00,ylim=c(xy1[3],xy1[4]),xlim=c(xy1[1],xy1[2]),col=tmp.col1[1],type="p", pch=20, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         dev.off()

       }
       xy<-legLocate(x.tmp,y.tmp)
       png(pngf,width=720,height=560)
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         if (b>1) par(new=TRUE)
         plot(y0~x0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="p", pch=px[b], ylab=yb.utf8[yi], xlab=sb.utf8[k])
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8, pch=px[1:ncg],bty="n",col=tmp.col)
       dev.off()

       pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         if (b>1) par(new=TRUE)
         plot(y0~x0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="p", pch=px[b], ylab=yb[yi], xlab=sb[k])
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8, pch=px[1:ncg],bty="n",col=tmp.col)
       dev.off()

       png(pngf0,width=720,height=560)
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         tmp.ord<-order(x0); x00<-x0[tmp.ord]; y00<-y0[tmp.ord];
         if (b>1) par(new=TRUE)
         plot(y00~x00,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         rm(tmp.ord,x00,y00)
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
       dev.off()

       pdf(pdff0,width=pdfwd, height=pdfht, family="GB1");
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         tmp.ord<-order(x0); x00<-x0[tmp.ord]; y00<-y0[tmp.ord];
         if (b>1) par(new=TRUE)
         plot(y00~x00,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         rm(tmp.ord,x00,y00)
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
       dev.off()

       if (!is.null(y0.tmp)) {
	     xy0<-legLocate(c(x.tmp,x.tmp),c(y0.low,y0.upp))
		 if (mdl$family[2]=="logit") {yb0<-paste("Log OR of", yb.utf8[yi]);} else {yb0<-paste("Log RR of", yb.utf8[yi]);}
	     png(pngf2,width=720,height=560)
         for (b in (1:ncg)) {
           y0<-y0.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
           tmp.ord<-order(x0); x00<-x0[tmp.ord];  y00<-y0[tmp.ord];
           if (b>1) par(new=TRUE)
           plot(y00~x00,ylim=c(xy0[3],xy0[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb0, xlab=sb.utf8[k])
           rm(tmp.ord,x00,y00)
         }
         legend(xy0[5],xy0[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
         dev.off()
	   }

     }
     gg<-c(gg,paste("<td>",yb[yi]," vs. ",sb[k],"</br><a href=\"",pngf,"\" target=_BLANK><img src=\"",pngf,"?time=\" width=320,height=320></a></td>",sep=""))
  }
  return(gg)
}
adjmean<-function(mdl, yi, xi) {
  if (length(xvname)>0) {allvname<- c(xvname,svname);  all.lv<-c(xlv,slv);} else {allvname<-c(svname); all.lv<-slv;}
  if (!is.na(colvname)) allvname<-c(allvname,colvname)
  nv = length(allvname)
  xi.lv <- levels(factor(WD[,svname[xi]]))
  newd0 <- matrix(0,ncol=nv,nrow=length(xi.lv))
  colnames(newd0)<-allvname
  for (b in 1:length(all.lv)) {
    if (all.lv[b]==0) {
      newd0[,b]<-mean(WD[,allvname[b]],na.rm=TRUE)
    } else {
      uniqv <- unique(WD[,allvname[b]])
      newd0[,b]<-uniqv[which.max(tabulate(match(WD[,allvname[b]], uniqv)))]
    }
  }
  newd0[,svname[xi]]<-as.numeric(xi.lv)
  if (!is.na(colvname)) {
    for (b in 1:ncg) {
      newd1<-newd0; newd1[,colvname]<-as.numeric(colv.lv[b]); 
      if (b==1) {newD<-newd1;} else {newD<-rbind(newD,newd1);}
    }
    f<-table(wd[,svname[xi]],wd[,colvname])
  } else {
    newD<-newd0; ncg<-1;
    f<-table(wd[,svname[xi]])
  }
  pred<-predict(mdl, data.frame(newD),se.fit=TRUE)
  meany.pop0<-tapply(wd[,yvname[yi]],wd[,svname[xi]],function(z) mean(z,na.rm=TRUE))[1]
  if (ylink[yi]=="logit") meany.pop0<-log(meany.pop0/(1-meany.pop0));
  if (ylink[yi]=="log") meany.pop0<-log(meany.pop0);
  shift<-meany.pop0-pred$fit[1]
  y.fit <- pred$fit+shift 
  y.low <- pred$fit+shift-pred$se.fit*1.96
  y.upp <- pred$fit+shift+pred$se.fit*1.96
  y.pred<- cbind(y.fit,y.low,y.upp); 
  tmp.ylab<-paste("Mean of", yb.utf8[yi]);
  cname.pred<-c("Mean","Mean.low","Mean.upp")
  if (ylink[yi]=="logit") {
    y.pred<-exp(y.pred); y.pred<-y.pred/(1+y.pred);
    tmp.ylab<-paste("% of", yb.utf8[yi]);
    cname.pred<-c("Rate","Rate.low","Rate.upp")  
  }
  if (ylink[yi]=="log") {
    y.pred<-exp(y.pred); tmp.ylab<-paste("% of", yb.utf8[yi]);
    cname.pred<-c("Rate","Rate.low","Rate.upp")  
  }
  y.fit<-y.pred[,1]; y.low<-y.pred[,2]; y.upp<-y.pred[,3]
  if (length(xvname)>0) {
    tmp.ylab<-paste("Adjusted", tmp.ylab)
    cname.pred<-paste("adj.", cname.pred, sep="")
  }
  y.pred<-round(y.pred,dec)
  colnames(y.pred)<-cname.pred
  y.pred <-cbind(newD[,svname[xi]],y.pred); colnames(y.pred)[1]<-svname[xi]
  if (!is.na(colvname)) {
    y.pred<-cbind(newD[,colvname],y.pred); colnames(y.pred)[1]<-colvname
  }
  y.pred<-rbind(colnames(y.pred),y.pred)
  px<-c(20,1:9)
  if (chk==1) {tmp.col<-rainbow(ncg);tmp.col1<-c("red","blue");
  } else {tmp.col<-rep("black",ncg);tmp.col1<-tmp.col}

  pngf<-paste(ofname,yvname[yi],svname[xi],"adjmean.png",sep="_")
  pdff<-paste(ofname,yvname[yi],svname[xi],"adjmean.pdf",sep="_")

  if (ncg>1) {
    xy<-legLocate(newD[,svname[xi]],c(y.fit))
    png(pngf,width=720,height=560)
    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      if (b>1) par(new=TRUE)
      plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b],type="b", pch=px[b],
        ylab=tmp.ylab, xlab=sb.utf8[xi])
    }
    legend(xy[5],xy[6],colv.lb,title=colvb,pch=px[1:ncg],bty="n",col=tmp.col)
    dev.off()

    pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      if (b>1) par(new=TRUE)
      plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b],type="b", pch=px[b],
        ylab=tmp.ylab, xlab=sb.utf8[xi])
    }
    legend(xy[5],xy[6],colv.lb,title=colvb,pch=px[1:ncg],bty="n",col=tmp.col)
    dev.off()

    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      y.lci<-y.low[newD[,colvname]==colv.lv[b]]
      y.uci<-y.upp[newD[,colvname]==colv.lv[b]]
      xy<-legLocate(c(x.tmp,x.tmp),c(y.lci,y.uci))
      png(paste(ofname,yvname[yi],svname[xi],colvname,colv.lv[b],"CI.png",sep="_"),width=720,height=560)
      plotCI(x.tmp,y=y.tmp,li=y.lci,ui=y.uci,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main=paste(colvb.utf8, colv.lb.utf8[b],sep=": "))
      lines(x.tmp,y.tmp,lty=2)
      dev.off()

      pdf(paste(ofname,yvname[yi],svname[xi],colvname,colv.lv[b],"CI.pdf",sep="_"),width=pdfwd, height=pdfht, family="GB1");
      plotCI(x.tmp,y=y.tmp,li=y.lci,ui=y.uci,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main=paste(colvb.utf8, colv.lb.utf8[b],sep=": "))
      lines(x.tmp,y.tmp,lty=2)
      dev.off()

    }
  } else {
    x.tmp<-newD[,svname[xi]]
    xy<-legLocate(c(x.tmp,x.tmp),c(y.low,y.upp))
    png(pngf,width=720,height=560)
    plotCI(x.tmp,y=y.fit,li=y.low,ui=y.upp,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main="Adjusted mean & 95% CI")
    lines(x.tmp,y.fit,lty=2)
    dev.off()

    pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
    plotCI(x.tmp,y=y.fit,li=y.low,ui=y.upp,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main="Adjusted mean & 95% CI")
    lines(x.tmp,y.fit,lty=2)
    dev.off()
  }
  oo<-c("</br>Adjusted mean ",yb[yi], " by ", sb[xi], "<table border=3>",mat2htmltable(y.pred),"</table>")
  gg<-paste("<td>",yb[yi]," vs. ",sb[xi],"</br><a href=\"",pngf,"\" target=_BLANK><img src=\"",pngf,"?time=\" width=320,height=320></a></td>", sep="")
  return(list(oo,gg))
}
if (!is.na(weights.var)) {weights<-WD[,weights.var];} else {weights<-1;}
WD<-cbind(WD,weights);
vlabelN<-(substr(vlabel,1,1)==" ");
vlabelZ<-vlabel[vlabelN];vlabelV<-vlabel[!vlabelN]
vnameV<-vname[!vlabelN];vnameZ<-vname[vlabelN]
ny<-length(yvname); yb<-vlabelV[match(yvname,vnameV)]; yb[is.na(yb)]<-yvname[is.na(yb)]
ns<-length(svname); sb<-vlabelV[match(svname,vnameV)]; sb[is.na(sb)]<-svname[is.na(sb)]
yb.utf8<-yb; Encoding(yb.utf8)<-"UTF-8"
sb.utf8<-sb; Encoding(sb.utf8)<-"UTF-8"
ssf<-rep(",fx=FALSE", ns); ssf[sdf>0]<-paste(",k=",sdf[sdf>0],sep="")
sxStr<-paste("s(",svname,ssf,sep="")
sxStr[slv>0]<-paste("factor(",svname[slv>0],")",sep="")
sxx<-paste("s(",svname,")",sep="")
sxx[slv>0]<-paste("factor(",svname[slv>0],")",sep="")
sxx<-matrix(sxx,ncol=1)
if (!is.na(colvname)) {
  sxStr[slv==0]<-paste(sxStr[slv==0],",by=factor(", colvname, ")",sep="")
  sxStr[slv>0]<-paste(sxStr[slv>0],"*factor(", colvname, ")",sep="")
  colv.lv<-levels(factor(WD[,colvname])); ncg<-length(colv.lv); colvb<-vlabel[vname==colvname]; 
  colv.lb<-vlabelZ[match(paste(colvname,colv.lv,sep="."),vnameZ)]
  colv.lb[is.na(colv.lb)]<-colv.lv[is.na(colv.lb)]
  colvb<-vlabelV[match(colvname,vnameV)]; if (is.na(colvb)) colvb<-colvname;
  colvb.utf8<-colvb; Encoding(colvb.utf8)<-"UTF-8"; colv.lb.utf8<-colv.lb; Encoding(colv.lb.utf8)<-"UTF-8"
  sxplots<-NA; sxterms<-NA
  for (i in (1:ns)) {
    sxplots<-c(sxplots,paste(svname[i],"_",colvname,colv.lv,sep=""));
    sxterms<-rbind(sxterms,paste(sxx[i,],":factor(",colvname,")",colv.lv,sep=""))
  }
  sxplots<-sxplots[-1]; sxterms<-matrix(sxterms[-1,],ncol=ncg)
  sxterms<-cbind(paste("factor(",colvname,")",sep=""),sxterms)
} else {ncg<-1;sxplots<-svname; sxterms<-sxx;}
sxStr[slv==0]<-paste(sxStr[slv==0],")",sep="")
nx<-0
if (length(xvname)>0) {
  if (!is.na(colvname)) {xvname<-xvname[xvname!=colvname];}
  nx<-length(xvname); 
}
if (nx>0) {
  xb<-vlabelV[match(xvname,vnameV)]; xb[is.na(xb)]<-xvname[is.na(xb)];
  xvv<-xvname; xvv[xlv>2]<-paste("factor(",xvname[xlv>2],")",sep="")
  if (!is.na(colvname)) {
    xvv[sxf=="S" | sxf=="s"]<-paste("factor(",colvname,")*",xvv[sxf=="S" | sxf=="s"],sep="")
  } 
  xv1<-paste(xvv,collapse="+")
}
if (is.na(par1)) par1<-1
if (ny!=ns & par1==2) par1<-1
if (par1==3) {nterms=ns*ncg*15+nx;} else {nterms=ncg*15+nx;}
w<-c("<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'></head><body>")
w<-c(w,paste("<h2>", title, "</h2>"))
wtab<-"</br></br>Generalize additive models</br>"
wpng<-"</br><table>";
for (i in (1:ny)) {
  if (par1!=3) {
    wtmp<-"";
    if (par1==2) {jstart<-i; jstop<-i;} else {jstart<-1; jstop<-ns;}
    for (j in (jstart:jstop)) {
      tmp.xx<-c(yvname[i],svname[j])
      if (nx>0) tmp.xx<-c(tmp.xx,xvname)
      if (!is.na(colvname[1])) tmp.xx<-c(tmp.xx,colvname)
      wd<-WD[,tmp.xx];
      wd<-wd[apply(is.na(wd),1,sum)==0,]
      fml<-paste(yvname[i],"~",sxStr[j],sep="")
      if (!is.na(colvname)) fml<-paste(fml,"+factor(",colvname,")",sep="")
      if (nx>0) fml<-paste(fml,"+",xv1,sep="")
      tmp.gam<-setgam(fml,i)
      wtab<-c(wtab,paste("</br></br>Outcome:",yb[i]))
      wtab<-c(wtab,paste("</br>Exposure:",sb[j]))
      wtab<-c(wtab,gam2htmltable(tmp.gam))
      if (slv[j]==0) {
        wtmp<-c(wtmp,gam2pngs(tmp.gam,i,j))
      } else {
        stmp<-adjmean(tmp.gam,i,j)
        wtmp<-c(wtmp,stmp[[2]])
        wtab<-c(wtab,stmp[[1]])
      }
    }
    wpng<-c(wpng,"<tr>",wtmp,"</tr>")
  } else {
      tmp.xx<-c(yvname[i],svname);
      if (nx>0) tmp.xx<-c(tmp.xx,xvname)
      if (!is.na(colvname[1])) tmp.xx<-c(tmp.xx,colvname)
      wd<-WD[,tmp.xx];
      wd<-wd[apply(is.na(wd),1,sum)==0,]
      fml<-paste(yvname[i],"~",paste(sxStr,collapse="+"),sep="")
      if (!is.na(colvname)) fml<-paste(fml,"+factor(",colvname,")",sep="")
      if (nx>0) fml<-paste(fml,"+",xv1,sep="")
      tmp.gam<-setgam(fml,i)
      wtab<-c(wtab,paste("</br></br>Outcome:",yb[i]))
      wtab<-c(wtab,gam2htmltable(tmp.gam))
      if (sum(slv==0)>0) {   
        wpng<-c(wpng,"<tr>",gam2pngs(tmp.gam,i,0),"</tr>")
      }
      if (sum(slv>0)>0) {   
        for (k in 1:ns) {
          if (slv[k]>0) {
            stmp<-adjmean(tmp.gam,i,k)
            wtab<-c(wtab,stmp[[1]])
            wpng<-c(wpng,stmp[[2]])
          }
        }
      }
  }
}
wpng<-c(wpng,"</table>")
w<-c(w,wpng,wtab)
w<-c(w,wd.subset)
w<-c(w,paste("</br>Created by EmpowerStats (www.empowerstats.com) and R on",Sys.Date()))
w<-c(w,"</body></html>")
fileConn<-file(paste(ofname,".htm",sep="")); writeLines(w, fileConn)

# Smooth-curve-fitting model：Caudate ncleus head_Right
Sys.setlocale(category = 'LC_ALL', locale = 'English_United States.1252'); 
.libPaths(file.path(R.home(),'library')); 
library(doBy); 
options(timeout=600); 
library(plotrix); 
library(stringi); 
library(stringr); 
library(survival); 
library(rms); 
library(nnet); 
library(car); 
library(mgcv); 
pdfwd<-6; pdfht<-6; 
load('C:/EmpowerXYS/A123/Analysis/data/UAstandardfibnegative.Rdata'); 
if (length(which(ls()=='EmpowerStatsR'))==0) EmpowerStatsR<-get(ls()[1]); 
names(EmpowerStatsR)<-toupper(names(EmpowerStatsR)); 
originalVNAME<-names(EmpowerStatsR); 
ofname<-'data_2_tbl'; 
vname<-c(NA,'AGE','SEX','SEX.0','SEX.1','CAUDATEHEAD.RIGHT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','ETAT.CLIBRE.0','ETAT.CLIBRE.1','ETAT.CLIBRE.2','ETAT.CLIBRE.3','HT','HT.0','HT.1','DM','DM.0','DM.1','DL','DL.0','DL.1')[-1]; 
vlabel<-c(NA,'AGE','SEX','  0','  1','CAUDATEHEAD.RIGHT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','  0','  1','  2','  3','HT','  0','  1','DM','  0','  1','DL','  0','  1')[-1]; 
varused4this <- c('AGE','SEX','CAUDATEHEAD.RIGHT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','HT','DM','DL'); 
pkgs<-c('mgcv','gdata'); 
for (g in pkgs) {  
if (!(g %in% rownames(installed.packages()))) install.packages(g,repos='https://cloud.r-project.org'); 
}
library(mgcv); 
library(gdata); 
WD <- EmpowerStatsR; rm(EmpowerStatsR); gc(); 
title<-'Spline smoothing plot'; 
wd.subset=''; 
weights.var<- NA; 
yvname<-c('CAUDATEHEAD.RIGHT_ROI'); 
ydist<-c('gaussian'); 
ylink<-c('identity'); 
ylv<-c(0); 
par1<-1; 
xvname<-c('AGE','SEX','BMI','EGFR','HT','DM','DL','ETAT.CLIBRE'); 
sxf<-c(0,0,0,0,0,0,0,0); 
xlv<-c(0,2,0,0,2,2,2,4); 
svname<-c('UA'); 
sdf<-c(0); 
slv<-c(0); 
par3<-1; 
timevar<- NA; 
vname.start<- NA; 
vname.stop<- NA; 
subjvname<- NA; 
colvname<- NA; 
chk<- 1; 
dec<-4; 
##R package## mgcv gdata ##R package##;
vec2shift<-function(vnew,vorg,f, opt) {
  if (is.na(f[1])) {
    mean1<-mean(vorg)
    if (opt=="logit") mean1<-log(mean1/(1-mean1))
    if (opt=="log") mean1<-log(mean1)
    vnew<-vnew+(mean1-mean(vnew))
  } else {
    mean1<-tapply(vorg,factor(f),mean)
    if (opt=="logit") mean1<-log(mean1/(1-mean1))
    if (opt=="log") mean1<-log(mean1)
    mean2<-tapply(vnew,factor(f),mean);  meand<-mean1-mean2;  lvf<-levels(factor(f))
    for (z in (1:length(lvf))) {vnew[factor(f)==lvf[z]]<-vnew[factor(f)==lvf[z]]+meand[z]; }
  }
  return(vnew)
}
getNumber<-function(str, n) {
  str<-substr(str,2,nchar(str)-1)
  for (i in (1:nchar(str))) {if (substr(str,i,i)==",") {p=i; break}; }
  ifelse(n==1,return(substr(str,1,p-1)),return(substr(str,p+1,nchar(str))))
}
legLocate<-function(x,y) {
  x[is.infinite(y)]<-NA
  y[is.infinite(y)]<-NA
  xmin<-min(x,na.rm=TRUE); xmax<-max(x,na.rm=TRUE)
  ymin<-min(y,na.rm=TRUE); ymax<-max(y,na.rm=TRUE)
  yoff<-(ymax-ymin); tmp<-table(cut(x,3),cut(y,4))
  tmp.r=which.min(tmp[,4]);tmp.c=4
  if (tmp[2,1]==0) {tmp.r=2;tmp.c=1}
  if (tmp[1,1]==0) {tmp.r=1;tmp.c=1}
  if (tmp[3,1]==0) {tmp.r=3;tmp.c=1}
  if (tmp[2,4]==0) {tmp.r=2;tmp.c=4}
  if (tmp[1,4]==0) {tmp.r=1;tmp.c=4}
  if (tmp[3,4]==0) {tmp.r=3;tmp.c=4}
  pos.y<-colnames(tmp)[tmp.c];  pos.x<-rownames(tmp)[tmp.r];  pct<-0.15
  if (tmp.c==4) {
     if (min(tmp[,4])>0) {pct<-0.3}
     ymax<-ymax+yoff*pct; legy<-ymax;ymin<-ymin-yoff*0.1
  } 
  if (tmp.c==1) {
     if (min(tmp[,1])>0) {pct<-0.3}
     legy<-as.numeric(getNumber(pos.y,2));ymin<-ymin-yoff*pct;ymax=ymax+yoff*0.1
  } 
  legx<-as.numeric(getNumber(pos.x,1))
  return(cbind(xmin,xmax,ymin,ymax,legx,legy))
}
mat2htmltable<-function(mat) {
  t1<- apply(mat,1,function(z) paste(z,collapse="</td><td>"))
  t2<- paste("<tr><td>",t1,"</td></tr>")
  return(paste(t2,collapse=" "))
}
setgam<-function(fml,yi) {
  if (ydist[yi]=="") ydist[yi]<-"gaussian"
  if (ydist[yi]=="exact") ydist[yi]<-"binomial"
  if (ydist[yi]=="breslow") ydist[yi]<-"binomial"
  if (ydist[yi]=="gaussian") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=gaussian(link="identity"))
  if (ydist[yi]=="binomial") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=binomial(link="logit"))
  if (ydist[yi]=="poisson") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=poisson(link="log"))
  if (ydist[yi]=="gamma") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=Gamma(link="inverse"))
  if (ydist[yi]=="negbin") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=negbin(c(1,10), link="log"))
  return(mdl)
}
gam2htmltable<-function(mdl) {
  gs<-summary(mdl)
  np<-length(gs$p.coeff)
  coe<-gs$p.table
  if (gs$family[[2]]=="log" | gs$family[[2]]=="logit") {
    cnames<-c(colnames(coe),"exp(est)","95%CI low","95%CI upp")
    coe<- cbind(coe, exp(coe[,1]), exp(coe[,1]-1.96*coe[,2]), exp(coe[,1]+1.96*coe[,2]))
  }
  if (gs$family[[2]]=="identity") {
    cnames<-c(colnames(coe),"95%CI low","95%CI upp")
    coe<- cbind(coe, coe[,1]-1.96*coe[,2], coe[,1]+1.96*coe[,2])
  }
  oo1<-cbind(c("",rownames(coe)),rbind(cnames,round(coe,dec)))
  oo<-c("</br>Linear terms effect<table border=3>",mat2htmltable(oo1),"</table>")
  if (!is.null(gs$pTerms.table)) {
    xsq<-gs$pTerms.table
    oo2<-cbind(c("",rownames(xsq)),rbind(colnames(xsq),round(xsq,dec)))
    oo<-c(oo, "</br>Chi-square tests for linear terms<table border=3>",mat2htmltable(oo2),"</table>")
  }
  if (!is.null(gs$s.table)) {
   stb<-gs$s.table
   oo3<-cbind(c("",rownames(stb)),rbind(colnames(stb),round(stb,dec)))
   oo<-c(oo, "</br>Approximate significance of smooth terms<table border=3>",mat2htmltable(oo3),"</table>")
  }
  p0<-c("N:", gs$n)
  p1<-c("Adj. r-square:", round(gs$r.sq,4))
  p2<-c("Deviance explained:", round(gs$dev.expl,4))
  p3<-c("UBRE score (sp.criterion):", round(gs$sp.criterion,4))
  p4<-c("Scale estimate:", gs$scale)
  p5<-c("family:", gs$family[[1]])
  p6<-c("link function:", gs$family[[2]])
  oo4<-rbind(p0,p1,p2,p3,p4,p5,p6)
  oo<-c(oo, "</br>Model statistics<table border=3>",mat2htmltable(oo4),"</table>")
  return(oo)
}
gam2pngs<-function(mdl,yi,xi) {
  pred<-predict.gam(mdl,type="terms",se.fit=TRUE)
  mfit<-NA; sfit<-NA; tmp.cname<-NA; kk0<-NA; mfit0<-NA;  ww0<-NULL
  if (xi==0) {kb=1; ke=ns;} else {kb=xi; ke=xi;}
  for (k in (kb:ke)) {
    if (slv[k]==0) {
      mfit<-cbind(mfit,apply(cbind(0,pred$fit[,sxterms[k,]]),1,sum));
      sfit<-cbind(sfit,apply(cbind(0,pred$se.fit[,sxterms[k,]]),1,sum));
      tmp.cname<-c(tmp.cname,svname[k])
      kk0<-c(kk0,k)
    }
  }
  tmp.cname<-tmp.cname[-1];  kk0<-kk0[-1]
  mfit<-matrix(mfit[,-1],ncol=length(tmp.cname)); colnames(mfit)<-paste(tmp.cname,".fit",sep="");
  sfit<-matrix(sfit[,-1],ncol=length(tmp.cname)); colnames(sfit)<-paste(tmp.cname,".se",sep="");
  if (!is.na(colvname)) {tmpfac<-wd[,colvname];} else {tmpfac<-NA;} 
  if (mdl$family[2]=="logit") {
    mfit0<-mfit;  mfit0.low<-mfit0-1.96*sfit;  mfit0.upp<-mfit0+1.96*sfit;
    colnames(mfit0.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit0.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit0)<-paste(tmp.cname,".fit",sep="");
    ww0<-cbind(wd,mfit0); if (is.na(colvname)) ww0<-cbind(ww0,mfit0.low,mfit0.upp)
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac,"logit"))
    mfit.low<-mfit-1.96*sfit; mfit.low<-matrix(exp(mfit.low)/(1+exp(mfit.low)),ncol=length(tmp.cname))
    mfit.upp<-mfit+1.96*sfit; mfit.upp<-matrix(exp(mfit.upp)/(1+exp(mfit.upp)),ncol=length(tmp.cname))
    mfit<-matrix(exp(mfit)/(1+exp(mfit)),ncol=length(tmp.cname))
    colnames(mfit.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit)<-paste(tmp.cname,".fit",sep="");
    ww<-cbind(wd,mfit); if (is.na(colvname)) ww<-cbind(ww,mfit.low,mfit.upp)
  } else if (mdl$family[2]=="log") {
    mfit0<-mfit;  mfit0.low<-mfit0-1.96*sfit;  mfit0.upp<-mfit0+1.96*sfit;
    colnames(mfit0.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit0.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit0)<-paste(tmp.cname,".fit",sep="");
    ww0<-cbind(wd,mfit0); if (is.na(colvname)) ww0<-cbind(ww0,mfit0.low,mfit0.upp)
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac,"log"))
    mfit.low<-mfit-1.96*sfit; mfit.low<-matrix(exp(mfit.low),ncol=length(tmp.cname))
    mfit.upp<-mfit+1.96*sfit; mfit.upp<-matrix(exp(mfit.upp),ncol=length(tmp.cname))
    mfit<-matrix(exp(mfit),ncol=length(tmp.cname))
    colnames(mfit.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit)<-paste(tmp.cname,".fit",sep="");
    ww<-cbind(wd,mfit); if (is.na(colvname)) ww<-cbind(ww,mfit.low,mfit.upp)
  } else if (mdl$family[2]=="identity") {
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac," "))
    ww<-cbind(wd,mfit,sfit)
  } else {
    ww<-cbind(wd,mfit,sfit)
  }
  if (xi!=0) {xf<-paste(ofname,yvname[yi],svname[xi],"gam.xls",sep="_");} else {xf<-paste(ofname,yvname[yi], "gam.xls", sep="_");}
  write.table(ww,file=xf,row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE,quote=FALSE)
  if (!is.null(ww0)) {
    if (xi!=0) {xf<-paste(ofname,yvname[yi],svname[xi],"gam0.xls",sep="_");} else {xf<-paste(ofname,yvname[yi], "gam0.xls", sep="_");}
    write.table(ww0,file=xf,row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE,quote=FALSE)
  }
  px<-c(20,1:9); gg<-"";
  for (k in kk0) {
     cname1<-paste(svname[k],".fit",sep="");  y.tmp<-mfit[,cname1]; y0.tmp<-NULL
     if (mdl$family[2]=="logit" | mdl$family[2]=="log") {
       cname2<-paste(svname[k],".low",sep=""); y.low<-mfit.low[,cname2]
       cname3<-paste(svname[k],".upp",sep=""); y.upp<-mfit.upp[,cname3]
	   y0.tmp<-mfit0[,cname1]; y0.low<-mfit0.low[,cname2]; y0.upp<-mfit0.upp[,cname3]
     } else {
       cname2<-paste(svname[k],".se",sep=""); se.tmp<-sfit[,cname2]; 
       y.low<-y.tmp-1.96*se.tmp; y.upp<-y.tmp+1.96*se.tmp
     }
     x.tmp<-wd[,svname[k]]; 
     pngf<-paste(ofname,yvname[yi],svname[k],"smooth.png",sep="_")
     pdff<-paste(ofname,yvname[yi],svname[k],"smooth.pdf",sep="_")
     pngf0<-paste(ofname,yvname[yi],svname[k],"smooth1.png",sep="_")
     pdff0<-paste(ofname,yvname[yi],svname[k],"smooth1.pdf",sep="_")
     pngf2<-paste(ofname,yvname[yi],svname[k],"smooth2.png",sep="_")
     pdff2<-paste(ofname,yvname[yi],svname[k],"smooth2.pdf",sep="_")
     if (is.na(colvname)) {
       if (chk==1) {tmp.col<-c("red","blue");} else {tmp.col<-rep("black",2);}
       xy<-legLocate(c(x.tmp,x.tmp),c(y.low,y.upp))
       png(pngf,width=720,height=560)
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       dev.off()

       pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       dev.off()

       png(pngf0,width=720,height=560)
       tmp.ord<-order(x.tmp);  x.tmp0<-x.tmp[tmp.ord]; 
       y.tmp0<-y.tmp[tmp.ord]; y.low0<-y.low[tmp.ord]; y.upp0<-y.upp[tmp.ord]
       plot(y.tmp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       rug(x.tmp0)
       dev.off()

       pdf(pdff0,width=pdfwd, height=pdfht, family="GB1");
       tmp.ord<-order(x.tmp);  x.tmp0<-x.tmp[tmp.ord]; 
       y.tmp0<-y.tmp[tmp.ord]; y.low0<-y.low[tmp.ord]; y.upp0<-y.upp[tmp.ord]
       plot(y.tmp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       rug(x.tmp0)
       dev.off()
	   
       if (!is.null(y0.tmp)) {
	     xy0<-legLocate(c(x.tmp,x.tmp),c(y0.low,y0.upp))
		 if (mdl$family[2]=="logit") {yb0<-paste("Log OR of", yb.utf8[yi]);} else {yb0<-paste("Log RR of", yb.utf8[yi]);}
	     png(pngf2,width=720,height=560)
         tmp.ord<-order(x.tmp);    x.tmp0<-x.tmp[tmp.ord]; 
         y0.tmp0<-y0.tmp[tmp.ord]; y0.low0<-y0.low[tmp.ord]; y0.upp0<-y0.upp[tmp.ord]
         plot(y0.tmp0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
         par(new=TRUE); 
         plot(y0.low0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
         par(new=TRUE); 
         plot(y0.upp0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb0, xlab=sb[k])
         rug(x.tmp0)
         dev.off()
	   }

       rm(tmp.ord,y.tmp0,x.tmp0,y.low0,y.upp0)
       xy<-legLocate(x.tmp,wd[,1])
       pngf1<-paste(ofname,yvname[yi],svname[k],"scatter.png",sep="_")
       pdff1<-paste(ofname,yvname[yi],svname[k],"scatter.pdf",sep="_")

       png(pngf1,width=720,height=560)
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1], type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(wd[,1]~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),type="p",pch=1,cex=0.5, ylab=yb.utf8[yi], xlab=sb.utf8[k])    
       dev.off()

       pdf(pdff1,width=pdfwd, height=pdfht, family="GB1");
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1], type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(wd[,1]~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),type="p",pch=1,cex=0.5, ylab=yb.utf8[yi], xlab=sb.utf8[k])    
       dev.off()

     } else {
       if (chk==1) {tmp.col<-rainbow(ncg);tmp.col1<-c("red","blue")} else {tmp.col<-rep("black",ncg);tmp.col1<-c("black","black")}
       for (b in (1:ncg)) {
         y00<-y.tmp[wd[,colvname]==colv.lv[b]]; x00<-x.tmp[wd[,colvname]==colv.lv[b]]
         xy1<-legLocate(x00,y00)
         pngf1<-paste(ofname,yvname[yi],svname[k],colvname,colv.lv[b],"smooth.png",sep="_")
         pdff1<-paste(ofname,yvname[yi],svname[k],colvname,colv.lv[b],"smooth.pdf",sep="_")

         png(pngf1,width=720,height=560)
         plot(y00~x00,ylim=c(xy1[3],xy1[4]),xlim=c(xy1[1],xy1[2]),col=tmp.col1[1],type="p", pch=20, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         dev.off()

         pdf(pdff1,width=pdfwd, height=pdfht, family="GB1");
         plot(y00~x00,ylim=c(xy1[3],xy1[4]),xlim=c(xy1[1],xy1[2]),col=tmp.col1[1],type="p", pch=20, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         dev.off()

       }
       xy<-legLocate(x.tmp,y.tmp)
       png(pngf,width=720,height=560)
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         if (b>1) par(new=TRUE)
         plot(y0~x0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="p", pch=px[b], ylab=yb.utf8[yi], xlab=sb.utf8[k])
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8, pch=px[1:ncg],bty="n",col=tmp.col)
       dev.off()

       pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         if (b>1) par(new=TRUE)
         plot(y0~x0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="p", pch=px[b], ylab=yb[yi], xlab=sb[k])
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8, pch=px[1:ncg],bty="n",col=tmp.col)
       dev.off()

       png(pngf0,width=720,height=560)
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         tmp.ord<-order(x0); x00<-x0[tmp.ord]; y00<-y0[tmp.ord];
         if (b>1) par(new=TRUE)
         plot(y00~x00,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         rm(tmp.ord,x00,y00)
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
       dev.off()

       pdf(pdff0,width=pdfwd, height=pdfht, family="GB1");
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         tmp.ord<-order(x0); x00<-x0[tmp.ord]; y00<-y0[tmp.ord];
         if (b>1) par(new=TRUE)
         plot(y00~x00,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         rm(tmp.ord,x00,y00)
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
       dev.off()

       if (!is.null(y0.tmp)) {
	     xy0<-legLocate(c(x.tmp,x.tmp),c(y0.low,y0.upp))
		 if (mdl$family[2]=="logit") {yb0<-paste("Log OR of", yb.utf8[yi]);} else {yb0<-paste("Log RR of", yb.utf8[yi]);}
	     png(pngf2,width=720,height=560)
         for (b in (1:ncg)) {
           y0<-y0.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
           tmp.ord<-order(x0); x00<-x0[tmp.ord];  y00<-y0[tmp.ord];
           if (b>1) par(new=TRUE)
           plot(y00~x00,ylim=c(xy0[3],xy0[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb0, xlab=sb.utf8[k])
           rm(tmp.ord,x00,y00)
         }
         legend(xy0[5],xy0[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
         dev.off()
	   }

     }
     gg<-c(gg,paste("<td>",yb[yi]," vs. ",sb[k],"</br><a href=\"",pngf,"\" target=_BLANK><img src=\"",pngf,"?time=\" width=320,height=320></a></td>",sep=""))
  }
  return(gg)
}
adjmean<-function(mdl, yi, xi) {
  if (length(xvname)>0) {allvname<- c(xvname,svname);  all.lv<-c(xlv,slv);} else {allvname<-c(svname); all.lv<-slv;}
  if (!is.na(colvname)) allvname<-c(allvname,colvname)
  nv = length(allvname)
  xi.lv <- levels(factor(WD[,svname[xi]]))
  newd0 <- matrix(0,ncol=nv,nrow=length(xi.lv))
  colnames(newd0)<-allvname
  for (b in 1:length(all.lv)) {
    if (all.lv[b]==0) {
      newd0[,b]<-mean(WD[,allvname[b]],na.rm=TRUE)
    } else {
      uniqv <- unique(WD[,allvname[b]])
      newd0[,b]<-uniqv[which.max(tabulate(match(WD[,allvname[b]], uniqv)))]
    }
  }
  newd0[,svname[xi]]<-as.numeric(xi.lv)
  if (!is.na(colvname)) {
    for (b in 1:ncg) {
      newd1<-newd0; newd1[,colvname]<-as.numeric(colv.lv[b]); 
      if (b==1) {newD<-newd1;} else {newD<-rbind(newD,newd1);}
    }
    f<-table(wd[,svname[xi]],wd[,colvname])
  } else {
    newD<-newd0; ncg<-1;
    f<-table(wd[,svname[xi]])
  }
  pred<-predict(mdl, data.frame(newD),se.fit=TRUE)
  meany.pop0<-tapply(wd[,yvname[yi]],wd[,svname[xi]],function(z) mean(z,na.rm=TRUE))[1]
  if (ylink[yi]=="logit") meany.pop0<-log(meany.pop0/(1-meany.pop0));
  if (ylink[yi]=="log") meany.pop0<-log(meany.pop0);
  shift<-meany.pop0-pred$fit[1]
  y.fit <- pred$fit+shift 
  y.low <- pred$fit+shift-pred$se.fit*1.96
  y.upp <- pred$fit+shift+pred$se.fit*1.96
  y.pred<- cbind(y.fit,y.low,y.upp); 
  tmp.ylab<-paste("Mean of", yb.utf8[yi]);
  cname.pred<-c("Mean","Mean.low","Mean.upp")
  if (ylink[yi]=="logit") {
    y.pred<-exp(y.pred); y.pred<-y.pred/(1+y.pred);
    tmp.ylab<-paste("% of", yb.utf8[yi]);
    cname.pred<-c("Rate","Rate.low","Rate.upp")  
  }
  if (ylink[yi]=="log") {
    y.pred<-exp(y.pred); tmp.ylab<-paste("% of", yb.utf8[yi]);
    cname.pred<-c("Rate","Rate.low","Rate.upp")  
  }
  y.fit<-y.pred[,1]; y.low<-y.pred[,2]; y.upp<-y.pred[,3]
  if (length(xvname)>0) {
    tmp.ylab<-paste("Adjusted", tmp.ylab)
    cname.pred<-paste("adj.", cname.pred, sep="")
  }
  y.pred<-round(y.pred,dec)
  colnames(y.pred)<-cname.pred
  y.pred <-cbind(newD[,svname[xi]],y.pred); colnames(y.pred)[1]<-svname[xi]
  if (!is.na(colvname)) {
    y.pred<-cbind(newD[,colvname],y.pred); colnames(y.pred)[1]<-colvname
  }
  y.pred<-rbind(colnames(y.pred),y.pred)
  px<-c(20,1:9)
  if (chk==1) {tmp.col<-rainbow(ncg);tmp.col1<-c("red","blue");
  } else {tmp.col<-rep("black",ncg);tmp.col1<-tmp.col}

  pngf<-paste(ofname,yvname[yi],svname[xi],"adjmean.png",sep="_")
  pdff<-paste(ofname,yvname[yi],svname[xi],"adjmean.pdf",sep="_")

  if (ncg>1) {
    xy<-legLocate(newD[,svname[xi]],c(y.fit))
    png(pngf,width=720,height=560)
    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      if (b>1) par(new=TRUE)
      plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b],type="b", pch=px[b],
        ylab=tmp.ylab, xlab=sb.utf8[xi])
    }
    legend(xy[5],xy[6],colv.lb,title=colvb,pch=px[1:ncg],bty="n",col=tmp.col)
    dev.off()

    pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      if (b>1) par(new=TRUE)
      plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b],type="b", pch=px[b],
        ylab=tmp.ylab, xlab=sb.utf8[xi])
    }
    legend(xy[5],xy[6],colv.lb,title=colvb,pch=px[1:ncg],bty="n",col=tmp.col)
    dev.off()

    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      y.lci<-y.low[newD[,colvname]==colv.lv[b]]
      y.uci<-y.upp[newD[,colvname]==colv.lv[b]]
      xy<-legLocate(c(x.tmp,x.tmp),c(y.lci,y.uci))
      png(paste(ofname,yvname[yi],svname[xi],colvname,colv.lv[b],"CI.png",sep="_"),width=720,height=560)
      plotCI(x.tmp,y=y.tmp,li=y.lci,ui=y.uci,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main=paste(colvb.utf8, colv.lb.utf8[b],sep=": "))
      lines(x.tmp,y.tmp,lty=2)
      dev.off()

      pdf(paste(ofname,yvname[yi],svname[xi],colvname,colv.lv[b],"CI.pdf",sep="_"),width=pdfwd, height=pdfht, family="GB1");
      plotCI(x.tmp,y=y.tmp,li=y.lci,ui=y.uci,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main=paste(colvb.utf8, colv.lb.utf8[b],sep=": "))
      lines(x.tmp,y.tmp,lty=2)
      dev.off()

    }
  } else {
    x.tmp<-newD[,svname[xi]]
    xy<-legLocate(c(x.tmp,x.tmp),c(y.low,y.upp))
    png(pngf,width=720,height=560)
    plotCI(x.tmp,y=y.fit,li=y.low,ui=y.upp,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main="Adjusted mean & 95% CI")
    lines(x.tmp,y.fit,lty=2)
    dev.off()

    pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
    plotCI(x.tmp,y=y.fit,li=y.low,ui=y.upp,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main="Adjusted mean & 95% CI")
    lines(x.tmp,y.fit,lty=2)
    dev.off()
  }
  oo<-c("</br>Adjusted mean ",yb[yi], " by ", sb[xi], "<table border=3>",mat2htmltable(y.pred),"</table>")
  gg<-paste("<td>",yb[yi]," vs. ",sb[xi],"</br><a href=\"",pngf,"\" target=_BLANK><img src=\"",pngf,"?time=\" width=320,height=320></a></td>", sep="")
  return(list(oo,gg))
}
if (!is.na(weights.var)) {weights<-WD[,weights.var];} else {weights<-1;}
WD<-cbind(WD,weights);
vlabelN<-(substr(vlabel,1,1)==" ");
vlabelZ<-vlabel[vlabelN];vlabelV<-vlabel[!vlabelN]
vnameV<-vname[!vlabelN];vnameZ<-vname[vlabelN]
ny<-length(yvname); yb<-vlabelV[match(yvname,vnameV)]; yb[is.na(yb)]<-yvname[is.na(yb)]
ns<-length(svname); sb<-vlabelV[match(svname,vnameV)]; sb[is.na(sb)]<-svname[is.na(sb)]
yb.utf8<-yb; Encoding(yb.utf8)<-"UTF-8"
sb.utf8<-sb; Encoding(sb.utf8)<-"UTF-8"
ssf<-rep(",fx=FALSE", ns); ssf[sdf>0]<-paste(",k=",sdf[sdf>0],sep="")
sxStr<-paste("s(",svname,ssf,sep="")
sxStr[slv>0]<-paste("factor(",svname[slv>0],")",sep="")
sxx<-paste("s(",svname,")",sep="")
sxx[slv>0]<-paste("factor(",svname[slv>0],")",sep="")
sxx<-matrix(sxx,ncol=1)
if (!is.na(colvname)) {
  sxStr[slv==0]<-paste(sxStr[slv==0],",by=factor(", colvname, ")",sep="")
  sxStr[slv>0]<-paste(sxStr[slv>0],"*factor(", colvname, ")",sep="")
  colv.lv<-levels(factor(WD[,colvname])); ncg<-length(colv.lv); colvb<-vlabel[vname==colvname]; 
  colv.lb<-vlabelZ[match(paste(colvname,colv.lv,sep="."),vnameZ)]
  colv.lb[is.na(colv.lb)]<-colv.lv[is.na(colv.lb)]
  colvb<-vlabelV[match(colvname,vnameV)]; if (is.na(colvb)) colvb<-colvname;
  colvb.utf8<-colvb; Encoding(colvb.utf8)<-"UTF-8"; colv.lb.utf8<-colv.lb; Encoding(colv.lb.utf8)<-"UTF-8"
  sxplots<-NA; sxterms<-NA
  for (i in (1:ns)) {
    sxplots<-c(sxplots,paste(svname[i],"_",colvname,colv.lv,sep=""));
    sxterms<-rbind(sxterms,paste(sxx[i,],":factor(",colvname,")",colv.lv,sep=""))
  }
  sxplots<-sxplots[-1]; sxterms<-matrix(sxterms[-1,],ncol=ncg)
  sxterms<-cbind(paste("factor(",colvname,")",sep=""),sxterms)
} else {ncg<-1;sxplots<-svname; sxterms<-sxx;}
sxStr[slv==0]<-paste(sxStr[slv==0],")",sep="")
nx<-0
if (length(xvname)>0) {
  if (!is.na(colvname)) {xvname<-xvname[xvname!=colvname];}
  nx<-length(xvname); 
}
if (nx>0) {
  xb<-vlabelV[match(xvname,vnameV)]; xb[is.na(xb)]<-xvname[is.na(xb)];
  xvv<-xvname; xvv[xlv>2]<-paste("factor(",xvname[xlv>2],")",sep="")
  if (!is.na(colvname)) {
    xvv[sxf=="S" | sxf=="s"]<-paste("factor(",colvname,")*",xvv[sxf=="S" | sxf=="s"],sep="")
  } 
  xv1<-paste(xvv,collapse="+")
}
if (is.na(par1)) par1<-1
if (ny!=ns & par1==2) par1<-1
if (par1==3) {nterms=ns*ncg*15+nx;} else {nterms=ncg*15+nx;}
w<-c("<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'></head><body>")
w<-c(w,paste("<h2>", title, "</h2>"))
wtab<-"</br></br>Generalize additive models</br>"
wpng<-"</br><table>";
for (i in (1:ny)) {
  if (par1!=3) {
    wtmp<-"";
    if (par1==2) {jstart<-i; jstop<-i;} else {jstart<-1; jstop<-ns;}
    for (j in (jstart:jstop)) {
      tmp.xx<-c(yvname[i],svname[j])
      if (nx>0) tmp.xx<-c(tmp.xx,xvname)
      if (!is.na(colvname[1])) tmp.xx<-c(tmp.xx,colvname)
      wd<-WD[,tmp.xx];
      wd<-wd[apply(is.na(wd),1,sum)==0,]
      fml<-paste(yvname[i],"~",sxStr[j],sep="")
      if (!is.na(colvname)) fml<-paste(fml,"+factor(",colvname,")",sep="")
      if (nx>0) fml<-paste(fml,"+",xv1,sep="")
      tmp.gam<-setgam(fml,i)
      wtab<-c(wtab,paste("</br></br>Outcome:",yb[i]))
      wtab<-c(wtab,paste("</br>Exposure:",sb[j]))
      wtab<-c(wtab,gam2htmltable(tmp.gam))
      if (slv[j]==0) {
        wtmp<-c(wtmp,gam2pngs(tmp.gam,i,j))
      } else {
        stmp<-adjmean(tmp.gam,i,j)
        wtmp<-c(wtmp,stmp[[2]])
        wtab<-c(wtab,stmp[[1]])
      }
    }
    wpng<-c(wpng,"<tr>",wtmp,"</tr>")
  } else {
      tmp.xx<-c(yvname[i],svname);
      if (nx>0) tmp.xx<-c(tmp.xx,xvname)
      if (!is.na(colvname[1])) tmp.xx<-c(tmp.xx,colvname)
      wd<-WD[,tmp.xx];
      wd<-wd[apply(is.na(wd),1,sum)==0,]
      fml<-paste(yvname[i],"~",paste(sxStr,collapse="+"),sep="")
      if (!is.na(colvname)) fml<-paste(fml,"+factor(",colvname,")",sep="")
      if (nx>0) fml<-paste(fml,"+",xv1,sep="")
      tmp.gam<-setgam(fml,i)
      wtab<-c(wtab,paste("</br></br>Outcome:",yb[i]))
      wtab<-c(wtab,gam2htmltable(tmp.gam))
      if (sum(slv==0)>0) {   
        wpng<-c(wpng,"<tr>",gam2pngs(tmp.gam,i,0),"</tr>")
      }
      if (sum(slv>0)>0) {   
        for (k in 1:ns) {
          if (slv[k]>0) {
            stmp<-adjmean(tmp.gam,i,k)
            wtab<-c(wtab,stmp[[1]])
            wpng<-c(wpng,stmp[[2]])
          }
        }
      }
  }
}
wpng<-c(wpng,"</table>")
w<-c(w,wpng,wtab)
w<-c(w,wd.subset)
w<-c(w,paste("</br>Created by EmpowerStats (www.empowerstats.com) and R on",Sys.Date()))
w<-c(w,"</body></html>")
fileConn<-file(paste(ofname,".htm",sep="")); writeLines(w, fileConn)

# Smooth-curve-fitting model：Putamen_Left
Sys.setlocale(category = 'LC_ALL', locale = 'English_United States.1252'); 
.libPaths(file.path(R.home(),'library')); 
library(doBy); 
options(timeout=600); 
library(plotrix); 
library(stringi); 
library(stringr); 
library(survival); 
library(rms); 
library(nnet); 
library(car); 
library(mgcv); 
pdfwd<-6; pdfht<-6; 
load('C:/EmpowerXYS/A123/Analysis/data/UAstandardfibnegative.Rdata'); 
if (length(which(ls()=='EmpowerStatsR'))==0) EmpowerStatsR<-get(ls()[1]); 
names(EmpowerStatsR)<-toupper(names(EmpowerStatsR)); 
originalVNAME<-names(EmpowerStatsR); 
ofname<-'data_3_tbl'; 
vname<-c(NA,'AGE','SEX','SEX.0','SEX.1','PUTAMEN.LEFT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','ETAT.CLIBRE.0','ETAT.CLIBRE.1','ETAT.CLIBRE.2','ETAT.CLIBRE.3','HT','HT.0','HT.1','DM','DM.0','DM.1','DL','DL.0','DL.1')[-1]; 
vlabel<-c(NA,'AGE','SEX','  0','  1','PUTAMEN.LEFT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','  0','  1','  2','  3','HT','  0','  1','DM','  0','  1','DL','  0','  1')[-1]; 
varused4this <- c('AGE','SEX','PUTAMEN.LEFT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','HT','DM','DL'); 
pkgs<-c('mgcv','gdata'); 
for (g in pkgs) {  
if (!(g %in% rownames(installed.packages()))) install.packages(g,repos='https://cloud.r-project.org'); 
}
library(mgcv); 
library(gdata); 
WD <- EmpowerStatsR; rm(EmpowerStatsR); gc(); 
title<-'Spline smoothing plot'; 
wd.subset=''; 
weights.var<- NA; 
yvname<-c('PUTAMEN.LEFT_ROI'); 
ydist<-c('gaussian'); 
ylink<-c('identity'); 
ylv<-c(0); 
par1<-1; 
xvname<-c('AGE','SEX','BMI','EGFR','HT','DM','DL','ETAT.CLIBRE'); 
sxf<-c(0,0,0,0,0,0,0,0); 
xlv<-c(0,2,0,0,2,2,2,4); 
svname<-c('UA'); 
sdf<-c(0); 
slv<-c(0); 
par3<-1; 
timevar<- NA; 
vname.start<- NA; 
vname.stop<- NA; 
subjvname<- NA; 
colvname<- NA; 
chk<- 1; 
dec<-4; 
##R package## mgcv gdata ##R package##;
vec2shift<-function(vnew,vorg,f, opt) {
  if (is.na(f[1])) {
    mean1<-mean(vorg)
    if (opt=="logit") mean1<-log(mean1/(1-mean1))
    if (opt=="log") mean1<-log(mean1)
    vnew<-vnew+(mean1-mean(vnew))
  } else {
    mean1<-tapply(vorg,factor(f),mean)
    if (opt=="logit") mean1<-log(mean1/(1-mean1))
    if (opt=="log") mean1<-log(mean1)
    mean2<-tapply(vnew,factor(f),mean);  meand<-mean1-mean2;  lvf<-levels(factor(f))
    for (z in (1:length(lvf))) {vnew[factor(f)==lvf[z]]<-vnew[factor(f)==lvf[z]]+meand[z]; }
  }
  return(vnew)
}
getNumber<-function(str, n) {
  str<-substr(str,2,nchar(str)-1)
  for (i in (1:nchar(str))) {if (substr(str,i,i)==",") {p=i; break}; }
  ifelse(n==1,return(substr(str,1,p-1)),return(substr(str,p+1,nchar(str))))
}
legLocate<-function(x,y) {
  x[is.infinite(y)]<-NA
  y[is.infinite(y)]<-NA
  xmin<-min(x,na.rm=TRUE); xmax<-max(x,na.rm=TRUE)
  ymin<-min(y,na.rm=TRUE); ymax<-max(y,na.rm=TRUE)
  yoff<-(ymax-ymin); tmp<-table(cut(x,3),cut(y,4))
  tmp.r=which.min(tmp[,4]);tmp.c=4
  if (tmp[2,1]==0) {tmp.r=2;tmp.c=1}
  if (tmp[1,1]==0) {tmp.r=1;tmp.c=1}
  if (tmp[3,1]==0) {tmp.r=3;tmp.c=1}
  if (tmp[2,4]==0) {tmp.r=2;tmp.c=4}
  if (tmp[1,4]==0) {tmp.r=1;tmp.c=4}
  if (tmp[3,4]==0) {tmp.r=3;tmp.c=4}
  pos.y<-colnames(tmp)[tmp.c];  pos.x<-rownames(tmp)[tmp.r];  pct<-0.15
  if (tmp.c==4) {
     if (min(tmp[,4])>0) {pct<-0.3}
     ymax<-ymax+yoff*pct; legy<-ymax;ymin<-ymin-yoff*0.1
  } 
  if (tmp.c==1) {
     if (min(tmp[,1])>0) {pct<-0.3}
     legy<-as.numeric(getNumber(pos.y,2));ymin<-ymin-yoff*pct;ymax=ymax+yoff*0.1
  } 
  legx<-as.numeric(getNumber(pos.x,1))
  return(cbind(xmin,xmax,ymin,ymax,legx,legy))
}
mat2htmltable<-function(mat) {
  t1<- apply(mat,1,function(z) paste(z,collapse="</td><td>"))
  t2<- paste("<tr><td>",t1,"</td></tr>")
  return(paste(t2,collapse=" "))
}
setgam<-function(fml,yi) {
  if (ydist[yi]=="") ydist[yi]<-"gaussian"
  if (ydist[yi]=="exact") ydist[yi]<-"binomial"
  if (ydist[yi]=="breslow") ydist[yi]<-"binomial"
  if (ydist[yi]=="gaussian") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=gaussian(link="identity"))
  if (ydist[yi]=="binomial") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=binomial(link="logit"))
  if (ydist[yi]=="poisson") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=poisson(link="log"))
  if (ydist[yi]=="gamma") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=Gamma(link="inverse"))
  if (ydist[yi]=="negbin") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=negbin(c(1,10), link="log"))
  return(mdl)
}
gam2htmltable<-function(mdl) {
  gs<-summary(mdl)
  np<-length(gs$p.coeff)
  coe<-gs$p.table
  if (gs$family[[2]]=="log" | gs$family[[2]]=="logit") {
    cnames<-c(colnames(coe),"exp(est)","95%CI low","95%CI upp")
    coe<- cbind(coe, exp(coe[,1]), exp(coe[,1]-1.96*coe[,2]), exp(coe[,1]+1.96*coe[,2]))
  }
  if (gs$family[[2]]=="identity") {
    cnames<-c(colnames(coe),"95%CI low","95%CI upp")
    coe<- cbind(coe, coe[,1]-1.96*coe[,2], coe[,1]+1.96*coe[,2])
  }
  oo1<-cbind(c("",rownames(coe)),rbind(cnames,round(coe,dec)))
  oo<-c("</br>Linear terms effect<table border=3>",mat2htmltable(oo1),"</table>")
  if (!is.null(gs$pTerms.table)) {
    xsq<-gs$pTerms.table
    oo2<-cbind(c("",rownames(xsq)),rbind(colnames(xsq),round(xsq,dec)))
    oo<-c(oo, "</br>Chi-square tests for linear terms<table border=3>",mat2htmltable(oo2),"</table>")
  }
  if (!is.null(gs$s.table)) {
   stb<-gs$s.table
   oo3<-cbind(c("",rownames(stb)),rbind(colnames(stb),round(stb,dec)))
   oo<-c(oo, "</br>Approximate significance of smooth terms<table border=3>",mat2htmltable(oo3),"</table>")
  }
  p0<-c("N:", gs$n)
  p1<-c("Adj. r-square:", round(gs$r.sq,4))
  p2<-c("Deviance explained:", round(gs$dev.expl,4))
  p3<-c("UBRE score (sp.criterion):", round(gs$sp.criterion,4))
  p4<-c("Scale estimate:", gs$scale)
  p5<-c("family:", gs$family[[1]])
  p6<-c("link function:", gs$family[[2]])
  oo4<-rbind(p0,p1,p2,p3,p4,p5,p6)
  oo<-c(oo, "</br>Model statistics<table border=3>",mat2htmltable(oo4),"</table>")
  return(oo)
}
gam2pngs<-function(mdl,yi,xi) {
  pred<-predict.gam(mdl,type="terms",se.fit=TRUE)
  mfit<-NA; sfit<-NA; tmp.cname<-NA; kk0<-NA; mfit0<-NA;  ww0<-NULL
  if (xi==0) {kb=1; ke=ns;} else {kb=xi; ke=xi;}
  for (k in (kb:ke)) {
    if (slv[k]==0) {
      mfit<-cbind(mfit,apply(cbind(0,pred$fit[,sxterms[k,]]),1,sum));
      sfit<-cbind(sfit,apply(cbind(0,pred$se.fit[,sxterms[k,]]),1,sum));
      tmp.cname<-c(tmp.cname,svname[k])
      kk0<-c(kk0,k)
    }
  }
  tmp.cname<-tmp.cname[-1];  kk0<-kk0[-1]
  mfit<-matrix(mfit[,-1],ncol=length(tmp.cname)); colnames(mfit)<-paste(tmp.cname,".fit",sep="");
  sfit<-matrix(sfit[,-1],ncol=length(tmp.cname)); colnames(sfit)<-paste(tmp.cname,".se",sep="");
  if (!is.na(colvname)) {tmpfac<-wd[,colvname];} else {tmpfac<-NA;} 
  if (mdl$family[2]=="logit") {
    mfit0<-mfit;  mfit0.low<-mfit0-1.96*sfit;  mfit0.upp<-mfit0+1.96*sfit;
    colnames(mfit0.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit0.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit0)<-paste(tmp.cname,".fit",sep="");
    ww0<-cbind(wd,mfit0); if (is.na(colvname)) ww0<-cbind(ww0,mfit0.low,mfit0.upp)
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac,"logit"))
    mfit.low<-mfit-1.96*sfit; mfit.low<-matrix(exp(mfit.low)/(1+exp(mfit.low)),ncol=length(tmp.cname))
    mfit.upp<-mfit+1.96*sfit; mfit.upp<-matrix(exp(mfit.upp)/(1+exp(mfit.upp)),ncol=length(tmp.cname))
    mfit<-matrix(exp(mfit)/(1+exp(mfit)),ncol=length(tmp.cname))
    colnames(mfit.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit)<-paste(tmp.cname,".fit",sep="");
    ww<-cbind(wd,mfit); if (is.na(colvname)) ww<-cbind(ww,mfit.low,mfit.upp)
  } else if (mdl$family[2]=="log") {
    mfit0<-mfit;  mfit0.low<-mfit0-1.96*sfit;  mfit0.upp<-mfit0+1.96*sfit;
    colnames(mfit0.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit0.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit0)<-paste(tmp.cname,".fit",sep="");
    ww0<-cbind(wd,mfit0); if (is.na(colvname)) ww0<-cbind(ww0,mfit0.low,mfit0.upp)
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac,"log"))
    mfit.low<-mfit-1.96*sfit; mfit.low<-matrix(exp(mfit.low),ncol=length(tmp.cname))
    mfit.upp<-mfit+1.96*sfit; mfit.upp<-matrix(exp(mfit.upp),ncol=length(tmp.cname))
    mfit<-matrix(exp(mfit),ncol=length(tmp.cname))
    colnames(mfit.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit)<-paste(tmp.cname,".fit",sep="");
    ww<-cbind(wd,mfit); if (is.na(colvname)) ww<-cbind(ww,mfit.low,mfit.upp)
  } else if (mdl$family[2]=="identity") {
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac," "))
    ww<-cbind(wd,mfit,sfit)
  } else {
    ww<-cbind(wd,mfit,sfit)
  }
  if (xi!=0) {xf<-paste(ofname,yvname[yi],svname[xi],"gam.xls",sep="_");} else {xf<-paste(ofname,yvname[yi], "gam.xls", sep="_");}
  write.table(ww,file=xf,row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE,quote=FALSE)
  if (!is.null(ww0)) {
    if (xi!=0) {xf<-paste(ofname,yvname[yi],svname[xi],"gam0.xls",sep="_");} else {xf<-paste(ofname,yvname[yi], "gam0.xls", sep="_");}
    write.table(ww0,file=xf,row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE,quote=FALSE)
  }
  px<-c(20,1:9); gg<-"";
  for (k in kk0) {
     cname1<-paste(svname[k],".fit",sep="");  y.tmp<-mfit[,cname1]; y0.tmp<-NULL
     if (mdl$family[2]=="logit" | mdl$family[2]=="log") {
       cname2<-paste(svname[k],".low",sep=""); y.low<-mfit.low[,cname2]
       cname3<-paste(svname[k],".upp",sep=""); y.upp<-mfit.upp[,cname3]
	   y0.tmp<-mfit0[,cname1]; y0.low<-mfit0.low[,cname2]; y0.upp<-mfit0.upp[,cname3]
     } else {
       cname2<-paste(svname[k],".se",sep=""); se.tmp<-sfit[,cname2]; 
       y.low<-y.tmp-1.96*se.tmp; y.upp<-y.tmp+1.96*se.tmp
     }
     x.tmp<-wd[,svname[k]]; 
     pngf<-paste(ofname,yvname[yi],svname[k],"smooth.png",sep="_")
     pdff<-paste(ofname,yvname[yi],svname[k],"smooth.pdf",sep="_")
     pngf0<-paste(ofname,yvname[yi],svname[k],"smooth1.png",sep="_")
     pdff0<-paste(ofname,yvname[yi],svname[k],"smooth1.pdf",sep="_")
     pngf2<-paste(ofname,yvname[yi],svname[k],"smooth2.png",sep="_")
     pdff2<-paste(ofname,yvname[yi],svname[k],"smooth2.pdf",sep="_")
     if (is.na(colvname)) {
       if (chk==1) {tmp.col<-c("red","blue");} else {tmp.col<-rep("black",2);}
       xy<-legLocate(c(x.tmp,x.tmp),c(y.low,y.upp))
       png(pngf,width=720,height=560)
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       dev.off()

       pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       dev.off()

       png(pngf0,width=720,height=560)
       tmp.ord<-order(x.tmp);  x.tmp0<-x.tmp[tmp.ord]; 
       y.tmp0<-y.tmp[tmp.ord]; y.low0<-y.low[tmp.ord]; y.upp0<-y.upp[tmp.ord]
       plot(y.tmp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       rug(x.tmp0)
       dev.off()

       pdf(pdff0,width=pdfwd, height=pdfht, family="GB1");
       tmp.ord<-order(x.tmp);  x.tmp0<-x.tmp[tmp.ord]; 
       y.tmp0<-y.tmp[tmp.ord]; y.low0<-y.low[tmp.ord]; y.upp0<-y.upp[tmp.ord]
       plot(y.tmp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       rug(x.tmp0)
       dev.off()
	   
       if (!is.null(y0.tmp)) {
	     xy0<-legLocate(c(x.tmp,x.tmp),c(y0.low,y0.upp))
		 if (mdl$family[2]=="logit") {yb0<-paste("Log OR of", yb.utf8[yi]);} else {yb0<-paste("Log RR of", yb.utf8[yi]);}
	     png(pngf2,width=720,height=560)
         tmp.ord<-order(x.tmp);    x.tmp0<-x.tmp[tmp.ord]; 
         y0.tmp0<-y0.tmp[tmp.ord]; y0.low0<-y0.low[tmp.ord]; y0.upp0<-y0.upp[tmp.ord]
         plot(y0.tmp0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
         par(new=TRUE); 
         plot(y0.low0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
         par(new=TRUE); 
         plot(y0.upp0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb0, xlab=sb[k])
         rug(x.tmp0)
         dev.off()
	   }

       rm(tmp.ord,y.tmp0,x.tmp0,y.low0,y.upp0)
       xy<-legLocate(x.tmp,wd[,1])
       pngf1<-paste(ofname,yvname[yi],svname[k],"scatter.png",sep="_")
       pdff1<-paste(ofname,yvname[yi],svname[k],"scatter.pdf",sep="_")

       png(pngf1,width=720,height=560)
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1], type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(wd[,1]~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),type="p",pch=1,cex=0.5, ylab=yb.utf8[yi], xlab=sb.utf8[k])    
       dev.off()

       pdf(pdff1,width=pdfwd, height=pdfht, family="GB1");
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1], type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(wd[,1]~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),type="p",pch=1,cex=0.5, ylab=yb.utf8[yi], xlab=sb.utf8[k])    
       dev.off()

     } else {
       if (chk==1) {tmp.col<-rainbow(ncg);tmp.col1<-c("red","blue")} else {tmp.col<-rep("black",ncg);tmp.col1<-c("black","black")}
       for (b in (1:ncg)) {
         y00<-y.tmp[wd[,colvname]==colv.lv[b]]; x00<-x.tmp[wd[,colvname]==colv.lv[b]]
         xy1<-legLocate(x00,y00)
         pngf1<-paste(ofname,yvname[yi],svname[k],colvname,colv.lv[b],"smooth.png",sep="_")
         pdff1<-paste(ofname,yvname[yi],svname[k],colvname,colv.lv[b],"smooth.pdf",sep="_")

         png(pngf1,width=720,height=560)
         plot(y00~x00,ylim=c(xy1[3],xy1[4]),xlim=c(xy1[1],xy1[2]),col=tmp.col1[1],type="p", pch=20, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         dev.off()

         pdf(pdff1,width=pdfwd, height=pdfht, family="GB1");
         plot(y00~x00,ylim=c(xy1[3],xy1[4]),xlim=c(xy1[1],xy1[2]),col=tmp.col1[1],type="p", pch=20, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         dev.off()

       }
       xy<-legLocate(x.tmp,y.tmp)
       png(pngf,width=720,height=560)
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         if (b>1) par(new=TRUE)
         plot(y0~x0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="p", pch=px[b], ylab=yb.utf8[yi], xlab=sb.utf8[k])
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8, pch=px[1:ncg],bty="n",col=tmp.col)
       dev.off()

       pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         if (b>1) par(new=TRUE)
         plot(y0~x0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="p", pch=px[b], ylab=yb[yi], xlab=sb[k])
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8, pch=px[1:ncg],bty="n",col=tmp.col)
       dev.off()

       png(pngf0,width=720,height=560)
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         tmp.ord<-order(x0); x00<-x0[tmp.ord]; y00<-y0[tmp.ord];
         if (b>1) par(new=TRUE)
         plot(y00~x00,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         rm(tmp.ord,x00,y00)
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
       dev.off()

       pdf(pdff0,width=pdfwd, height=pdfht, family="GB1");
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         tmp.ord<-order(x0); x00<-x0[tmp.ord]; y00<-y0[tmp.ord];
         if (b>1) par(new=TRUE)
         plot(y00~x00,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         rm(tmp.ord,x00,y00)
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
       dev.off()

       if (!is.null(y0.tmp)) {
	     xy0<-legLocate(c(x.tmp,x.tmp),c(y0.low,y0.upp))
		 if (mdl$family[2]=="logit") {yb0<-paste("Log OR of", yb.utf8[yi]);} else {yb0<-paste("Log RR of", yb.utf8[yi]);}
	     png(pngf2,width=720,height=560)
         for (b in (1:ncg)) {
           y0<-y0.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
           tmp.ord<-order(x0); x00<-x0[tmp.ord];  y00<-y0[tmp.ord];
           if (b>1) par(new=TRUE)
           plot(y00~x00,ylim=c(xy0[3],xy0[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb0, xlab=sb.utf8[k])
           rm(tmp.ord,x00,y00)
         }
         legend(xy0[5],xy0[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
         dev.off()
	   }

     }
     gg<-c(gg,paste("<td>",yb[yi]," vs. ",sb[k],"</br><a href=\"",pngf,"\" target=_BLANK><img src=\"",pngf,"?time=\" width=320,height=320></a></td>",sep=""))
  }
  return(gg)
}
adjmean<-function(mdl, yi, xi) {
  if (length(xvname)>0) {allvname<- c(xvname,svname);  all.lv<-c(xlv,slv);} else {allvname<-c(svname); all.lv<-slv;}
  if (!is.na(colvname)) allvname<-c(allvname,colvname)
  nv = length(allvname)
  xi.lv <- levels(factor(WD[,svname[xi]]))
  newd0 <- matrix(0,ncol=nv,nrow=length(xi.lv))
  colnames(newd0)<-allvname
  for (b in 1:length(all.lv)) {
    if (all.lv[b]==0) {
      newd0[,b]<-mean(WD[,allvname[b]],na.rm=TRUE)
    } else {
      uniqv <- unique(WD[,allvname[b]])
      newd0[,b]<-uniqv[which.max(tabulate(match(WD[,allvname[b]], uniqv)))]
    }
  }
  newd0[,svname[xi]]<-as.numeric(xi.lv)
  if (!is.na(colvname)) {
    for (b in 1:ncg) {
      newd1<-newd0; newd1[,colvname]<-as.numeric(colv.lv[b]); 
      if (b==1) {newD<-newd1;} else {newD<-rbind(newD,newd1);}
    }
    f<-table(wd[,svname[xi]],wd[,colvname])
  } else {
    newD<-newd0; ncg<-1;
    f<-table(wd[,svname[xi]])
  }
  pred<-predict(mdl, data.frame(newD),se.fit=TRUE)
  meany.pop0<-tapply(wd[,yvname[yi]],wd[,svname[xi]],function(z) mean(z,na.rm=TRUE))[1]
  if (ylink[yi]=="logit") meany.pop0<-log(meany.pop0/(1-meany.pop0));
  if (ylink[yi]=="log") meany.pop0<-log(meany.pop0);
  shift<-meany.pop0-pred$fit[1]
  y.fit <- pred$fit+shift 
  y.low <- pred$fit+shift-pred$se.fit*1.96
  y.upp <- pred$fit+shift+pred$se.fit*1.96
  y.pred<- cbind(y.fit,y.low,y.upp); 
  tmp.ylab<-paste("Mean of", yb.utf8[yi]);
  cname.pred<-c("Mean","Mean.low","Mean.upp")
  if (ylink[yi]=="logit") {
    y.pred<-exp(y.pred); y.pred<-y.pred/(1+y.pred);
    tmp.ylab<-paste("% of", yb.utf8[yi]);
    cname.pred<-c("Rate","Rate.low","Rate.upp")  
  }
  if (ylink[yi]=="log") {
    y.pred<-exp(y.pred); tmp.ylab<-paste("% of", yb.utf8[yi]);
    cname.pred<-c("Rate","Rate.low","Rate.upp")  
  }
  y.fit<-y.pred[,1]; y.low<-y.pred[,2]; y.upp<-y.pred[,3]
  if (length(xvname)>0) {
    tmp.ylab<-paste("Adjusted", tmp.ylab)
    cname.pred<-paste("adj.", cname.pred, sep="")
  }
  y.pred<-round(y.pred,dec)
  colnames(y.pred)<-cname.pred
  y.pred <-cbind(newD[,svname[xi]],y.pred); colnames(y.pred)[1]<-svname[xi]
  if (!is.na(colvname)) {
    y.pred<-cbind(newD[,colvname],y.pred); colnames(y.pred)[1]<-colvname
  }
  y.pred<-rbind(colnames(y.pred),y.pred)
  px<-c(20,1:9)
  if (chk==1) {tmp.col<-rainbow(ncg);tmp.col1<-c("red","blue");
  } else {tmp.col<-rep("black",ncg);tmp.col1<-tmp.col}

  pngf<-paste(ofname,yvname[yi],svname[xi],"adjmean.png",sep="_")
  pdff<-paste(ofname,yvname[yi],svname[xi],"adjmean.pdf",sep="_")

  if (ncg>1) {
    xy<-legLocate(newD[,svname[xi]],c(y.fit))
    png(pngf,width=720,height=560)
    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      if (b>1) par(new=TRUE)
      plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b],type="b", pch=px[b],
        ylab=tmp.ylab, xlab=sb.utf8[xi])
    }
    legend(xy[5],xy[6],colv.lb,title=colvb,pch=px[1:ncg],bty="n",col=tmp.col)
    dev.off()

    pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      if (b>1) par(new=TRUE)
      plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b],type="b", pch=px[b],
        ylab=tmp.ylab, xlab=sb.utf8[xi])
    }
    legend(xy[5],xy[6],colv.lb,title=colvb,pch=px[1:ncg],bty="n",col=tmp.col)
    dev.off()

    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      y.lci<-y.low[newD[,colvname]==colv.lv[b]]
      y.uci<-y.upp[newD[,colvname]==colv.lv[b]]
      xy<-legLocate(c(x.tmp,x.tmp),c(y.lci,y.uci))
      png(paste(ofname,yvname[yi],svname[xi],colvname,colv.lv[b],"CI.png",sep="_"),width=720,height=560)
      plotCI(x.tmp,y=y.tmp,li=y.lci,ui=y.uci,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main=paste(colvb.utf8, colv.lb.utf8[b],sep=": "))
      lines(x.tmp,y.tmp,lty=2)
      dev.off()

      pdf(paste(ofname,yvname[yi],svname[xi],colvname,colv.lv[b],"CI.pdf",sep="_"),width=pdfwd, height=pdfht, family="GB1");
      plotCI(x.tmp,y=y.tmp,li=y.lci,ui=y.uci,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main=paste(colvb.utf8, colv.lb.utf8[b],sep=": "))
      lines(x.tmp,y.tmp,lty=2)
      dev.off()

    }
  } else {
    x.tmp<-newD[,svname[xi]]
    xy<-legLocate(c(x.tmp,x.tmp),c(y.low,y.upp))
    png(pngf,width=720,height=560)
    plotCI(x.tmp,y=y.fit,li=y.low,ui=y.upp,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main="Adjusted mean & 95% CI")
    lines(x.tmp,y.fit,lty=2)
    dev.off()

    pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
    plotCI(x.tmp,y=y.fit,li=y.low,ui=y.upp,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main="Adjusted mean & 95% CI")
    lines(x.tmp,y.fit,lty=2)
    dev.off()
  }
  oo<-c("</br>Adjusted mean ",yb[yi], " by ", sb[xi], "<table border=3>",mat2htmltable(y.pred),"</table>")
  gg<-paste("<td>",yb[yi]," vs. ",sb[xi],"</br><a href=\"",pngf,"\" target=_BLANK><img src=\"",pngf,"?time=\" width=320,height=320></a></td>", sep="")
  return(list(oo,gg))
}
if (!is.na(weights.var)) {weights<-WD[,weights.var];} else {weights<-1;}
WD<-cbind(WD,weights);
vlabelN<-(substr(vlabel,1,1)==" ");
vlabelZ<-vlabel[vlabelN];vlabelV<-vlabel[!vlabelN]
vnameV<-vname[!vlabelN];vnameZ<-vname[vlabelN]
ny<-length(yvname); yb<-vlabelV[match(yvname,vnameV)]; yb[is.na(yb)]<-yvname[is.na(yb)]
ns<-length(svname); sb<-vlabelV[match(svname,vnameV)]; sb[is.na(sb)]<-svname[is.na(sb)]
yb.utf8<-yb; Encoding(yb.utf8)<-"UTF-8"
sb.utf8<-sb; Encoding(sb.utf8)<-"UTF-8"
ssf<-rep(",fx=FALSE", ns); ssf[sdf>0]<-paste(",k=",sdf[sdf>0],sep="")
sxStr<-paste("s(",svname,ssf,sep="")
sxStr[slv>0]<-paste("factor(",svname[slv>0],")",sep="")
sxx<-paste("s(",svname,")",sep="")
sxx[slv>0]<-paste("factor(",svname[slv>0],")",sep="")
sxx<-matrix(sxx,ncol=1)
if (!is.na(colvname)) {
  sxStr[slv==0]<-paste(sxStr[slv==0],",by=factor(", colvname, ")",sep="")
  sxStr[slv>0]<-paste(sxStr[slv>0],"*factor(", colvname, ")",sep="")
  colv.lv<-levels(factor(WD[,colvname])); ncg<-length(colv.lv); colvb<-vlabel[vname==colvname]; 
  colv.lb<-vlabelZ[match(paste(colvname,colv.lv,sep="."),vnameZ)]
  colv.lb[is.na(colv.lb)]<-colv.lv[is.na(colv.lb)]
  colvb<-vlabelV[match(colvname,vnameV)]; if (is.na(colvb)) colvb<-colvname;
  colvb.utf8<-colvb; Encoding(colvb.utf8)<-"UTF-8"; colv.lb.utf8<-colv.lb; Encoding(colv.lb.utf8)<-"UTF-8"
  sxplots<-NA; sxterms<-NA
  for (i in (1:ns)) {
    sxplots<-c(sxplots,paste(svname[i],"_",colvname,colv.lv,sep=""));
    sxterms<-rbind(sxterms,paste(sxx[i,],":factor(",colvname,")",colv.lv,sep=""))
  }
  sxplots<-sxplots[-1]; sxterms<-matrix(sxterms[-1,],ncol=ncg)
  sxterms<-cbind(paste("factor(",colvname,")",sep=""),sxterms)
} else {ncg<-1;sxplots<-svname; sxterms<-sxx;}
sxStr[slv==0]<-paste(sxStr[slv==0],")",sep="")
nx<-0
if (length(xvname)>0) {
  if (!is.na(colvname)) {xvname<-xvname[xvname!=colvname];}
  nx<-length(xvname); 
}
if (nx>0) {
  xb<-vlabelV[match(xvname,vnameV)]; xb[is.na(xb)]<-xvname[is.na(xb)];
  xvv<-xvname; xvv[xlv>2]<-paste("factor(",xvname[xlv>2],")",sep="")
  if (!is.na(colvname)) {
    xvv[sxf=="S" | sxf=="s"]<-paste("factor(",colvname,")*",xvv[sxf=="S" | sxf=="s"],sep="")
  } 
  xv1<-paste(xvv,collapse="+")
}
if (is.na(par1)) par1<-1
if (ny!=ns & par1==2) par1<-1
if (par1==3) {nterms=ns*ncg*15+nx;} else {nterms=ncg*15+nx;}
w<-c("<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'></head><body>")
w<-c(w,paste("<h2>", title, "</h2>"))
wtab<-"</br></br>Generalize additive models</br>"
wpng<-"</br><table>";
for (i in (1:ny)) {
  if (par1!=3) {
    wtmp<-"";
    if (par1==2) {jstart<-i; jstop<-i;} else {jstart<-1; jstop<-ns;}
    for (j in (jstart:jstop)) {
      tmp.xx<-c(yvname[i],svname[j])
      if (nx>0) tmp.xx<-c(tmp.xx,xvname)
      if (!is.na(colvname[1])) tmp.xx<-c(tmp.xx,colvname)
      wd<-WD[,tmp.xx];
      wd<-wd[apply(is.na(wd),1,sum)==0,]
      fml<-paste(yvname[i],"~",sxStr[j],sep="")
      if (!is.na(colvname)) fml<-paste(fml,"+factor(",colvname,")",sep="")
      if (nx>0) fml<-paste(fml,"+",xv1,sep="")
      tmp.gam<-setgam(fml,i)
      wtab<-c(wtab,paste("</br></br>Outcome:",yb[i]))
      wtab<-c(wtab,paste("</br>Exposure:",sb[j]))
      wtab<-c(wtab,gam2htmltable(tmp.gam))
      if (slv[j]==0) {
        wtmp<-c(wtmp,gam2pngs(tmp.gam,i,j))
      } else {
        stmp<-adjmean(tmp.gam,i,j)
        wtmp<-c(wtmp,stmp[[2]])
        wtab<-c(wtab,stmp[[1]])
      }
    }
    wpng<-c(wpng,"<tr>",wtmp,"</tr>")
  } else {
      tmp.xx<-c(yvname[i],svname);
      if (nx>0) tmp.xx<-c(tmp.xx,xvname)
      if (!is.na(colvname[1])) tmp.xx<-c(tmp.xx,colvname)
      wd<-WD[,tmp.xx];
      wd<-wd[apply(is.na(wd),1,sum)==0,]
      fml<-paste(yvname[i],"~",paste(sxStr,collapse="+"),sep="")
      if (!is.na(colvname)) fml<-paste(fml,"+factor(",colvname,")",sep="")
      if (nx>0) fml<-paste(fml,"+",xv1,sep="")
      tmp.gam<-setgam(fml,i)
      wtab<-c(wtab,paste("</br></br>Outcome:",yb[i]))
      wtab<-c(wtab,gam2htmltable(tmp.gam))
      if (sum(slv==0)>0) {   
        wpng<-c(wpng,"<tr>",gam2pngs(tmp.gam,i,0),"</tr>")
      }
      if (sum(slv>0)>0) {   
        for (k in 1:ns) {
          if (slv[k]>0) {
            stmp<-adjmean(tmp.gam,i,k)
            wtab<-c(wtab,stmp[[1]])
            wpng<-c(wpng,stmp[[2]])
          }
        }
      }
  }
}
wpng<-c(wpng,"</table>")
w<-c(w,wpng,wtab)
w<-c(w,wd.subset)
w<-c(w,paste("</br>Created by EmpowerStats (www.empowerstats.com) and R on",Sys.Date()))
w<-c(w,"</body></html>")
fileConn<-file(paste(ofname,".htm",sep="")); writeLines(w, fileConn)

# Smooth-curve-fitting model：Putamen_Right
Sys.setlocale(category = 'LC_ALL', locale = 'English_United States.1252'); 
.libPaths(file.path(R.home(),'library')); 
library(doBy); 
options(timeout=600); 
library(plotrix); 
library(stringi); 
library(stringr); 
library(survival); 
library(rms); 
library(nnet); 
library(car); 
library(mgcv); 
pdfwd<-6; pdfht<-6; 
load('C:/EmpowerXYS/A123/Analysis/data/UAstandardfibnegative.Rdata'); 
if (length(which(ls()=='EmpowerStatsR'))==0) EmpowerStatsR<-get(ls()[1]); 
names(EmpowerStatsR)<-toupper(names(EmpowerStatsR)); 
originalVNAME<-names(EmpowerStatsR); 
ofname<-'data_4_tbl'; 
vname<-c(NA,'AGE','SEX','SEX.0','SEX.1','PUTAMEN.RIGHT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','ETAT.CLIBRE.0','ETAT.CLIBRE.1','ETAT.CLIBRE.2','ETAT.CLIBRE.3','HT','HT.0','HT.1','DM','DM.0','DM.1','DL','DL.0','DL.1')[-1]; 
vlabel<-c(NA,'AGE','SEX','  0','  1','PUTAMEN.RIGHT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','  0','  1','  2','  3','HT','  0','  1','DM','  0','  1','DL','  0','  1')[-1]; 
varused4this <- c('AGE','SEX','PUTAMEN.RIGHT_ROI','BMI','EGFR','UA','ETAT.CLIBRE','HT','DM','DL'); 
pkgs<-c('mgcv','gdata'); 
for (g in pkgs) {  
if (!(g %in% rownames(installed.packages()))) install.packages(g,repos='https://cloud.r-project.org'); 
}
library(mgcv); 
library(gdata); 
WD <- EmpowerStatsR; rm(EmpowerStatsR); gc(); 
title<-'Spline smoothing plot'; 
wd.subset=''; 
weights.var<- NA; 
yvname<-c('PUTAMEN.RIGHT_ROI'); 
ydist<-c('gaussian'); 
ylink<-c('identity'); 
ylv<-c(0); 
par1<-1; 
xvname<-c('AGE','SEX','BMI','EGFR','HT','DM','DL','ETAT.CLIBRE'); 
sxf<-c(0,0,0,0,0,0,0,0); 
xlv<-c(0,2,0,0,2,2,2,4); 
svname<-c('UA'); 
sdf<-c(0); 
slv<-c(0); 
par3<-1; 
timevar<- NA; 
vname.start<- NA; 
vname.stop<- NA; 
subjvname<- NA; 
colvname<- NA; 
chk<- 1; 
dec<-4; 
##R package## mgcv gdata ##R package##;
vec2shift<-function(vnew,vorg,f, opt) {
  if (is.na(f[1])) {
    mean1<-mean(vorg)
    if (opt=="logit") mean1<-log(mean1/(1-mean1))
    if (opt=="log") mean1<-log(mean1)
    vnew<-vnew+(mean1-mean(vnew))
  } else {
    mean1<-tapply(vorg,factor(f),mean)
    if (opt=="logit") mean1<-log(mean1/(1-mean1))
    if (opt=="log") mean1<-log(mean1)
    mean2<-tapply(vnew,factor(f),mean);  meand<-mean1-mean2;  lvf<-levels(factor(f))
    for (z in (1:length(lvf))) {vnew[factor(f)==lvf[z]]<-vnew[factor(f)==lvf[z]]+meand[z]; }
  }
  return(vnew)
}
getNumber<-function(str, n) {
  str<-substr(str,2,nchar(str)-1)
  for (i in (1:nchar(str))) {if (substr(str,i,i)==",") {p=i; break}; }
  ifelse(n==1,return(substr(str,1,p-1)),return(substr(str,p+1,nchar(str))))
}
legLocate<-function(x,y) {
  x[is.infinite(y)]<-NA
  y[is.infinite(y)]<-NA
  xmin<-min(x,na.rm=TRUE); xmax<-max(x,na.rm=TRUE)
  ymin<-min(y,na.rm=TRUE); ymax<-max(y,na.rm=TRUE)
  yoff<-(ymax-ymin); tmp<-table(cut(x,3),cut(y,4))
  tmp.r=which.min(tmp[,4]);tmp.c=4
  if (tmp[2,1]==0) {tmp.r=2;tmp.c=1}
  if (tmp[1,1]==0) {tmp.r=1;tmp.c=1}
  if (tmp[3,1]==0) {tmp.r=3;tmp.c=1}
  if (tmp[2,4]==0) {tmp.r=2;tmp.c=4}
  if (tmp[1,4]==0) {tmp.r=1;tmp.c=4}
  if (tmp[3,4]==0) {tmp.r=3;tmp.c=4}
  pos.y<-colnames(tmp)[tmp.c];  pos.x<-rownames(tmp)[tmp.r];  pct<-0.15
  if (tmp.c==4) {
     if (min(tmp[,4])>0) {pct<-0.3}
     ymax<-ymax+yoff*pct; legy<-ymax;ymin<-ymin-yoff*0.1
  } 
  if (tmp.c==1) {
     if (min(tmp[,1])>0) {pct<-0.3}
     legy<-as.numeric(getNumber(pos.y,2));ymin<-ymin-yoff*pct;ymax=ymax+yoff*0.1
  } 
  legx<-as.numeric(getNumber(pos.x,1))
  return(cbind(xmin,xmax,ymin,ymax,legx,legy))
}
mat2htmltable<-function(mat) {
  t1<- apply(mat,1,function(z) paste(z,collapse="</td><td>"))
  t2<- paste("<tr><td>",t1,"</td></tr>")
  return(paste(t2,collapse=" "))
}
setgam<-function(fml,yi) {
  if (ydist[yi]=="") ydist[yi]<-"gaussian"
  if (ydist[yi]=="exact") ydist[yi]<-"binomial"
  if (ydist[yi]=="breslow") ydist[yi]<-"binomial"
  if (ydist[yi]=="gaussian") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=gaussian(link="identity"))
  if (ydist[yi]=="binomial") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=binomial(link="logit"))
  if (ydist[yi]=="poisson") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=poisson(link="log"))
  if (ydist[yi]=="gamma") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=Gamma(link="inverse"))
  if (ydist[yi]=="negbin") mdl<-gam(formula(fml),weights=wd$weights,data=wd, family=negbin(c(1,10), link="log"))
  return(mdl)
}
gam2htmltable<-function(mdl) {
  gs<-summary(mdl)
  np<-length(gs$p.coeff)
  coe<-gs$p.table
  if (gs$family[[2]]=="log" | gs$family[[2]]=="logit") {
    cnames<-c(colnames(coe),"exp(est)","95%CI low","95%CI upp")
    coe<- cbind(coe, exp(coe[,1]), exp(coe[,1]-1.96*coe[,2]), exp(coe[,1]+1.96*coe[,2]))
  }
  if (gs$family[[2]]=="identity") {
    cnames<-c(colnames(coe),"95%CI low","95%CI upp")
    coe<- cbind(coe, coe[,1]-1.96*coe[,2], coe[,1]+1.96*coe[,2])
  }
  oo1<-cbind(c("",rownames(coe)),rbind(cnames,round(coe,dec)))
  oo<-c("</br>Linear terms effect<table border=3>",mat2htmltable(oo1),"</table>")
  if (!is.null(gs$pTerms.table)) {
    xsq<-gs$pTerms.table
    oo2<-cbind(c("",rownames(xsq)),rbind(colnames(xsq),round(xsq,dec)))
    oo<-c(oo, "</br>Chi-square tests for linear terms<table border=3>",mat2htmltable(oo2),"</table>")
  }
  if (!is.null(gs$s.table)) {
   stb<-gs$s.table
   oo3<-cbind(c("",rownames(stb)),rbind(colnames(stb),round(stb,dec)))
   oo<-c(oo, "</br>Approximate significance of smooth terms<table border=3>",mat2htmltable(oo3),"</table>")
  }
  p0<-c("N:", gs$n)
  p1<-c("Adj. r-square:", round(gs$r.sq,4))
  p2<-c("Deviance explained:", round(gs$dev.expl,4))
  p3<-c("UBRE score (sp.criterion):", round(gs$sp.criterion,4))
  p4<-c("Scale estimate:", gs$scale)
  p5<-c("family:", gs$family[[1]])
  p6<-c("link function:", gs$family[[2]])
  oo4<-rbind(p0,p1,p2,p3,p4,p5,p6)
  oo<-c(oo, "</br>Model statistics<table border=3>",mat2htmltable(oo4),"</table>")
  return(oo)
}
gam2pngs<-function(mdl,yi,xi) {
  pred<-predict.gam(mdl,type="terms",se.fit=TRUE)
  mfit<-NA; sfit<-NA; tmp.cname<-NA; kk0<-NA; mfit0<-NA;  ww0<-NULL
  if (xi==0) {kb=1; ke=ns;} else {kb=xi; ke=xi;}
  for (k in (kb:ke)) {
    if (slv[k]==0) {
      mfit<-cbind(mfit,apply(cbind(0,pred$fit[,sxterms[k,]]),1,sum));
      sfit<-cbind(sfit,apply(cbind(0,pred$se.fit[,sxterms[k,]]),1,sum));
      tmp.cname<-c(tmp.cname,svname[k])
      kk0<-c(kk0,k)
    }
  }
  tmp.cname<-tmp.cname[-1];  kk0<-kk0[-1]
  mfit<-matrix(mfit[,-1],ncol=length(tmp.cname)); colnames(mfit)<-paste(tmp.cname,".fit",sep="");
  sfit<-matrix(sfit[,-1],ncol=length(tmp.cname)); colnames(sfit)<-paste(tmp.cname,".se",sep="");
  if (!is.na(colvname)) {tmpfac<-wd[,colvname];} else {tmpfac<-NA;} 
  if (mdl$family[2]=="logit") {
    mfit0<-mfit;  mfit0.low<-mfit0-1.96*sfit;  mfit0.upp<-mfit0+1.96*sfit;
    colnames(mfit0.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit0.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit0)<-paste(tmp.cname,".fit",sep="");
    ww0<-cbind(wd,mfit0); if (is.na(colvname)) ww0<-cbind(ww0,mfit0.low,mfit0.upp)
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac,"logit"))
    mfit.low<-mfit-1.96*sfit; mfit.low<-matrix(exp(mfit.low)/(1+exp(mfit.low)),ncol=length(tmp.cname))
    mfit.upp<-mfit+1.96*sfit; mfit.upp<-matrix(exp(mfit.upp)/(1+exp(mfit.upp)),ncol=length(tmp.cname))
    mfit<-matrix(exp(mfit)/(1+exp(mfit)),ncol=length(tmp.cname))
    colnames(mfit.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit)<-paste(tmp.cname,".fit",sep="");
    ww<-cbind(wd,mfit); if (is.na(colvname)) ww<-cbind(ww,mfit.low,mfit.upp)
  } else if (mdl$family[2]=="log") {
    mfit0<-mfit;  mfit0.low<-mfit0-1.96*sfit;  mfit0.upp<-mfit0+1.96*sfit;
    colnames(mfit0.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit0.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit0)<-paste(tmp.cname,".fit",sep="");
    ww0<-cbind(wd,mfit0); if (is.na(colvname)) ww0<-cbind(ww0,mfit0.low,mfit0.upp)
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac,"log"))
    mfit.low<-mfit-1.96*sfit; mfit.low<-matrix(exp(mfit.low),ncol=length(tmp.cname))
    mfit.upp<-mfit+1.96*sfit; mfit.upp<-matrix(exp(mfit.upp),ncol=length(tmp.cname))
    mfit<-matrix(exp(mfit),ncol=length(tmp.cname))
    colnames(mfit.low)<-paste(tmp.cname,".low",sep="");
    colnames(mfit.upp)<-paste(tmp.cname,".upp",sep="");
    colnames(mfit)<-paste(tmp.cname,".fit",sep="");
    ww<-cbind(wd,mfit); if (is.na(colvname)) ww<-cbind(ww,mfit.low,mfit.upp)
  } else if (mdl$family[2]=="identity") {
    mfit<-apply(mfit,2,function(z) vec2shift(z,mdl$fitted.value,tmpfac," "))
    ww<-cbind(wd,mfit,sfit)
  } else {
    ww<-cbind(wd,mfit,sfit)
  }
  if (xi!=0) {xf<-paste(ofname,yvname[yi],svname[xi],"gam.xls",sep="_");} else {xf<-paste(ofname,yvname[yi], "gam.xls", sep="_");}
  write.table(ww,file=xf,row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE,quote=FALSE)
  if (!is.null(ww0)) {
    if (xi!=0) {xf<-paste(ofname,yvname[yi],svname[xi],"gam0.xls",sep="_");} else {xf<-paste(ofname,yvname[yi], "gam0.xls", sep="_");}
    write.table(ww0,file=xf,row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE,quote=FALSE)
  }
  px<-c(20,1:9); gg<-"";
  for (k in kk0) {
     cname1<-paste(svname[k],".fit",sep="");  y.tmp<-mfit[,cname1]; y0.tmp<-NULL
     if (mdl$family[2]=="logit" | mdl$family[2]=="log") {
       cname2<-paste(svname[k],".low",sep=""); y.low<-mfit.low[,cname2]
       cname3<-paste(svname[k],".upp",sep=""); y.upp<-mfit.upp[,cname3]
	   y0.tmp<-mfit0[,cname1]; y0.low<-mfit0.low[,cname2]; y0.upp<-mfit0.upp[,cname3]
     } else {
       cname2<-paste(svname[k],".se",sep=""); se.tmp<-sfit[,cname2]; 
       y.low<-y.tmp-1.96*se.tmp; y.upp<-y.tmp+1.96*se.tmp
     }
     x.tmp<-wd[,svname[k]]; 
     pngf<-paste(ofname,yvname[yi],svname[k],"smooth.png",sep="_")
     pdff<-paste(ofname,yvname[yi],svname[k],"smooth.pdf",sep="_")
     pngf0<-paste(ofname,yvname[yi],svname[k],"smooth1.png",sep="_")
     pdff0<-paste(ofname,yvname[yi],svname[k],"smooth1.pdf",sep="_")
     pngf2<-paste(ofname,yvname[yi],svname[k],"smooth2.png",sep="_")
     pdff2<-paste(ofname,yvname[yi],svname[k],"smooth2.pdf",sep="_")
     if (is.na(colvname)) {
       if (chk==1) {tmp.col<-c("red","blue");} else {tmp.col<-rep("black",2);}
       xy<-legLocate(c(x.tmp,x.tmp),c(y.low,y.upp))
       png(pngf,width=720,height=560)
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       dev.off()

       pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="p", pch=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       dev.off()

       png(pngf0,width=720,height=560)
       tmp.ord<-order(x.tmp);  x.tmp0<-x.tmp[tmp.ord]; 
       y.tmp0<-y.tmp[tmp.ord]; y.low0<-y.low[tmp.ord]; y.upp0<-y.upp[tmp.ord]
       plot(y.tmp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       rug(x.tmp0)
       dev.off()

       pdf(pdff0,width=pdfwd, height=pdfht, family="GB1");
       tmp.ord<-order(x.tmp);  x.tmp0<-x.tmp[tmp.ord]; 
       y.tmp0<-y.tmp[tmp.ord]; y.low0<-y.low[tmp.ord]; y.upp0<-y.upp[tmp.ord]
       plot(y.tmp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.low0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
       par(new=TRUE); 
       plot(y.upp0~x.tmp0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb.utf8[yi], xlab=sb.utf8[k])
       rug(x.tmp0)
       dev.off()
	   
       if (!is.null(y0.tmp)) {
	     xy0<-legLocate(c(x.tmp,x.tmp),c(y0.low,y0.upp))
		 if (mdl$family[2]=="logit") {yb0<-paste("Log OR of", yb.utf8[yi]);} else {yb0<-paste("Log RR of", yb.utf8[yi]);}
	     png(pngf2,width=720,height=560)
         tmp.ord<-order(x.tmp);    x.tmp0<-x.tmp[tmp.ord]; 
         y0.tmp0<-y0.tmp[tmp.ord]; y0.low0<-y0.low[tmp.ord]; y0.upp0<-y0.upp[tmp.ord]
         plot(y0.tmp0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[1],type="l", lty=1, lwd=2, ylab="", xlab="")
         par(new=TRUE); 
         plot(y0.low0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab="", xlab="")
         par(new=TRUE); 
         plot(y0.upp0~x.tmp0,ylim=c(xy0[3],xy0[4]),xlim=c(xy0[1],xy0[2]),col=tmp.col[2], type="l", lty=3, lwd=1, ylab=yb0, xlab=sb[k])
         rug(x.tmp0)
         dev.off()
	   }

       rm(tmp.ord,y.tmp0,x.tmp0,y.low0,y.upp0)
       xy<-legLocate(x.tmp,wd[,1])
       pngf1<-paste(ofname,yvname[yi],svname[k],"scatter.png",sep="_")
       pdff1<-paste(ofname,yvname[yi],svname[k],"scatter.pdf",sep="_")

       png(pngf1,width=720,height=560)
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1], type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(wd[,1]~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),type="p",pch=1,cex=0.5, ylab=yb.utf8[yi], xlab=sb.utf8[k])    
       dev.off()

       pdf(pdff1,width=pdfwd, height=pdfht, family="GB1");
       plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[1], type="p", pch=20, ylab="", xlab="")
       par(new=TRUE); 
       plot(wd[,1]~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),type="p",pch=1,cex=0.5, ylab=yb.utf8[yi], xlab=sb.utf8[k])    
       dev.off()

     } else {
       if (chk==1) {tmp.col<-rainbow(ncg);tmp.col1<-c("red","blue")} else {tmp.col<-rep("black",ncg);tmp.col1<-c("black","black")}
       for (b in (1:ncg)) {
         y00<-y.tmp[wd[,colvname]==colv.lv[b]]; x00<-x.tmp[wd[,colvname]==colv.lv[b]]
         xy1<-legLocate(x00,y00)
         pngf1<-paste(ofname,yvname[yi],svname[k],colvname,colv.lv[b],"smooth.png",sep="_")
         pdff1<-paste(ofname,yvname[yi],svname[k],colvname,colv.lv[b],"smooth.pdf",sep="_")

         png(pngf1,width=720,height=560)
         plot(y00~x00,ylim=c(xy1[3],xy1[4]),xlim=c(xy1[1],xy1[2]),col=tmp.col1[1],type="p", pch=20, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         dev.off()

         pdf(pdff1,width=pdfwd, height=pdfht, family="GB1");
         plot(y00~x00,ylim=c(xy1[3],xy1[4]),xlim=c(xy1[1],xy1[2]),col=tmp.col1[1],type="p", pch=20, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         dev.off()

       }
       xy<-legLocate(x.tmp,y.tmp)
       png(pngf,width=720,height=560)
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         if (b>1) par(new=TRUE)
         plot(y0~x0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="p", pch=px[b], ylab=yb.utf8[yi], xlab=sb.utf8[k])
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8, pch=px[1:ncg],bty="n",col=tmp.col)
       dev.off()

       pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         if (b>1) par(new=TRUE)
         plot(y0~x0,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="p", pch=px[b], ylab=yb[yi], xlab=sb[k])
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8, pch=px[1:ncg],bty="n",col=tmp.col)
       dev.off()

       png(pngf0,width=720,height=560)
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         tmp.ord<-order(x0); x00<-x0[tmp.ord]; y00<-y0[tmp.ord];
         if (b>1) par(new=TRUE)
         plot(y00~x00,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         rm(tmp.ord,x00,y00)
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
       dev.off()

       pdf(pdff0,width=pdfwd, height=pdfht, family="GB1");
       for (b in (1:ncg)) {
         y0<-y.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
         tmp.ord<-order(x0); x00<-x0[tmp.ord]; y00<-y0[tmp.ord];
         if (b>1) par(new=TRUE)
         plot(y00~x00,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb.utf8[yi], xlab=sb.utf8[k])
         rm(tmp.ord,x00,y00)
       }
       legend(xy[5],xy[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
       dev.off()

       if (!is.null(y0.tmp)) {
	     xy0<-legLocate(c(x.tmp,x.tmp),c(y0.low,y0.upp))
		 if (mdl$family[2]=="logit") {yb0<-paste("Log OR of", yb.utf8[yi]);} else {yb0<-paste("Log RR of", yb.utf8[yi]);}
	     png(pngf2,width=720,height=560)
         for (b in (1:ncg)) {
           y0<-y0.tmp[wd[,colvname]==colv.lv[b]]; x0<-x.tmp[wd[,colvname]==colv.lv[b]]
           tmp.ord<-order(x0); x00<-x0[tmp.ord];  y00<-y0[tmp.ord];
           if (b>1) par(new=TRUE)
           plot(y00~x00,ylim=c(xy0[3],xy0[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b], type="l", lty=b, lwd=2, ylab=yb0, xlab=sb.utf8[k])
           rm(tmp.ord,x00,y00)
         }
         legend(xy0[5],xy0[6],colv.lb.utf8,title=colvb.utf8,lty=(1:ncg),bty="n",col=tmp.col)
         dev.off()
	   }

     }
     gg<-c(gg,paste("<td>",yb[yi]," vs. ",sb[k],"</br><a href=\"",pngf,"\" target=_BLANK><img src=\"",pngf,"?time=\" width=320,height=320></a></td>",sep=""))
  }
  return(gg)
}
adjmean<-function(mdl, yi, xi) {
  if (length(xvname)>0) {allvname<- c(xvname,svname);  all.lv<-c(xlv,slv);} else {allvname<-c(svname); all.lv<-slv;}
  if (!is.na(colvname)) allvname<-c(allvname,colvname)
  nv = length(allvname)
  xi.lv <- levels(factor(WD[,svname[xi]]))
  newd0 <- matrix(0,ncol=nv,nrow=length(xi.lv))
  colnames(newd0)<-allvname
  for (b in 1:length(all.lv)) {
    if (all.lv[b]==0) {
      newd0[,b]<-mean(WD[,allvname[b]],na.rm=TRUE)
    } else {
      uniqv <- unique(WD[,allvname[b]])
      newd0[,b]<-uniqv[which.max(tabulate(match(WD[,allvname[b]], uniqv)))]
    }
  }
  newd0[,svname[xi]]<-as.numeric(xi.lv)
  if (!is.na(colvname)) {
    for (b in 1:ncg) {
      newd1<-newd0; newd1[,colvname]<-as.numeric(colv.lv[b]); 
      if (b==1) {newD<-newd1;} else {newD<-rbind(newD,newd1);}
    }
    f<-table(wd[,svname[xi]],wd[,colvname])
  } else {
    newD<-newd0; ncg<-1;
    f<-table(wd[,svname[xi]])
  }
  pred<-predict(mdl, data.frame(newD),se.fit=TRUE)
  meany.pop0<-tapply(wd[,yvname[yi]],wd[,svname[xi]],function(z) mean(z,na.rm=TRUE))[1]
  if (ylink[yi]=="logit") meany.pop0<-log(meany.pop0/(1-meany.pop0));
  if (ylink[yi]=="log") meany.pop0<-log(meany.pop0);
  shift<-meany.pop0-pred$fit[1]
  y.fit <- pred$fit+shift 
  y.low <- pred$fit+shift-pred$se.fit*1.96
  y.upp <- pred$fit+shift+pred$se.fit*1.96
  y.pred<- cbind(y.fit,y.low,y.upp); 
  tmp.ylab<-paste("Mean of", yb.utf8[yi]);
  cname.pred<-c("Mean","Mean.low","Mean.upp")
  if (ylink[yi]=="logit") {
    y.pred<-exp(y.pred); y.pred<-y.pred/(1+y.pred);
    tmp.ylab<-paste("% of", yb.utf8[yi]);
    cname.pred<-c("Rate","Rate.low","Rate.upp")  
  }
  if (ylink[yi]=="log") {
    y.pred<-exp(y.pred); tmp.ylab<-paste("% of", yb.utf8[yi]);
    cname.pred<-c("Rate","Rate.low","Rate.upp")  
  }
  y.fit<-y.pred[,1]; y.low<-y.pred[,2]; y.upp<-y.pred[,3]
  if (length(xvname)>0) {
    tmp.ylab<-paste("Adjusted", tmp.ylab)
    cname.pred<-paste("adj.", cname.pred, sep="")
  }
  y.pred<-round(y.pred,dec)
  colnames(y.pred)<-cname.pred
  y.pred <-cbind(newD[,svname[xi]],y.pred); colnames(y.pred)[1]<-svname[xi]
  if (!is.na(colvname)) {
    y.pred<-cbind(newD[,colvname],y.pred); colnames(y.pred)[1]<-colvname
  }
  y.pred<-rbind(colnames(y.pred),y.pred)
  px<-c(20,1:9)
  if (chk==1) {tmp.col<-rainbow(ncg);tmp.col1<-c("red","blue");
  } else {tmp.col<-rep("black",ncg);tmp.col1<-tmp.col}

  pngf<-paste(ofname,yvname[yi],svname[xi],"adjmean.png",sep="_")
  pdff<-paste(ofname,yvname[yi],svname[xi],"adjmean.pdf",sep="_")

  if (ncg>1) {
    xy<-legLocate(newD[,svname[xi]],c(y.fit))
    png(pngf,width=720,height=560)
    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      if (b>1) par(new=TRUE)
      plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b],type="b", pch=px[b],
        ylab=tmp.ylab, xlab=sb.utf8[xi])
    }
    legend(xy[5],xy[6],colv.lb,title=colvb,pch=px[1:ncg],bty="n",col=tmp.col)
    dev.off()

    pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      if (b>1) par(new=TRUE)
      plot(y.tmp~x.tmp,ylim=c(xy[3],xy[4]),xlim=c(xy[1],xy[2]),col=tmp.col[b],type="b", pch=px[b],
        ylab=tmp.ylab, xlab=sb.utf8[xi])
    }
    legend(xy[5],xy[6],colv.lb,title=colvb,pch=px[1:ncg],bty="n",col=tmp.col)
    dev.off()

    for (b in 1:ncg) {
      x.tmp<-newD[newD[,colvname]==colv.lv[b],svname[xi]]
      y.tmp<-y.fit[newD[,colvname]==colv.lv[b]]
      y.lci<-y.low[newD[,colvname]==colv.lv[b]]
      y.uci<-y.upp[newD[,colvname]==colv.lv[b]]
      xy<-legLocate(c(x.tmp,x.tmp),c(y.lci,y.uci))
      png(paste(ofname,yvname[yi],svname[xi],colvname,colv.lv[b],"CI.png",sep="_"),width=720,height=560)
      plotCI(x.tmp,y=y.tmp,li=y.lci,ui=y.uci,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main=paste(colvb.utf8, colv.lb.utf8[b],sep=": "))
      lines(x.tmp,y.tmp,lty=2)
      dev.off()

      pdf(paste(ofname,yvname[yi],svname[xi],colvname,colv.lv[b],"CI.pdf",sep="_"),width=pdfwd, height=pdfht, family="GB1");
      plotCI(x.tmp,y=y.tmp,li=y.lci,ui=y.uci,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main=paste(colvb.utf8, colv.lb.utf8[b],sep=": "))
      lines(x.tmp,y.tmp,lty=2)
      dev.off()

    }
  } else {
    x.tmp<-newD[,svname[xi]]
    xy<-legLocate(c(x.tmp,x.tmp),c(y.low,y.upp))
    png(pngf,width=720,height=560)
    plotCI(x.tmp,y=y.fit,li=y.low,ui=y.upp,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main="Adjusted mean & 95% CI")
    lines(x.tmp,y.fit,lty=2)
    dev.off()

    pdf(pdff,width=pdfwd, height=pdfht, family="GB1");
    plotCI(x.tmp,y=y.fit,li=y.low,ui=y.upp,pch=20,lwd=1,col=tmp.col[1], xlim=c(xy[1],xy[2]),ylim=c(xy[3],xy[4]),
         ylab=tmp.ylab, xlab=sb.utf8[xi], main="Adjusted mean & 95% CI")
    lines(x.tmp,y.fit,lty=2)
    dev.off()
  }
  oo<-c("</br>Adjusted mean ",yb[yi], " by ", sb[xi], "<table border=3>",mat2htmltable(y.pred),"</table>")
  gg<-paste("<td>",yb[yi]," vs. ",sb[xi],"</br><a href=\"",pngf,"\" target=_BLANK><img src=\"",pngf,"?time=\" width=320,height=320></a></td>", sep="")
  return(list(oo,gg))
}
if (!is.na(weights.var)) {weights<-WD[,weights.var];} else {weights<-1;}
WD<-cbind(WD,weights);
vlabelN<-(substr(vlabel,1,1)==" ");
vlabelZ<-vlabel[vlabelN];vlabelV<-vlabel[!vlabelN]
vnameV<-vname[!vlabelN];vnameZ<-vname[vlabelN]
ny<-length(yvname); yb<-vlabelV[match(yvname,vnameV)]; yb[is.na(yb)]<-yvname[is.na(yb)]
ns<-length(svname); sb<-vlabelV[match(svname,vnameV)]; sb[is.na(sb)]<-svname[is.na(sb)]
yb.utf8<-yb; Encoding(yb.utf8)<-"UTF-8"
sb.utf8<-sb; Encoding(sb.utf8)<-"UTF-8"
ssf<-rep(",fx=FALSE", ns); ssf[sdf>0]<-paste(",k=",sdf[sdf>0],sep="")
sxStr<-paste("s(",svname,ssf,sep="")
sxStr[slv>0]<-paste("factor(",svname[slv>0],")",sep="")
sxx<-paste("s(",svname,")",sep="")
sxx[slv>0]<-paste("factor(",svname[slv>0],")",sep="")
sxx<-matrix(sxx,ncol=1)
if (!is.na(colvname)) {
  sxStr[slv==0]<-paste(sxStr[slv==0],",by=factor(", colvname, ")",sep="")
  sxStr[slv>0]<-paste(sxStr[slv>0],"*factor(", colvname, ")",sep="")
  colv.lv<-levels(factor(WD[,colvname])); ncg<-length(colv.lv); colvb<-vlabel[vname==colvname]; 
  colv.lb<-vlabelZ[match(paste(colvname,colv.lv,sep="."),vnameZ)]
  colv.lb[is.na(colv.lb)]<-colv.lv[is.na(colv.lb)]
  colvb<-vlabelV[match(colvname,vnameV)]; if (is.na(colvb)) colvb<-colvname;
  colvb.utf8<-colvb; Encoding(colvb.utf8)<-"UTF-8"; colv.lb.utf8<-colv.lb; Encoding(colv.lb.utf8)<-"UTF-8"
  sxplots<-NA; sxterms<-NA
  for (i in (1:ns)) {
    sxplots<-c(sxplots,paste(svname[i],"_",colvname,colv.lv,sep=""));
    sxterms<-rbind(sxterms,paste(sxx[i,],":factor(",colvname,")",colv.lv,sep=""))
  }
  sxplots<-sxplots[-1]; sxterms<-matrix(sxterms[-1,],ncol=ncg)
  sxterms<-cbind(paste("factor(",colvname,")",sep=""),sxterms)
} else {ncg<-1;sxplots<-svname; sxterms<-sxx;}
sxStr[slv==0]<-paste(sxStr[slv==0],")",sep="")
nx<-0
if (length(xvname)>0) {
  if (!is.na(colvname)) {xvname<-xvname[xvname!=colvname];}
  nx<-length(xvname); 
}
if (nx>0) {
  xb<-vlabelV[match(xvname,vnameV)]; xb[is.na(xb)]<-xvname[is.na(xb)];
  xvv<-xvname; xvv[xlv>2]<-paste("factor(",xvname[xlv>2],")",sep="")
  if (!is.na(colvname)) {
    xvv[sxf=="S" | sxf=="s"]<-paste("factor(",colvname,")*",xvv[sxf=="S" | sxf=="s"],sep="")
  } 
  xv1<-paste(xvv,collapse="+")
}
if (is.na(par1)) par1<-1
if (ny!=ns & par1==2) par1<-1
if (par1==3) {nterms=ns*ncg*15+nx;} else {nterms=ncg*15+nx;}
w<-c("<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'></head><body>")
w<-c(w,paste("<h2>", title, "</h2>"))
wtab<-"</br></br>Generalize additive models</br>"
wpng<-"</br><table>";
for (i in (1:ny)) {
  if (par1!=3) {
    wtmp<-"";
    if (par1==2) {jstart<-i; jstop<-i;} else {jstart<-1; jstop<-ns;}
    for (j in (jstart:jstop)) {
      tmp.xx<-c(yvname[i],svname[j])
      if (nx>0) tmp.xx<-c(tmp.xx,xvname)
      if (!is.na(colvname[1])) tmp.xx<-c(tmp.xx,colvname)
      wd<-WD[,tmp.xx];
      wd<-wd[apply(is.na(wd),1,sum)==0,]
      fml<-paste(yvname[i],"~",sxStr[j],sep="")
      if (!is.na(colvname)) fml<-paste(fml,"+factor(",colvname,")",sep="")
      if (nx>0) fml<-paste(fml,"+",xv1,sep="")
      tmp.gam<-setgam(fml,i)
      wtab<-c(wtab,paste("</br></br>Outcome:",yb[i]))
      wtab<-c(wtab,paste("</br>Exposure:",sb[j]))
      wtab<-c(wtab,gam2htmltable(tmp.gam))
      if (slv[j]==0) {
        wtmp<-c(wtmp,gam2pngs(tmp.gam,i,j))
      } else {
        stmp<-adjmean(tmp.gam,i,j)
        wtmp<-c(wtmp,stmp[[2]])
        wtab<-c(wtab,stmp[[1]])
      }
    }
    wpng<-c(wpng,"<tr>",wtmp,"</tr>")
  } else {
      tmp.xx<-c(yvname[i],svname);
      if (nx>0) tmp.xx<-c(tmp.xx,xvname)
      if (!is.na(colvname[1])) tmp.xx<-c(tmp.xx,colvname)
      wd<-WD[,tmp.xx];
      wd<-wd[apply(is.na(wd),1,sum)==0,]
      fml<-paste(yvname[i],"~",paste(sxStr,collapse="+"),sep="")
      if (!is.na(colvname)) fml<-paste(fml,"+factor(",colvname,")",sep="")
      if (nx>0) fml<-paste(fml,"+",xv1,sep="")
      tmp.gam<-setgam(fml,i)
      wtab<-c(wtab,paste("</br></br>Outcome:",yb[i]))
      wtab<-c(wtab,gam2htmltable(tmp.gam))
      if (sum(slv==0)>0) {   
        wpng<-c(wpng,"<tr>",gam2pngs(tmp.gam,i,0),"</tr>")
      }
      if (sum(slv>0)>0) {   
        for (k in 1:ns) {
          if (slv[k]>0) {
            stmp<-adjmean(tmp.gam,i,k)
            wtab<-c(wtab,stmp[[1]])
            wpng<-c(wpng,stmp[[2]])
          }
        }
      }
  }
}
wpng<-c(wpng,"</table>")
w<-c(w,wpng,wtab)
w<-c(w,wd.subset)
w<-c(w,paste("</br>Created by EmpowerStats (www.empowerstats.com) and R on",Sys.Date()))
w<-c(w,"</body></html>")
fileConn<-file(paste(ofname,".htm",sep="")); writeLines(w, fileConn)
