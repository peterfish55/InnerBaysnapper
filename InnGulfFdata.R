# library(R2admb)
setwd("C:/admb/admb101-gcc452-win64/examples/admb/SnapInner2015")
area="nw"      # HERE set area="ea" or "nw"  or "sw"
plot.out=0     # HERE set =0 or 1 for screen or =2 for tiff output 
#  ====================================================================================================
dat <- read.csv('InnerSBsnapperAges1997_2013.csv', header = T)
dat
xx=dat[(dat$Area==area) & (as.numeric(as.character(dat$Cohort)))>0,c('Year', 'LCFmm', 'TLmm', 'Cohorts')] 
yy=expand.grid(Year=1997:2013, Cohorts=1:31)
IGages=merge(yy,xx,by=c("Year","Cohorts"),all.x=TRUE)
IGages
IGages$weight=0.000000148*IGages$LCFmm^2.6703    # add a coulmn of weights
IGages
sy=aggregate(weight ~ Cohorts,mean,data=IGages)
ss=expand.grid(Cohorts=1:31)
IGweight=merge(ss,sy,by=c("Cohorts"),all.x=TRUE)
IGweight[is.na(IGweight)]=0
IGweight[1,2]=0.1
IGweight
IGweight=IGweight[,'weight']
IGweight

maxwt <- ifelse(area=='eg', 6.1, ifelse(area=='nw', 6.5, 6.8))
corrwt <- ifelse(area=='eg', 6.0, ifelse(area=='nw', 5.5, 5.7))
plot(1:31,IGweight,"o")
IGweight[17]<-corrwt
IGweight[IGweight==0]<-maxwt
plot(1:31,IGweight,"o")
for (i in 8:30) IGweight[i]=(IGweight[i-1]+IGweight[i]+IGweight[i+1])/3
plot(1:31,IGweight,"o")
for (i in 15:30) IGweight[i]=(IGweight[i-1]+IGweight[i]+IGweight[i+1])/3
plot(1:31,IGweight,"o")

# determine age composition
xx$Year<-factor(xx$Year,levels=1997:2013)  # add in years with no data
xx$Cohorts<-factor(xx$Cohorts,levels=1:31)  # add in ages 1
IGagecomp=tapply(xx[,1],list(xx$Year,xx$Cohorts),length)
IGagecomp[is.na(IGagecomp)]=0
names(IGagecomp)=NULL
Nages=rowSums(IGagecomp)
names(Nages)=NULL
## Write a comment here
nIGagecomp <- matrix('\n', ncol=ncol(IGagecomp)+1, nrow=nrow(IGagecomp))
nIGagecomp[1:nrow(IGagecomp),1:ncol(IGagecomp)] <- IGagecomp
nIGagecomp <- t(nIGagecomp)
      
# =================================  bring in CATCH data  ===========================
options(stringsAsFactors = FALSE)
dat <- read.csv('CA_20150513.csv', header = T)
dat
dat$CAM=as.numeric(as.character(dat$CAM))
catch=subset(dat, Area==area, select=c('Year', 'CAL', 'CAM', 'CAH','CAC','CAcv'))
catch$CAT=catch$CAM+catch$CAC
catch$CAsd=catch$CAM*catch$CAcv
Ocatch=catch[,'CAT']
Ocatchsd=catch[,'CAsd']
Ocatchsd[Ocatchsd==0]=1
# ===================Bring in  DEPM data  ============================================
options(stringsAsFactors = FALSE)
datt <- read.csv('DEPMdata - Zdaily 0.3 20150513.csv', header = T)
depm=subset(datt, Area==area, select=c('Area','DEPML','DEPMM','DEPMU','DEPMsd','DEPMcv'))
depm$lnDEPM=log(depm$DEPMM)
depm$lnDEPM[depm$DEPMM==0]=0
DEPML=depm[,'DEPML']
DEPMM=depm[,'DEPMM']
DEPMU=depm[,'DEPMU']
DEPMsd=depm[,'DEPMsd']
DEPMcv=depm[,'DEPMcv']
lnDEPMM=depm[,'lnDEPM']
lnDEPMsd=depm[,'DEPMcv']

# Pool Age-Composition data over all years;
options(stringsAsFactors = FALSE)

dat <- read.csv('InnerSBsnapperAges1997_2013.csv', header = T)
dd=dat[(dat$Area==area) & (as.numeric(as.character(dat$Cohorts)))>0,c('TLmm', 'Cohorts')] 
dd$TLcm=round((dd$TLmm+5)/10,0)
zz=dd[c('TLcm', 'Cohorts')]
zz
zz$TLcm<-factor(zz$TLcm,levels=1:90)  # add in years with no data
zz$Cohorts<-factor(zz$Cohorts,levels=1:31)  # add in ages 0,1
AgeLen=tapply(zz[,1],list(zz$TLcm,zz$Cohorts),length)
AgeLen[is.na(AgeLen)]=0
# Determine of TL's in each age class 
aa<-colSums(AgeLen)
aa[aa==0]<-0.1
ageprop<-t(round(t(AgeLen)/aa,5))
ageprop
# determine proportions retained pre year 2000 
zz=ageprop[1:44,]
zz[1,1:2]=1
ageprop45=1-colSums(zz)
plot(1:31,ageprop45,"o")
ageprop45[2]=0
plot(1:31,ageprop45,"o")
names(ageprop45)=NULL
ageprop45
yy=ageprop[50:70,]
ageprop5070=colSums(yy)
plot(1:31,ageprop5070,"o")
ageprop5070[ageprop5070==1]=0
plot(1:31,ageprop5070,"o")
ageprop5070[14] <- ifelse(area=='nw', 0.4)
plot(1:31,ageprop5070,"o")
for (i in 8:20) ageprop5070[i]=(ageprop5070[i-1]+ageprop5070[i]+ageprop5070[i+1])/3
plot(1:31,ageprop5070,"o")
for (i in 8:20) ageprop5070[i]=(ageprop5070[i-1]+ageprop5070[i]+ageprop5070[i+1])/3
plot(1:31,ageprop5070,"o")

names(ageprop5070)=NULL
# ==================== sexual maturity =======================================
# ====  VonB parameters circa 01/06/2011 supplied by gary 14/05/2014 =========
# ======  maturity parameters supplied by Gary  14/05/2014 ===================
K <- c(0.174,	0.151,	0.169)
t0 <- c(-0.022,	-0.084,	0.100)
Linf <- c(755,	728,	770)
LM50 <- c(350,400,420)
LM95 <- c(480,600,570)
LMslope=log(19)/(LM95-LM50)
# =========================================================================
pos <- ifelse(area=='eg', 1, ifelse(area=='nw', 2, 3))
LL=rep(NA, length(1:31))
SexMat=LL
LL = Linf[pos]*(1-exp(-K[pos]*(1:31+0.5-t0[pos])))
SexMat=1/(1+exp(-LMslope[pos]*(LL-LM50[pos])))
SexMat[SexMat>0.99] <- 1
# =======make the ADMB  .dat file ==================================
# =====  set future catch level  ===================================  
for (i in 32:38)  Ocatch[i]=10
for (i in 32:38)  Ocatchsd[i]=Ocatch[i]*catch$CAcv[i]



# old SexMat
# 0.0058,0.0544,0.2820,0.6637,0.8845,0.96,0.98,0.99,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,								
getwd()
fisharea=paste(area,"F.dat",sep="")
sink(fisharea)
cat("#Mode_quota\n",
1,10,
"\n#objective_function_weights\n",
1,1,1,1,1,																																		
"\n#phases_R*_Q_Ho_M_steep_select_dev_AnnualF\n",	
1,-1,-2,-2,-4,3,3,1,
"\n#fyear_dyear_lyear_pyear_fage_lage\n",	
1983,1997,2013,2020,1,31,																																						
"\n#legSize_PropF_Sdcorrect\n",																																					
3,0.5,0.69,																																					
"\n#Mprior_Mpriorsd_steepprior_steeppriorSD\n",																																				
0.12,0.005,0.75,0.02,																																				
"\n#qgrow\n",	
rep.int(1,2013-1983+1),
"\n#age\n",
1:31,
"\n#Sexmat\n",	
SexMat,								
"\n#Retention_pre2001\n",																																							
ageprop45,
"\n#Retention_post2000\n",
ageprop5070,
"\n#weightF\n",																																				
IGweight,
"\n#weightM\n",
IGweight,
"\n#obsCatch\n",																																			
Ocatch,
"\n#obsCatchsd\n",
Ocatchsd,
"\n#Ndepm\n",																																							
17,																																							
"\n#DEPM_SSB_1998_2013_LO_ME_HI_CV_SD\n",
DEPML,"\n",
DEPMM,"\n",
DEPMU,"\n",
DEPMcv,"\n",
DEPMsd,"\n",
"\n#lnNdepm\n",
lnDEPMM,
"\n#lnNdepm_sd\n",
lnDEPMsd,	
"\n#Nage_samples\n",
Nages,
"\n#ages\n",
nIGagecomp,
"\n#checknum\n",
1234)																																						
sink()

fisharea=paste(area,"F.pin",sep="")
sink(fisharea)
cat("#Rstar_Q_F0_M_H_A50_Aslope \n",
17.0,1.0,0.02,0.12,0.75,4.3,1.2, 
"\n#dev \n",
-0.08,0.11,-0.04,0.2,0.02,0.53,0.54,0.79,0.52,-0.09,-1.17,-1.08,-0.48,-0.4,-0.43,-0.15,0.31,0.36,0.02,0.35,-0.11,-0.06,0.46,0.61,0.48,-0.44,-0.44,-0.21,-0.09,-0.01,-0.01,0,0,0,0,0,0,0,
"\n#devSD_RdevSD\n", 
0.6,0.6,
"\n#F \n",
0.038,0.05,0.06,0.07,0.082,0.096,0.114,0.139,0.162,0.21,0.304,0.418,0.716,0.889,0.848,0.111,0.01,0.012,0.008,0.01,0.083,0.059,0.033,0.067,0.071,0.064,0.076,0.061,0.052,0.046,0.044,0.122,0.143,0.162,0.178,0.187,0.195,0.202,
"\n#EndofFile")
sink()
# =====================================================================================================
# =======================    RUN THE MODEL 8 TIMES TO GET PHI TO CONVERGE =============================
fisharea=paste(area,"F.bat",sep="")
shell(fisharea)

# ======================================================================================================
#  ==============  GET THE MODEL OUTPUTS SET UP  =======================================================
# rm(list=ls())
options(stringsAsFactors=FALSE)
library(zoo)
par(cex=0.8, cex.axis=0.8, cex.lab=0.8, cex.main=1)
allyears=1983:2020
ageyears=1997:2013
depmyears=1998:2013
startage=1
endage=31
model.name=area
predprop=read.table(paste("pred_prop_age_composition.dat",sep=""), header=FALSE)
head(predprop)
obsprop=read.table(paste("prop_age_composition.dat",sep=""), header=FALSE)
head(obsprop)
phi=read.table(paste("value_of_phi.dat",sep=""), header=FALSE)
phi
list.files()
std=read.table(paste(area,"F.std",sep=""), header=TRUE, fill=TRUE)
std[,5]=NULL
head(std)
unique(std$name)

lin=readLines(paste(area, "F.rep",sep=""))
obscatch=strsplit(lin[grep("Observed Catch", lin)+2], " ")[[1]]
obscatch=as.numeric(obscatch[obscatch!=""])
depm=strsplit(lin[grep("DEPM", lin)+1], " ")[[1]]
depm=as.numeric(depm[depm!=""])
depmlow=strsplit(lin[grep("DEPM", lin)+2], " ")[[1]]
depmlow=as.numeric(depmlow[depmlow!=""])
depmupp=strsplit(lin[grep("DEPM", lin)+3], " ")[[1]]
depmupp=as.numeric(depmupp[depmupp!=""])
obscatch=strsplit(lin[grep("Observed Catch", lin)+2], " ")[[1]]
obscatch=as.numeric(obscatch[obscatch!=""])
obscatch=strsplit(lin[grep("Observed Catch", lin)+2], " ")[[1]]
obscatch=as.numeric(obscatch[obscatch!=""])

retain1=strsplit(lin[grep("Retain1", lin)+1], " ")[[1]]
retain1=as.numeric(retain1[retain1!=""])
retain2=strsplit(lin[grep("Retain2", lin)+1], " ")[[1]]
retain2=as.numeric(retain2[retain2!=""])
Nsamp=strsplit(lin[grep("NageSamp", lin)+1], " ")[[1]]
Nsamp=as.numeric(Nsamp[Nsamp!=""])
SSage=strsplit(lin[grep("SSage", lin)+1], " ")[[1]]
SSage=as.numeric(SSage[SSage!=""])[1]
SSageWT=strsplit(lin[grep("Wt_age", lin)+1], " ")[[1]]
SSageWT=as.numeric(SSageWT[SSageWT!=""])[1]
phi=phi*SSageWT/SSage


SRRa=std[which(std$name=="SRRa"), 3]
SRRb=std[which(std$name=="SRRb"), 3]
VEstar=std[which(std$name=="VEstar"), 3]
Estar=std[which(std$name=="Estar"), 3]
Rstar=std[which(std$name=="Rstar"), 3]
SelA50=std[which(std$name=="SelA50"), 3]
SelSlope=std[which(std$name=="SelSlope"), 3]

##spawning biomass
unique(std$name)
ssb=std[which(std$name=="SSB"), 3:4]
ssb$uppCL=ssb$value+1.96*ssb$std
ssb$lowCL=ssb$value-1.96*ssb$std
ssb$year=allyears

##prop SSB
propSSB=std[which(std$name=="propSSB"), 3:4]
propSSB$uppCL=propSSB$value+1.96*propSSB$std
propSSB$lowCL=propSSB$value-1.96*propSSB$std
propSSB$year=allyears


##recruits
recdev=std[which(std$name=="dev" | std$name=="fut_dev"), 3:4]
recdev$uppCL=recdev$value+1.96*recdev$std
recdev$lowCL=recdev$value-1.96*recdev$std
recdev$year=allyears

unique(std$name)
rec=std[which(std$name=="Rec"), 3:4]
rec$uppCL=rec$value+1.96*rec$std
rec$lowCL=rec$value-1.96*rec$std
rec$year=allyears

##fish mort
fmort=std[which(std$name=="Fy"), 3:4]
fmort$uppCL=fmort$value+1.96*fmort$std
fmort$lowCL=fmort$value-1.96*fmort$std
fmort$year=allyears
##eff fish mort
eff.fmort=std[which(std$name=="Fy_eff"), 3:4]
eff.fmort$uppCL=eff.fmort$value+1.96*eff.fmort$std
eff.fmort$lowCL=eff.fmort$value-1.96*eff.fmort$std
eff.fmort$year=allyears

##catch
predcatch=std[which(std$name=="CA"), 3:4]
predcatch$uppCL=predcatch$value+1.96*predcatch$std
predcatch$lowCL=predcatch$value-1.96*predcatch$std
predcatch$year=allyears
##age comp
# =============================================================================
# ======  now plot model fit and save plots   =================================
# =============================================================================
#        save to single pdf file or individual tiff files
# if (plot.out==1){
#   ##plot to single pdf file
#   pdf(file=paste("Model Output ", model.name, ".pdf", sep=""), width=7, height=11, onefile=TRUE, paper="A4")
# }


if (plot.out==1 | plot.out==0){ par(mfrow=c(4,2), xaxs="i", yaxs="i", mar=c(4,4,2,2), mgp=c(2.5,1,0))
} else par(mfrow=c(1,1), xaxs="i", yaxs="i", mar=c(4,4,2,2), mgp=c(2.5,1,0))

# =============================================================================
## (1) plot SRR
# =============================================================================
SS=0:600
plot(SS, SS/(SRRa+SRRb*SS), ylim=c(0,1.05*max(SS/(SRRa+SRRb*SS))), type="l", xlab="Spawning Stock", ylab="Recruitment",main="Stock Recruitment Relationship", las=1)
points(VEstar, Rstar, col=2, cex=2, pch=16)
points(ssb$value[1], ssb$value[1]/(SRRa+SRRb*ssb$value[1]), col=4, cex=2, pch=4)

if (plot.out==2){dev.print(device=tiff, paste(model.name, " 1 Stock Recruitment Relationship", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}

# =============================================================================
## (2) plot SB
# =============================================================================

plot(allyears, ssb$value, type="l", ylim=c(0,1.05*max(ssb$uppCL, VEstar)), las=1, 
     xlab="Year", ylab="Spawning Biomass", main="Spawning Biomass")
polygon(c(allyears, rev(allyears)), c(ssb$lowCL, rev(ssb$uppCL)), col=8, border=NA)
lines(allyears, ssb$value, type="l", lwd=2)
lines(depmyears[depm>0], depm[depm>0], "p", pch=16)
arrows(depmyears[depm>0], depmlow[depm>0], depmyears[depm>0], depmupp[depm>0], code=3, angle=90, length=0.1, col=1, lty=1)
abline(h=c(0.2, 0.3, 0.4, 1)*VEstar, col=2, lwd=c(1,1,2, 1), lty=c(1, 1, 1, 2))
# mtext(paste("q=",round(depm_catchability, 2), sep=""), side=3, line=-1.5, adj=0.99, cex=0.8)
box()

# mtext(paste("H0=",round(initial_harvest_fraction, 2), sep=""), side=3, line=-2.8, adj=0.99, cex=0.8)

if (plot.out==2){dev.print(device=tiff, paste(model.name, " 2 Spawning Biomass", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}

plot(allyears, propSSB$value, type="l", ylim=c(0,100), las=1, 
     xlab="Year", ylab="Spawning Biomass (% Unfished)", main="Spawning Biomass (% Unfished)")
polygon(c(allyears, rev(allyears)), c(propSSB$lowCL, rev(propSSB$uppCL)), col=8, border=NA)
lines(allyears, propSSB$value, type="l", lwd=2)
lines(depmyears[depm>0], depm[depm>0]/VEstar*100, "p", pch=16)
arrows(depmyears[depm>0], depmlow[depm>0]/VEstar*100, depmyears[depm>0], depmupp[depm>0]/VEstar*100, code=3, angle=90, length=0.1, col=1, lty=1)
abline(h=c(0.2, 0.3, 0.4)*100, col=2, lwd=c(1,1,2), lty=c(1, 1, 1))
# mtext(paste("q=",round(depm_catchability, 2), sep=""), side=3, line=-1.5, adj=0.99, cex=0.8)
box()

# mtext(paste("H0=",round(initial_harvest_fraction, 2), sep=""), side=3, line=-2.8, adj=0.99, cex=0.8)

if (plot.out==2){dev.print(device=tiff, paste(model.name, " 2a Spawning Biomass Prop", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}

# =============================================================================
## (3) catch
# =============================================================================

plot(allyears, predcatch$value, type="l", ylim=c(0,ceiling(max(40/10))*10), las=1, xlab="Year", ylab="Annual Yield", main="Annual Yield")
polygon(c(allyears, rev(allyears)), c(predcatch$lowCL, rev(predcatch$uppCL)), col=8, border=NA)
lines(allyears, predcatch$value, type="l", lwd=2)
points(allyears, obscatch, pch=16, col=2)
box()

if (plot.out==2){dev.print(device=tiff, paste(model.name, " 3 Annual Yield", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}

## =============================================================================
## (5) Selectivity and Retention
## =============================================================================
ages=startage:endage
plot(ages, 1/ (1.0 + exp(-SelSlope * (ages - SelA50) )), type="l", lwd=2, xlim=c(0,endage), ylim=c(0,1), las=1,
     xlab="Age", ylab="Probability", main="Selectivity")
lines(ages, retain1[1:length(ages)], type="l", col=4)
lines(ages, retain2[1:length(ages)], type="l", col=2)
legend("right", c("Selectivity", "Retention to 2000", "Retention from 2001"), lty=1, lwd=c(2,1,1), col=c(1,4,2), bty="n", y.intersp=1.5)

if (plot.out==2){dev.print(device=tiff, paste(model.name, " 4 Selectivity and Retention", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}

# =============================================================================
## (5) fishing mortality
# =============================================================================
plot(allyears, fmort$value, type="l", ylim=c(0,1.0*max(fmort$uppCL)), las=1, xlab="Year", ylab="Fishing Mortality", main="Fully Selected Fishing Mortality")
polygon(c(allyears, rev(allyears)), c(fmort$lowCL, rev(fmort$uppCL)), col=8, border=NA)
lines(allyears, fmort$value, type="l", lwd=2)
abline(h=c(2/3,1,1.5)*0.12, col=2, lty=c(2,1,2))
box()
# mtext(paste("Q=",round(depm_catchability, 2), sep=""), side=3, line=-1, adj=0.99)
if (plot.out==2){dev.print(device=tiff, paste(model.name, " 5 Fishing Mortality", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}


# =============================================================================
## (6) effective fishing mortality
# =============================================================================
plot(allyears, eff.fmort$value, type="l", ylim=c(0,1.05*max(eff.fmort$uppCL)), las=1, 
     xlab="Year", ylab="Fishing Mortality", main="Effective Fishing Mortality")
polygon(c(allyears, rev(allyears)), c(eff.fmort$lowCL, rev(eff.fmort$uppCL)), col=8, border=NA)
lines(allyears, eff.fmort$value, type="l", lwd=2)
abline(h=c(2/3,1,1.5)*0.12, col=2, lty=c(2,1,2))
box()
# mtext(paste("Q=",round(depm_catchability, 2), sep=""), side=3, line=-1, adj=0.99)
if (plot.out==2){dev.print(device=tiff, paste(model.name, " 5a Effective Fishing Mortality", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}

## =============================================================================
## (5) recruitment
## =============================================================================


plot(allyears, recdev$value, type="l", ylim=c(-2.6,2.6), las=1, xlab="Year", ylab="Recruitment", main="Recruitment Deviation")
polygon(c(allyears, rev(allyears)), c(recdev$lowCL, rev(recdev$uppCL)), col=8, border=NA)
lines(allyears, recdev$value, type="l", lwd=2)
abline(h=0, col=1, lty=2)
abline(h=mean(recdev$value[recdev$year<=2013]), col=2)
box()

if (plot.out==2){dev.print(device=tiff, paste(model.name, " 7 Recruitment Deviation", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}

plot(allyears, rec$value, type="l", ylim=c(0, max(rec$uppCL)), las=1, xlab="Year", ylab="Recruitment", main="Recruitment")
polygon(c(allyears, rev(allyears)), c(rec$lowCL, rev(rec$uppCL)), col=8, border=NA)
lines(allyears, rec$value, type="l", lwd=2)
box()


if (plot.out==2){dev.print(device=tiff, paste(model.name, " 7a Recruitment", ".tiff", sep=""), width=170, height=120, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}



## =============================================================================
## (5) plot age composition data
## =============================================================================
##data
predprop
obsprop


par(mfrow=c(4,3), mar=c(4,3,2,2), xaxs="i", yaxs="i", mgp=c(2.5,1,0))
yr <- 2011
for (yr in ageyears){
  if (sum(obsprop[yr-min(ageyears)+1, ])>0){
    age_comp_index=which(ageyears==yr)
    data_index=which(allyears==yr)
    plot(startage:endage,obsprop[age_comp_index,], type="o", xlim=c(0,endage),
         ylim=c(0,max(predprop[age_comp_index,], obsprop[age_comp_index,])*1.05), 
         lty=2, main="", xlab="Age", ylab="", las=1)
    lines(startage:endage,predprop[age_comp_index,],col="red")
    mtext(yr, side=3, line=0.2, adj=0.5, cex=0.6, font=2)
    mtext(paste("N=",Nsamp[age_comp_index], sep=""), side=3, line=-1, adj=0.99, cex=0.5)
    mtext(paste("Neff=",round(phi*Nsamp[age_comp_index],0), sep=""), side=3, line=-2, adj=0.99, cex=0.5)
    legend("right", c("obs", "exp"), lty=c(NA,1), col=1:2, pch=c(1, NA), cex=0.8, bty="n")
  }
}

if (plot.out==2){dev.print(device=tiff, paste(model.name, " 8 Age Composition", ".tiff", sep=""), width=170, height=200, units="mm", 
                           res=300, pointsize = 14, compression ='lzw')}
if (plot.out==1) dev.off()
