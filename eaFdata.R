library(R2admb)

# setwd("M:/Fisheries Research/Stock Assessment & Data Analysis/Assessments/Completed Assessments/DEPM_SharkBay/Models/EstF_EG")
setwd("C:/admb/admb101-gcc452-win64/examples/admb/SnapInner2015")
area="ea" 
# ages area="ea"   catchDEPM area="ea"
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
IGweight=IGweight[,'weight']

plot(1:31,IGweight,"o")
for (i in 18:31) IGweight[i]=6.16
plot(1:31,IGweight,"o")
for (i in 8:20) IGweight[i]=(IGweight[i-1]+IGweight[i]+IGweight[i+1])/3
plot(1:31,IGweight,"o")

# determine age composition
xx$Year<-factor(xx$Year,levels=1997:2013)  # add in years with no data
xx$Cohorts<-factor(xx$Cohorts,levels=1:31)  # add in ages 1
IGagecomp=tapply(xx[,1],list(xx$Year,xx$Cohorts),length)
IGagecomp[is.na(IGagecomp)]=0
names(IGagecomp)=NULL
Nages=rowSums(IGagecomp)
names(Nages)=NULL

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
depm=subset(datt, Area==area, select=c('Area','DEPMM','DEPMsd','DEPMcv'))
depm$lnDEPM=log(depm$DEPMM)
depm$lnDEPM[depm$DEPMM==0]=0
DEPMM=depm[,'DEPMM']
DEPMsd=depm[,'DEPMsd']
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

yy=ageprop[50:70,]
ageprop5070=colSums(yy)
plot(1:31,ageprop5070,"o")
for (i in 8:20)
ageprop5070[i]=(ageprop5070[i-1]+ageprop5070[i]+ageprop5070[i+1])/3
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

fisharea=paste(area,"F.dat",sep="")
sink(fisharea)
# sink("eaF.dat")
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
DEPMM,"\n",
DEPMM,"\n",
DEPMM,"\n",
DEPMsd,"\n",
DEPMsd,"\n",
"\n#lnNdepm\n",
lnDEPMM,
"\n#lnNdepm_sd\n",
lnDEPMsd,	
"\n#Nage_samples\n",
Nages,
"\n#ages\n",
IGagecomp[1,],"\n",
IGagecomp[2,],"\n",
IGagecomp[3,],"\n",
IGagecomp[4,],"\n",
IGagecomp[5,],"\n",
IGagecomp[6,],"\n",
IGagecomp[7,],"\n",
IGagecomp[8,],"\n",
IGagecomp[9,],"\n",
IGagecomp[10,],"\n",
IGagecomp[11,],"\n",
IGagecomp[12,],"\n",
IGagecomp[13,],"\n",
IGagecomp[14,],"\n",
IGagecomp[15,],"\n",
IGagecomp[16,],"\n",
IGagecomp[17,],
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

fisharea=paste(area,"F.bat",sep="")
shell(fisharea)
