# For Ben

# Hi Ben
# There are several things included in this series of functions.
#They all rely on the data that I gave you. 
# The first part are summary statistics functions. 
#You may already have that done, but this breaks it down by catagory 
# It's futhermore split into time series (getXData) and total statistics (getxVals).

#I've used the time series information for figures. MIND You, these sets of functions 
# Unfortunatly don't include the means. I did that later mannually and am in the process
# Of reworking my code to account for figures in a more efficient way.
# This first part will take you from line 25-510

# Below that I get into 3 plotting functions that may be useful to you 
# The first one is doSpaghettii Plots, which show you the trend for the data long term.
# This is more useful on an individual model based scale (with how it's currently set up)
#but you can probably tweak it to make comparisons easier. 
#(I'm working on modifying a lot of of this figure code, so we can give each other feedback)

# The last two provide the code for the figures that I sent you.
#If you don't want the pdfs and just want to look at it in the console, comment out
#the pdf() and dev.off() 

# Make sure you set your working directories! I had one for each HCR that I specified 
#within R for ease of access to data and plot output

# Feel free to ask questions! 

require(matrixStats)
require(Rgb)
library("expss")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#            Limit Data 
#    MedOFL, MedTAC, 10% and 90% CIs
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
getLimitData<-function(data = data){
  
  data<-data
  
  # Extract the info
  TAC<-data$TAC
  OFL<-data$OFL
  
  # Get the Time Series Medians, and the 10 and 90th percentiles 
  
  #TAC
  TACMed<-apply(TAC,2,median)
  TAC10<-colQuantiles(TAC, probs = .10)
  TAC90<-colQuantiles(TAC, probs = .90)
  
  #medians and 90% and 10% CI
  OFLMed<-apply(OFL,2,median)
  OFL10<-colQuantiles(OFL, probs = .10)
  OFL90<-colQuantiles(OFL, probs = .90)
  
  
  LimitData<-list("TACMed"=TACMed, "TAC10" = TAC10, "TAC90" = TAC90, "OFLMed"=OFLMed, "OFL10"=OFL10, "OFL90"=OFL90)
  
  return(LimitData)
}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#        Bio Data 
# MedB0, MedBmsy, MedMMB, MedMFB, MedELMB, MedaveRec, MedRec, 10% and 90% CIs
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
getBioData<-function(data = data){
  
  data<- data
  
  #Extract the Info
  B0<-data$B0
  Bmsy<-data$Bmsy
  MMB<-data$MMB
  MFB<-data$MFB
  ELMB<-data$ELMB
  ELMB_State<-data$ELMB_State
  RecAve<-data$AveRec
  Rec<-data$Rec
  
  # Get the Time Series Medians, and the 10 and 90th percentiles 
  
  #B0
  B0Med<-apply(B0,2,median)
  B010<-colQuantiles(B0, probs = .10)
  B090<-colQuantiles(B0, probs = .90)
  
  #Bmsy
  BmsyMed<-apply(Bmsy,2,median)
  Bmsy10<-colQuantiles(Bmsy, probs = .10)
  Bmsy90<-colQuantiles(Bmsy, probs = .90)
  
  #MMB
  MMBMed<-apply(MMB,2,median)
  MMB10<-colQuantiles(MMB, probs = .10)
  MMB90<-colQuantiles(MMB, probs = .90)
  
  #MFB
  MFBMed<-apply(MFB,2,median)
  MFB10<-colQuantiles(MFB, probs = .10)
  MFB90<-colQuantiles(MFB, probs = .90)
  
  #ELMB (State Definition)
  ELMBMed<-apply(ELMB_State,2,median)
  ELMB10<-colQuantiles(ELMB_State, probs = .10)
  ELMB90<-colQuantiles(ELMB_State, probs = .90)
  
  #AveRec
  RecAveMed<-apply(RecAve,2,median)
  RecAve10<-colQuantiles(RecAve, probs = .10)
  RecAve90<-colQuantiles(RecAve, probs = .90)
  
  #Rec
  RecMed<-apply(Rec,2,median)
  Rec10<-colQuantiles(Rec, probs = .10)
  Rec90<-colQuantiles(Rec, probs = .90)
  
  BioData <-list(
    "B0Med" = B0Med,
    "B010" = B010,
    "B090" = B090,
    "BmsyMed" = BmsyMed,
    "Bmsy10" = Bmsy10,
    "Bmsy90" = Bmsy90,
    "MMBMed" = MMBMed,
    "MMB10" = MMB10,
    "MMB90" = MMB90,
    "MFBMed" = MFBMed,
    "MFB10" = MFB10,
    "MFB90" = MFB90,
    "ELMBMed" = ELMBMed,
    "ELMB10" = ELMB10,
    "ELMB90" = ELMB90,
    "RecAveMed" = RecAveMed,
    "RecAve10" = RecAve10,
    "RecAve90" = RecAve90,
    "RecMed" = RecMed,
    "Rec10" = Rec10,
    "Rec90" = Rec90)
  
  return(BioData)
}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## Catch Data 
# MedELMC, MedMMCB, MedIMCB, 10% and 90% CIs
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


getCatchData<-function(data = data){
  
  data<- data
  
  #Extract the Info
  ELMC<-data$ELMC
  MMCB<-data$MMCB
  IMCB<-data$IMCB
  
  # Get the Time Series Medians, and the 10 and 90th percentiles 
  
  #ELM Catches 
  ELMCMed<-apply(ELMC,2,median)
  ELMC10<-colQuantiles(ELMC, probs = .10)
  ELMC90<-colQuantiles(ELMC, probs = .90)
  
  #Mature Male Catch Biomass 
  MMCBMed<-apply(MMCB,2,median)
  MMCB10<-colQuantiles(MMCB, probs = .10)
  MMCB90<-colQuantiles(MMCB, probs = .90)
  
  #Immature Male Catch Biomass 
  IMCBMed<-apply(IMCB,2,median)
  IMCB10<-colQuantiles(IMCB, probs = .10)
  IMCB90<-colQuantiles(IMCB, probs = .90)
  
  
  CatchData <-list(
    "ELMCMed" = ELMCMed,
    "ELMC10" = ELMC10,
    "ELMC90" = ELMC90,
    "MMCBMed" = MMCBMed,
    "MMBC10" = MMCB10,
    "MMCB90" = MMCB90,
    "IMCBMed" = IMCBMed,
    "IMCB10" = IMCB10,
    "IMCB90" = IMCB90)
  
  
  return(CatchData)
}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## Discard Data 
# MedELMD, MedMMDB, MedIMDB, MedMFDB, MedIFDB, 10% and 90% CIs
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

getDiscardData<-function(data = data){
  
  data<- data
  
  #Extract the Info
  MFDB<-data$MFDB
  IFDB<-data$IFDB
  ELMD<-data$ELMD
  MMDB<-data$MMDB
  IMDB<-data$IMDB
  
  # Get the Time Series Medians, and the 10 and 90th percentiles 
  
  #Mature Female Discard Biomass  
  MFDBMed<-apply(MFDB,2,median)
  MFDB10<-colQuantiles(MFDB, probs = .10)
  MFDB90<-colQuantiles(MFDB, probs = .90)
  
  #Immature Female Discard Biomass 
  IFDBMed<-apply(IFDB,2,median)
  IFDB10<-colQuantiles(IFDB, probs = .10)
  IFDB90<-colQuantiles(IFDB, probs = .90)
  
  #Exploitable Legal Male Discard Biomass 
  ELMDMed<-apply(ELMD,2,median)
  ELMD10<-colQuantiles(ELMD, probs = .10)
  ELMD90<-colQuantiles(ELMD, probs = .90)
  
  # Mature Male Discard Biomass 
  MMDBMed<-apply(MMDB,2,median)
  MMDB10<-colQuantiles(MMDB, probs = .10)
  MMDB90<-colQuantiles(MMDB, probs = .90)
  
  #Exploitable Legal Male Discard Biomass 
  IMDBMed<-apply(IMDB,2,median)
  IMDB10<-colQuantiles(IMDB, probs = .10)
  IMDB90<-colQuantiles(IMDB, probs = .90)
  
  
  CatchData <-list(
    "MFDBMed" = MFDBMed,
    "MFDB10" = MFDB10,
    "MFDB90" = MFDB90,
    "IFDBMed" = IFDBMed,
    "IFDB10" = IFDB10,
    "IFDB90" = IFDB90,
    "ELMDMed"= ELMDMed,
    "ELMD10" = ELMD10,
    "ELMD90" = ELMD90,
    "MMDBMed" = MMCBMed,
    "MMDB10" = MMCB10,
    "MMDB90" = MMCB90,
    "IMDBMed" = IMCBMed,
    "IMDB10" = IMCB10,
    "IMDB90" = IMCB90)
  
  
  return(DiscardData)
}


#################################
#
###     (2) Full Averages 
#
################################

#Four Functions to extract information 

## Limit Data 
# AveOFL, AveClosedYears, Ave, min, and max TAC, all associated Errors 

getLimitVals<-function(data=data){
  
  data = data
  
  #Extract Data
  
  TACset<-data$`TAC Set`
  TAC<-data$TAC
  OFL<-data$OFL
  
  #Get Values 
  
  #Extract Averages and standard deviations
  
  #TAC_Set
  simNOTAC<-rowMeans(TACset)
  aveNOTAC<-mean(simNOTAC)
  sdNOTAC<-sd(simNOTAC)
  cvNOTAC<-sdNOTAC/aveNOTAC
  
  #TAC for all Years 
  simTacAve<-rowMeans(TAC)
  aveTAC<-mean(simTacAve)
  sdTAC<-sd(simTacAve)
  cvTAC<-sdTAC/aveTAC
  
  maxTAC<-max(TAC)
  minTAC<-min(TAC)
  
  # OFL for all Years 
  simOFLAve<-rowMeans(OFL)
  aveOFL<-mean(simOFLAve)
  sdOFL<-sd(simOFLAve)
  cvOFL<-sdOFL/aveOFL
  
  
  
  LimitVals<-c("aveNOTAC" = aveNOTAC,
               "cvNOTAC" = cvNOTAC,
               "aveTAC"= aveTAC,
               "minTAC" = minTAC,
               "maxTAC" = maxTAC,
               "cvTAC" = cvTAC,
               "aveOFL"= aveOFL,
               "cvOFL"= cvOFL)
  
  return(LimitVals)
  
}


## Bio Data 
# B0, Bmsy, aveMMB, aveMFB, aveELMB,aveRec, Rec, all associated Errors

getBioVals<-function(data=data){
  
  data= data
  
  # Extract info
  B0<-data$B0
  Bmsy<-data$Bmsy
  MMB<-data$MMB
  MFB<-data$MFB
  ELMB<-data$ELMB #preferred Males 
  ELMB_State<-data$ELMB_State
  RecAve<-data$AveRec
  Rec<-data$Rec
  
  #B0
  simB0<-rowMeans(B0)
  aveB0<-mean(simB0)
  sdB0<-sd(simB0)
  cvB0<-sdB0/aveB0
  
  #Bmsy
  simBmsy<-rowMeans(Bmsy)
  aveBmsy<-mean(simBmsy)
  sdBmsy<-sd(simBmsy)
  cvBmsy<-sdBmsy/aveBmsy
  
  #MMB
  simMMBAve<-rowMeans(MMB)
  aveMMB<-mean(simMMBAve)
  sdMMB<-sd(simMMBAve)
  cvMMB<-sdMMB/aveMMB
  
  #MFB
  simMFBAve<-rowMeans(MFB)
  aveMFB<-mean(simMFBAve)
  sdMFB<-sd(simMFBAve)
  cvMFB<-sdMFB/aveMFB
  
  #ELMB
  simELMBAve<-rowMeans(ELMB_State)
  aveELMB<-mean(simELMBAve)
  sdELMB<-sd(simELMBAve)
  cvELMB<-sdELMB/aveELMB
  
  #Average Recruitment
  simRecAve<-rowMeans(RecAve)
  aveRecAve<-mean(simRecAve)
  sdRec<-sd(simRecAve)
  cvRecAve<-sdRec/aveRecAve
  
  #Recruitment 
  simRec<-rowMeans(Rec)
  aveRec<-mean(simRec)
  sdRec<-sd(simRec)
  cvRec<-sdRec/aveRec
  
  
  getBioVals<-c("aveB0" = aveB0,
                "cvB0" = cvB0,
                "aveBmsy"= aveBmsy,
                "cvBmsy" = cvBmsy,
                "aveMMB" = aveMMB,
                "cvMMB" = cvMMB,
                "aveMFB" = aveMFB,
                "cvMFB" = cvMFB,
                "aveELMB" = aveELMB,
                "cvELMB" = cvELMB,
                "aveRecAve" = aveRecAve,
                "cvRecAve" = cvRecAve,
                "aveRec" = aveRec,
                "cvRec" = cvRec)
  
  return(getBioVals)
  
}


## Catch Data 
# aveELMC, aveMMCB, aveIMCB, all associated Errors 

getCatchVals<-function(data=data){
  
  data= data
  
  # get info 
  
  ELMC<-data$ELMC
  MMCB<-data$MMCB
  IMCB<-data$IMCB
  
  #Catches
  
  #Explot. Legal Male Catch Biomass 
  
  simELMCAve<-rowMeans(ELMC)
  aveELMC<-mean(simELMCAve)
  sdELMC<-sd(simELMCAve)
  cvELMC<-sdELMC/aveELMC
  
  #MMCB
  simMMCBAve<-rowMeans(MMCB)
  aveMMCB<-mean(simMMCBAve)
  sdMMCB<-sd(simMMCBAve)
  cvMMCB<-sdMMCB/aveMMCB
  
  #IMCB
  simIMCBAve<-rowMeans(IMCB)
  aveIMCB<-mean(simIMCBAve)
  sdIMCB<-sd(simIMCBAve)
  cvIMCB<-sdIMCB/aveIMCB
  
  CatchVals<-c("aveELMC" = aveELMC,
               "cvELMC" = cvELMC,
               "aveMMCB"= aveMMCB,
               "cvMMCB" = cvMMCB,
               "aveIMCB"= aveIMCB,
               "cvIMCB"= cvIMCB)
  
  return(CatchVals)
  
}


## Discard Data 
# aveELMD, aveMMDB, aveIMDB, aveMFDB, aveIFDB, all associated Errors 


getDiscardVals<-function(data=data){
  
  data= data
  
  # get info 
  
  MFDB<-data$MFDB
  IFDB<-data$IFDB
  ELMD<-data$ELMD
  MMDB<-data$MMDB
  IMDB<-data$IMDB
  
  #Discards
  
  #Explot. Legal Male Discard Biomass 
  simELMDAve<-rowMeans(ELMD)
  aveELMD<-mean(simELMDAve)
  sdELMD<-sd(simELMDAve)
  cvELMD<-sdELMD/aveELMD
  
  #MMDB
  simMMDBAve<-rowMeans(MMDB)
  aveMMDB<-mean(simMMDBAve)
  sdMMDB<-sd(simMMDBAve)
  cvMMDB<-sdMMDB/aveMMDB
  
  #IMDB
  simIMDBAve<-rowMeans(IMDB)
  aveIMDB<-mean(simIMDBAve)
  sdIMDB<-sd(simIMDBAve)
  cvIMDB<-sdIMDB/aveIMDB
  
  #MFDB
  simMFDBAve<-rowMeans(MFDB)
  aveMFDB<-mean(simMFDBAve)
  sdMFDB<-sd(simMFDBAve)
  cvMFDB<-sdMFDB/aveMFDB
  
  #IFDB
  simIFDBAve<-rowMeans(IFDB)
  aveIFDB<-mean(simIFDBAve)
  sdIFDB<-sd(simIFDBAve)
  cvIFDB<-sdIFDB/aveIFDB
  
  CatchVals<-c("aveELMD" = aveELMD,
               "cvELMD" = cvELMD,
               "aveMMDB"= aveMMDB,
               "cvMMDB" = cvMMDB,
               "aveIMDB"= aveIMDB,
               "cvIMDB"= cvIMDB,
               "aveMFDB"= aveMFDB,
               "cvMFDB" = cvMFDB,
               "aveIFDB"= aveIFDB,
               "cvIFDB"= cvIFDB)
  
  return(DiscardVals)
  
}


#####################
# Spaghetti Plots 
####################
doSpaghettiPlots<-function(data, dir = "dir", startYear, endYear, modTyp, sims, title = "main", fileout = "fileout"){
  
  dir<-dir
  #setwd(dir)
  
  data<-data
  TAC<-data$TAC
  OFL<-data$OFL
  Bmsy<-data$Bmsy
  MMB<-data$MMB
  MFB<-data$MFB
  ELMB<-data$ELMB_State
  Rec<-data$Rec
  Catch<-data$MMCB
  Discards<-data$IFDB + data$MFDB + data$MMDB + data$IMDB
  
  # Plot Details
  
  if(modTyp == 1){#estimation model 
    startYear<-startYear+1
    endYear<-endYear
  }else{ endYear<- endYear-1
  startYear<-StartYear
  }
  
  
  years<-seq(startYear,endYear,1)
  sims<-sims
  TACMAX<-max(TAC)
  OFLMAX<-max(OFL)
  BmsyMAX<-max(Bmsy)
  MMBMAX<-max(MMB)
  MFBMAX<-max(MFB)
  ELMBMAX<-max(ELMB)
  RecMAX<-max(Rec)
  CatchMAX<-max(Catch)
  DisMAX<-max(Discards)
  
  
  # Start Graphs 
  par(mfrow = c(3, 1))
  par(mar = c(2, 4, 2, 1), oma = c(2, 2, 2, 0))
  
  set.seed(1)
  col<-rgb(red = runif(24), green = runif(24), blue = runif(24), alpha=.8)
  
  
  #OFL
  plot(years, OFL[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,OFLMAX), main = "OFL", ylab = "1000's of Tons")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, OFL[i,], col = col[1])
  }
  
  #Bmsy
  plot(years, Bmsy[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,BmsyMAX), main = "Bmsy", ylab = "1000's of Tons")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, Bmsy[i,], col = col[12])
  }
  
  
  #Rec
  plot(years, Rec[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,RecMAX), main = "Recruitment", ylab = "1000's of individuals")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, Rec[i,], col = col[22])
  }
  
  mtext(title, side=3, line=0, outer = TRUE )
  
  
  #MMB
  plot(years, MMB[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,MMBMAX), main = "MMB", ylab = "1000's of Tons")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, MMB[i,], col = col[24])
  }
  
  #MFB
  plot(years, MFB[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,MMBMAX), main = "MFB", ylab = "1000's of Tons")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, MFB[i,], col = col[10])
  }
  
  #ELMB
  plot(years, ELMB[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,MMBMAX), main = "ELMB", ylab = "1000's of Tons")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, ELMB[i,], col = col[14])
  }
  
  mtext(title, side=3, line=0, outer = TRUE )
  
  
  #TAC
  
  plot(years, TAC[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,TACMAX), main = "TAC", ylab = "1000's of Tons")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, TAC[i,], col = col[3])
  }
  
  #Catch
  plot(years, Catch[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,CatchMAX), main = "Catch", ylab = "1000's of Tons")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, Catch[i,], col = col[16])
  }
  
  #Discards
  plot(years, Discards[1,], typ="l",xlim=c(startYear,endYear),ylim=c(0,TACMAX), main = "Discards", ylab = "1000's of Tons")
  box(col = "grey60")
  for(i in 2:sims){
    lines(years, Discards[i,], col = col[16])
  }
  
  mtext(title, side=3, line=0, outer = TRUE )
  
  
}

test<-doSpaghettiPlots(FemOnly,dir = NA, startYear = 2018, endYear = 2038, modTyp = 1, sims = 3, title = "Female Only Test", fileout = NA )


doTestPlots<-function(data=data, dir = "dir", HCR = "HCR", startYear, endYear, mod, file = "file"){
  
  data<-data #est
  HCR<-HCR
  setwd(dir)
  print(HCR)
  # Estimation Model = 1, operating model = 2
  # data<-readRDS("")
  
  # Model 1 = estimation model 
  if(mod==1){
    startYear<-startYear+1
    endYear<-endYear}else{
      startYear<-startYear
      endYear<-endYear-1}
  
  years<-seq(startYear, endYear) #estimation Mod
  xx<-c(years, rev(years))
  
  
  # Plot Index (One for OpMod, one for EstMod)
  # Limits and 2X1
  # TAC 
  # OFL
  
  #2X2
  # MMB
  # MFB
  # ELMB_State 
  # Rec
  
  #2X2
  # Top Row
  # ELMC
  # MMCB
  # Bottom Row
  # MFDB 
  # MMDB 
  
  
  ##################
  #  dataRead 
  #####################
  # Plot 1 
  #TAC
  TAC<-data$TAC
  TAC_Med<-apply(TAC, 2, median)
  TAC_Quants<- colQuantiles(TAC, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  #Polygons
  yyTAC5<-c(TAC_Quants[,1], rev(TAC_Med))
  yyTAC95<-c(TAC_Quants[,4], rev(TAC_Med))
  yyTAC25<-c(TAC_Quants[,2], rev(TAC_Med))
  yyTAC75<-c(TAC_Quants[,3], rev(TAC_Med))
  
  #OFL
  OFL<-data$OFL
  OFL_Med<-apply(OFL, 2, median)
  OFL_Quants<- colQuantiles(OFL, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  #Polygons
  yyOFL5<-c(OFL_Quants[,1], rev(OFL_Med))
  yyOFL95<-c(OFL_Quants[,4], rev(OFL_Med))
  yyOFL25<-c(OFL_Quants[,2], rev(OFL_Med))
  yyOFL75<-c(OFL_Quants[,3], rev(OFL_Med))
  
  
  # Plot 2
  #MMB 
  MMB<-data$MMB
  MMB_Med<-apply(MMB, 2, median)
  MMB_Quants<- colQuantiles(MMB, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  yyMMB5<-c(MMB_Quants[,1], rev(MMB_Med))
  yyMMB95<-c(MMB_Quants[,4], rev(MMB_Med))
  yyMMB25<-c(MMB_Quants[,2], rev(MMB_Med))
  yyMMB75<-c(MMB_Quants[,3], rev(MMB_Med))
  
  
  #MFB
  MFB<-data$MFB
  MFB_Med<-apply(MFB, 2, median)
  MFB_Quants<- colQuantiles(MFB, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  yyMFB5<-c(MFB_Quants[,1], rev(MFB_Med))
  yyMFB95<-c(MFB_Quants[,4], rev(MFB_Med))
  yyMFB25<-c(MFB_Quants[,2], rev(MFB_Med))
  yyMFB75<-c(MFB_Quants[,3], rev(MFB_Med))
  
  
  #ELMB
  ELMB<-data$ELMB_State
  ELMB_Med<-apply(ELMB, 2, median)
  ELMB_Quants<- colQuantiles(ELMB, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  yyELMB5<-c(ELMB_Quants[,1], rev(ELMB_Med))
  yyELMB95<-c(ELMB_Quants[,4], rev(ELMB_Med))
  yyELMB25<-c(ELMB_Quants[,2], rev(ELMB_Med))
  yyELMB75<-c(ELMB_Quants[,3], rev(ELMB_Med))
  
  #Rec
  Rec<-data$Rec
  Rec_Med<-apply(Rec, 2, median)
  Rec_Quants<- colQuantiles(Rec, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  yyRec5<-c(Rec_Quants[,1], rev(Rec_Med))
  yyRec95<-c(Rec_Quants[,4], rev(Rec_Med))
  yyRec25<-c(Rec_Quants[,2], rev(Rec_Med))
  yyRec75<-c(Rec_Quants[,3], rev(Rec_Med))
  
  
  # Plot 3 
  #ELMC
  ELMC<-data$ELMC
  ELMC_Med<-apply(ELMC, 2, median)
  ELMC_Quants<- colQuantiles(ELMC, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  yyELMC5<-c(ELMC_Quants[,1], rev(ELMC_Med))
  yyELMC95<-c(ELMC_Quants[,4], rev(ELMC_Med))
  yyELMC25<-c(ELMC_Quants[,2], rev(ELMC_Med))
  yyELMC75<-c(ELMC_Quants[,3], rev(ELMC_Med))
  
  #MMCB
  MMCB<-data$MMCB
  MMCB_Med<-apply(MMCB, 2, median)
  MMCB_Quants<- colQuantiles(MMCB, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  yyMMCB5<-c(MMCB_Quants[,1], rev(MMCB_Med))
  yyMMCB95<-c(MMCB_Quants[,4], rev(MMCB_Med))
  yyMMCB25<-c(MMCB_Quants[,2], rev(MMCB_Med))
  yyMMCB75<-c(MMCB_Quants[,3], rev(MMCB_Med))
  
  #MFDB
  MFDB<-data$MFDB
  MFDB_Med<-apply(MFDB, 2, median)
  MFDB_Quants<- colQuantiles(MFDB, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  yyMFDB5<-c(MFDB_Quants[,1], rev(MFDB_Med))
  yyMFDB95<-c(MFDB_Quants[,4], rev(MFDB_Med))
  yyMFDB25<-c(MFDB_Quants[,2], rev(MFDB_Med))
  yyMFDB75<-c(MFDB_Quants[,3], rev(MFDB_Med))
  
  
  #MMDB
  MMDB<-data$MMDB
  MMDB_Med<-apply(MMDB, 2, median)
  MMDB_Quants<- colQuantiles(MMDB, rows = NULL, cols = NULL, probs = c(.05, .25, .75, .95))
  yyMMDB5<-c(MMDB_Quants[,1], rev(MMDB_Med))
  yyMMDB95<-c(MMDB_Quants[,4], rev(MMDB_Med))
  yyMMDB25<-c(MMDB_Quants[,2], rev(MMDB_Med))
  yyMMDB75<-c(MMDB_Quants[,3], rev(MMDB_Med))
  
  ###########################  
  #years<-seq(2019, 2118)
  
  # Plot 1 
  
  pdf(file) 
  par(mfrow = c(2, 1))
  par(oma = c(0, 0, 2, 0))
  
  plot(years, TAC_Med, typ = "l", ylim = c(0,max(TAC_Quants[,4])), main = "TAC Median & Quantiles", ylab = "1000's of Tons")
  polygon(xx,yyTAC5, col="grey40", border = NA)
  polygon(xx,yyTAC95,col="grey40", border = NA)
  polygon(xx,yyTAC25, col="grey80", border = NA)
  polygon(xx,yyTAC75, col="grey80", border = NA)
  
  lines(years, TAC_Quants[,1], lty= 2, col = "grey20")
  lines(years, TAC_Quants[,2], lty= 2, col = "grey20")
  lines(years, TAC_Quants[,3], lty= 2, col = "grey20")
  lines(years, TAC_Quants[,4], lty= 2, col = "grey20")
  
  lines(years, TAC[1,], lty= 1, col = "navy")
  lines(years, TAC[2,], lty= 1, col = "blue")
  lines(years, TAC[3,], lty= 1)
  
  lines(years, TAC_Med, lty= 1, col = "red", lwd = 3)
  
  
  plot(years, OFL_Med, typ = "l", ylim = c(0,max(OFL_Quants[,4])), main = "OFL Median & Quantiles",ylab = "1000's of Tons")
  polygon(xx,yyOFL5, col="grey40", border = NA)
  polygon(xx,yyOFL95,col="grey40", border = NA)
  polygon(xx,yyOFL25, col="grey80", border = NA)
  polygon(xx,yyOFL75, col="grey80", border = NA)
  
  lines(years, OFL_Quants[,1], lty= 2, col = "grey60")
  lines(years, OFL_Quants[,2], lty= 2, col = "grey60")
  lines(years, OFL_Quants[,3], lty= 2, col = "grey60")
  lines(years, OFL_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, OFL[1,], lty= 1, col = "navy")
  lines(years, OFL[2,], lty= 1, col = "blue")
  lines(years, OFL[3,], lty= 1)
  
  lines(years, OFL_Med, lty= 1, col = "red", lwd = 3)
  
  mtext(HCR, side=3, line=1, outer = TRUE )
  
  
  
  #Plot 2 
  par(mfrow = c(2, 2))
  #par(mar = c(2, 1, 2, 1), oma = c(3, 3, 2, 1))
  
  plot(years, MMB_Med, typ = "l", ylim = c(0,max(MMB_Quants[,4])), main = "MMB Median & Quantiles",ylab = "1000's of Tons")
  polygon(xx,yyMMB5, col="grey40", border = NA)
  polygon(xx,yyMMB95,col="grey40", border = NA)
  polygon(xx,yyMMB25, col="grey80", border = NA)
  polygon(xx,yyMMB75, col="grey80", border = NA)
  lines(years, MMB_Quants[,1], lty= 2, col = "grey60")
  lines(years, MMB_Quants[,2], lty= 2, col = "grey60")
  lines(years, MMB_Quants[,3], lty= 2, col = "grey60")
  lines(years, MMB_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, MMB[1,], lty= 1, col = "navy")
  lines(years, MMB[2,], lty= 1, col = "blue")
  lines(years, MMB[3,], lty= 1)
  
  lines(years, MMB_Med, lty= 1, col = "red", lwd = 3)
  
  
  
  
  plot(years, MFB_Med, typ = "l", ylim = c(0,max(MFB_Quants[,4])), main = "MFB Median & Quantiles",ylab = "1000's of Tons")
  polygon(xx,yyMFB5, col="grey40", border = NA)
  polygon(xx,yyMFB95,col="grey40", border = NA)
  polygon(xx,yyMFB25, col="grey80", border = NA)
  polygon(xx,yyMFB75, col="grey80", border = NA)
  
  lines(years, MFB_Quants[,1], lty= 2, col = "grey60")
  lines(years, MFB_Quants[,2], lty= 2, col = "grey60")
  lines(years, MFB_Quants[,3], lty= 2, col = "grey60")
  lines(years, MFB_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, MFB[1,], lty= 1, col = "navy")
  lines(years, MFB[2,], lty= 1, col = "blue")
  lines(years, MFB[3,], lty= 1)
  
  lines(years, MFB_Med, lty= 1, col = "red", lwd = 3)
  
  
  
  
  plot(years, ELMB_Med, typ = "l", ylim = c(0,max(ELMB_Quants[,4])), main = "ELMB Median & Quantiles",ylab = "1000's of Tons")
  polygon(xx,yyELMB5, col="grey40", border = NA)
  polygon(xx,yyELMB95,col="grey40", border = NA)
  polygon(xx,yyELMB25, col="grey80", border = NA)
  polygon(xx,yyELMB75, col="grey80", border = NA)
  
  lines(years, ELMB_Quants[,1], lty= 2, col = "grey60")
  lines(years, ELMB_Quants[,2], lty= 2, col = "grey60")
  lines(years, ELMB_Quants[,3], lty= 2, col = "grey60")
  lines(years, ELMB_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, ELMB[1,], lty= 1, col = "navy")
  lines(years, ELMB[2,], lty= 1, col = "blue")
  lines(years, ELMB[3,], lty= 1)
  
  lines(years, ELMB_Med, lty= 1, col = "red", lwd = 3)
  
  plot(years, Rec_Med, typ = "l", ylim = c(0,max(Rec_Quants[,4])), main = "Rec Median & Quantiles",ylab = "1000's of crab")
  polygon(xx,yyRec5, col="grey40", border = NA)
  polygon(xx,yyRec95,col="grey40", border = NA)
  polygon(xx,yyRec25, col="grey80", border = NA)
  polygon(xx,yyRec75, col="grey80", border = NA)  
  
  
  lines(years, Rec_Quants[,1], lty= 2, col = "grey60")
  lines(years, Rec_Quants[,2], lty= 2, col = "grey60")
  lines(years, Rec_Quants[,3], lty= 2, col = "grey60")
  lines(years, Rec_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, Rec[1,], lty= 1, col = "navy")
  lines(years, Rec[2,], lty= 1, col = "blue")
  lines(years, Rec[3,], lty= 1)
  
  lines(years, Rec_Med, lty= 1, col = "red", lwd = 3)
  
  mtext(HCR, side=3, line=0, outer = TRUE )
  
  #Plot 3 
  
  plot(years, ELMC_Med, typ = "l", ylim = c(0,max(ELMC_Quants[,4])), main = "ELMC Median & Quantiles",ylab = "1000's of Tons")
  polygon(xx,yyELMC5, col="grey40", border = NA)
  polygon(xx,yyELMC95,col="grey40", border = NA)
  polygon(xx,yyELMC25, col="grey80", border = NA)
  polygon(xx,yyELMC75, col="grey80", border = NA) 
  
  
  lines(years, ELMC_Quants[,1], lty= 2, col = "grey60")
  lines(years, ELMC_Quants[,2], lty= 2, col = "grey60")
  lines(years, ELMC_Quants[,3], lty= 2, col = "grey60")
  lines(years, ELMC_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, ELMC[1,], lty= 1, col = "navy")
  lines(years, ELMC[2,], lty= 1, col = "blue")
  lines(years, ELMC[3,], lty= 1)
  
  lines(years, ELMC_Med, lty= 1, col = "red", lwd = 3)
  
  
  plot(years, MMCB_Med, typ = "l", ylim = c(0,max(MMCB_Quants[,4])), main = "MMCB Median & Quantiles",ylab = "1000's of Tons")
  polygon(xx,yyMMCB5, col="grey40", border = NA)
  polygon(xx,yyMMCB95,col="grey40", border = NA)
  polygon(xx,yyMMCB25, col="grey80", border = NA)
  polygon(xx,yyMMCB75, col="grey80", border = NA) 
  
  lines(years, MMCB_Quants[,1], lty= 2, col = "grey60")
  lines(years, MMCB_Quants[,2], lty= 2, col = "grey60")
  lines(years, MMCB_Quants[,3], lty= 2, col = "grey60")
  lines(years, MMCB_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, MMCB[1,], lty= 1, col = "navy")
  lines(years, MMCB[2,], lty= 1, col = "blue")
  lines(years, MMCB[3,], lty= 1)
  
  lines(years, MMCB_Med, lty= 1, col = "red", lwd = 3)
  
  
  plot(years, MMDB_Med, typ = "l", ylim = c(0,max(MMDB_Quants[,4])), main = "MMDB Median & Quantiles",ylab = "1000's of Tons")
  polygon(xx,yyMMDB5, col="grey40", border = NA)
  polygon(xx,yyMMDB95,col="grey40", border = NA)
  polygon(xx,yyMMDB25, col="grey80", border = NA)
  polygon(xx,yyMMDB75, col="grey80", border = NA) 
  
  lines(years, MMDB_Quants[,1], lty= 2, col = "grey60")
  lines(years, MMDB_Quants[,2], lty= 2, col = "grey60")
  lines(years, MMDB_Quants[,3], lty= 2, col = "grey60")
  lines(years, MMDB_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, MMDB[1,], lty= 1, col = "navy")
  lines(years, MMDB[2,], lty= 1, col = "blue")
  lines(years, MMDB[3,], lty= 1)
  
  lines(years, MMDB_Med, lty= 1, col = "red", lwd = 3)
  
  
  plot(years, MFDB_Med, typ = "l", ylim = c(0,max(MFDB_Quants[,4])), main = "MFDB Median & Quantiles",ylab = "1000's of Tons")
  polygon(xx,yyMFDB5, col="grey40", border = NA)
  polygon(xx,yyMFDB95,col="grey40", border = NA)
  polygon(xx,yyMFDB25, col="grey80", border = NA)
  polygon(xx,yyMFDB75, col="grey80", border = NA) 
  
  
  lines(years, MFDB_Quants[,1], lty= 2, col = "grey60")
  lines(years, MFDB_Quants[,2], lty= 2, col = "grey60")
  lines(years, MFDB_Quants[,3], lty= 2, col = "grey60")
  lines(years, MFDB_Quants[,4], lty= 2, col = "grey60")
  
  lines(years, MFDB[1,], lty= 1, col = "navy")
  lines(years, MFDB[2,], lty= 1, col = "blue")
  lines(years, MFDB[3,], lty= 1)
  
  lines(years, MFDB_Med, lty= 1, col = "red", lwd = 3)
  
  mtext(HCR, side=3, line=0, outer = TRUE )
  
  
  dev.off()
  
  
}



#2X2
#EstMod TAC vs OpMod Catch
#OFL, ABC from OpMod with TAC from EstMod
#Add MMB trajectories to previous plot 
#Add recruitment trajectories to previous plot 
#Keep y axis the same for all these plots

newPlots<-function(dataEst, dataOp, dir = "dir", HCR= "HCR", file ="file"){
  
  dataEst<-dataEst
  dataOp<-dataOp
  
  yearsEst<-seq(2019,2118)
  yearsOp<-seq(2018, 2117)
  
  #Data Read
  ################## 
  
  TAC<-dataEst$TAC
  TACAve<-colMeans(TAC)
  TACMed<-colQuantiles(TAC, probs = .5)
  
  catch<-dataOp$MMCB
  MMCBAve<-colMeans(catch)
  MMCBMed<-colQuantiles(catch, probs = .5)
  
  catch2<-dataOp$IMCB
  IMCBAve<-colMeans(catch2)
  IMCBMed<-colQuantiles(catch2, probs = .5)
  
  OFL<-dataOp$OFL
  OFLAve<-colMeans(OFL)
  OFLMed<-colQuantiles(OFL, probs = .5)
  OFL_Quants<-colQuantiles(OFL, probs<-c( .05, .25, .75, .95))
  
  
  ABC<-OFLMed*0.8
  
  Tac_OFL_ratio<-matrix(data=NA, nrow =100, ncol =100)
  Tac_ABC_ratio<-matrix(data=NA, nrow =100, ncol =100)
  
  
  for(j in 1:100){
    for(i in 1:100){
      
      Tac_OFL_ratio[i,j]<-TAC[i,j]/OFL[i,j]
      Tac_ABC_ratio[i,j]<-TAC[i,j]/(0.8*OFL[i,j])
      
    }
  }
  
  Tac_ABC_ratio_quants<-colQuantiles(Tac_ABC_ratio, probs = c(.05, .25, 0.5, .75, .95))
  Tac_OFL_ratio_quants<-colQuantiles(Tac_OFL_ratio, probs = c(.05, .25, 0.5, .75, .95))
  yy_Rabc_5<-c(Tac_ABC_ratio_quants[,1], rev(Tac_ABC_ratio_quants[,3]))
  yy_Rabc_95<-c(Tac_ABC_ratio_quants[,5], rev(Tac_ABC_ratio_quants[,3]))
  yy_Rabc_25<-c(Tac_ABC_ratio_quants[,2], rev(Tac_ABC_ratio_quants[,3]))
  yy_Rabc_75<-c(Tac_ABC_ratio_quants[,4], rev(Tac_ABC_ratio_quants[,3]))
  
  
  yy_Rofl_5<-c(Tac_OFL_ratio_quants[,1], rev(Tac_OFL_ratio_quants[,3]))
  yy_Rofl_95<-c(Tac_OFL_ratio_quants[,5], rev(Tac_OFL_ratio_quants[,3]))
  yy_Rofl_25<-c(Tac_OFL_ratio_quants[,2], rev(Tac_OFL_ratio_quants[,3]))
  yy_Rofl_75<-c(Tac_OFL_ratio_quants[,4], rev(Tac_OFL_ratio_quants[,3]))
  
  
  
  MMB<-dataOp$MMB
  MMBAve<-colMeans(MMB)
  MMBMed<-colQuantiles(MMB, probs = .5)
  
  Rec<-dataOp$Rec
  RecAve<-colMeans(Rec)
  RecMed<-colQuantiles(Rec, probs = .5)
  
  xxest<-c(yearsEst, rev(yearsEst))
  xxop<-c(yearsOp, rev(yearsOp))
  ######################  
  #Plots 
  ##################  
  
  setwd(dir)
  pdf(file)
  
  par(mfrow = c(2, 2))
  par(oma = c(0, 0, 2, 0 ))
  
  # TAC vs Catches 
  plot(yearsEst, TACMed, type = "l", lwd = 4, ylab = "1000's of Tons", ylim = c(0,25), main = "TAC vs Catches")
  lines(yearsOp, (MMCBMed+IMCBMed), col = "red")
  legend("topright",c("TAC Med (estMod)", "Catch Med (opMod)"),lty=c(1,1),col=c("black", "red"),lwd=c(4,1),cex=.8)
  
  #Ratios TAC:OFL
  plot(yearsOp, Tac_OFL_ratio_quants[,3], type = "l", lwd = 2, col= "red", ylab = "Ratio", ylim = c(0,2), main = "TAC:OFL")
  polygon(xxop,yy_Rofl_5, col="grey40", border = NA)
  polygon(xxest,yy_Rofl_95,col="grey40", border = NA)
  polygon(xxest,yy_Rofl_25, col="grey80", border = NA)
  polygon(xxest,yy_Rofl_75, col="grey80", border = NA) 
  
  abline(h=1, lty =2)
  lines(yearsOp, Tac_OFL_ratio_quants[,3],lwd = 2, col= "red")
  legend("topright","OFL median",lty=1,col="red",lwd=2,cex=.8)
  
  
  #TAC: ABC
  plot(yearsOp, Tac_ABC_ratio_quants[,3], type = "l", lwd = 2, col= "red", ylab = "Ratio", ylim = c(0,2), main = "TAC:ABC")
  
  polygon(xxop,yy_Rabc_5, col="grey40", border = NA)
  polygon(xxest,yy_Rabc_95,col="grey40", border = NA)
  polygon(xxest,yy_Rabc_25, col="grey80", border = NA)
  polygon(xxest,yy_Rabc_75, col="grey80", border = NA) 
  
  abline(h=1, lty =2)
  lines(yearsOp, Tac_ABC_ratio_quants[,3],lwd = 2, col= "red")
  legend("topright","ABC median",lty=1,col="red",lwd=2,cex=.8)
  
  # Other Plots
  plot(yearsOp, OFLMed, type = "l", lwd = 2, ylim = c(0,40),ylab = "1000's of Tons", main = "OFL,ABC,TAC med w/TAC sims")
  lines(yearsOp, TAC[1,], lty=3, lwd = 1, col = "grey40")
  lines(yearsOp, TAC[2,], lty=3, lwd = 1, col = "grey40")
  lines(yearsOp, TAC[3,], lty=3, lwd = 1, col = "grey40")
  lines(yearsOp, (OFLMed*.80), lty = 1, lwd =2, col = "navy")
  lines(yearsEst, TACMed, lty=1, lwd = 2, col = "dark red")
  
  
  mtext(HCR, side=3, line=1, outer = TRUE )
  
  dev.off()
  
  
  
  
}

# Example for run set-up
#ramp1 
setwd(dir)
HCR2_1_EstMod<-readRDS("HCR2_R1_EstMod")
HCR2_1_OpMod<-readRDS("HCR2_R1_OpMod")
doTestPlots(data= HCR2_1_EstMod, dir = plots, HCR = "HCR2_R1_EstMod", 2018, 2118, 1, file = "HCR2_R1_est.pdf")
doTestPlots(data= HCR2_1_OpMod, dir = plots, HCR = "HCR2_R1_OpMod", 2018, 2118, 2, file = "HCR2_R1_opmod.pdf")
newPlots(dataEst=HCR2_1_EstMod, dataOp=HCR2_1_OpMod, dir = plots, HCR= "HCR2_R1", file ="HCR2_R1_Comparisons.pdf")

