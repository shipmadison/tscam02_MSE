#### Plot Outlines 

# [1] EXTRACT THE DATA
########################
#### (1) Time Series ###
#######################
require(matrixStats)
require(Rgb)

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

########################################
#####   [2]    Make the Plots      #####
########################################


# Extract Raw Data Operating Model and show recruitment trends 


FemOnly<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "FemOnly_PlotTest", FileName= "FemOnly_Plots_EM")
# MaleR1<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "MaleOnlyR1_PlotTest", FileName= "MaleR1_Plots_EM")
# MaleR2<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "MaleOnlyR2_PlotTest", FileName= "MaleR2_Plots_EM")
# MaleR3<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "MaleOnlyR3_PlotTest", FileName= "MaleR3_Plots_EM")
# ELM50<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "ELM50_PlotTest", FileName= "ELM50_Plots_EM")
# Dimmer<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "Dimmer_PlotTest", FileName= "Dimmer_Plots_EM")
# Block<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "Block_PlotTest", FileName= "Block_Plots_EM")
# StatusQuo<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "StatusQuo_PlotTest", FileName= "StatusQuo_Plots_EM")


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


#########
#Recrutiment 

FemOnlyRec<-FemOnly$Rec
MaleR1Rec<-MaleR1$Rec
MaleR2Rec<-MaleR2$Rec
MaleR3Rec<-MaleR3$Rec
ELM50Rec<-ELM50$Rec
DimmerRec<-Dimmer$Rec
BlockRec<-Block$Rec
StatQuoRec<-StatusQuo$Rec


par(mfrow = c(1, 1))

plot(years, FemOnlyRec[1,], type = "l", main =" Recruitment Check")
lines(years,MaleR1Rec[1,], type = "l", col = col[1])
lines(years,MaleR2Rec[1,], type = "l", col = col[2])
lines(years,MaleR3Rec[1,], type = "l", col = col[3])
lines(years,ELM50Rec[1,], type = "l", col = col[4])
lines(years,DimmerRec[1,], type = "l", col = col[10])
lines(years,BlockRec[1,], type = "l", col = col[6])
lines(years,StatQuoRec[1,], type = "l", col = "grey20")

mtext("Years", side=1, line=1, outer = TRUE )
mtext("Millions of Individuals", side=2, line=1.25, outer = TRUE )


### Limit Plots --Time Series Process 




#Model Runs to be Plotted 

Fem_TimeSeries<-getLimitData(FemOnly)
MaleR1_TimeSeries<-getLimitData(MaleR1)
MaleR2_TimeSeries<-getLimitData(MaleR2)
MaleR3_TimeSeries<-getLimitData(MaleR3)

Dim_TimeSeries<-getLimitData(Dimmer)
Block_TimeSeries<-getLimitData(Block)

ELM50_TimeSeries<-getLimitData(ELM50)
Status_TimeSeries<-getLimitData(StatusQuo)

### TAC  Median####


#Preliminary Set-up

startYear<-2019
endYear<-2038

years<-seq(startYear, endYear)
xx<-c(years, rev(years))


## SECTION NOT DONE

#Fem
yyTAC10_fem<-c(Fem_TimeSeries$TAC10, rev(Fem_TimeSeries$TACMed))
yyTAC90_fem<-c(Fem_TimeSeries$TAC90, rev(Fem_TimeSeries$TACMed))

#MaleR1
yyTAC10_MaleR1<-c(MaleR1_TimeSeries$TAC10, rev(MaleR1_TimeSeries$TACMed))
yyTAC90_MaleR1<-c(MaleR1_TimeSeries$TAC90, rev(MaleR1_TimeSeries$TACMed))

#MaleR2
yyTAC10_MaleR2<-c(MaleR2_TimeSeries$TAC10, rev(MaleR2_TimeSeries$TACMed))
yyTAC90_MaleR2<-c(MaleR2_TimeSeries$TAC90, rev(MaleR2_TimeSeries$TACMed))

#MaleR3
yyTAC10_MaleR3<-c(MaleR3_TimeSeries$TAC10, rev(MaleR3_TimeSeries$TACMed))
yyTAC90_MaleR3<-c(MaleR3_TimeSeries$TAC90, rev(MaleR3_TimeSeries$TACMed))


#ELM50
yyTAC10_ELM50<-c(ELM50_TimeSeries$TAC10, rev(Dim_TimeSeries$TACMed))
yyTAC90_Dim<-c(Dim_TimeSeries$TAC90, rev(Dim_TimeSeries$TACMed))

#Block
yyTAC10_Block<-c(Block_TimeSeries$TAC10, rev(Block_TimeSeries$TACMed))
yyTAC90_Block<-c(Block_TimeSeries$TAC90, rev(Block_TimeSeries$TACMed))

#Dimmer
yyTAC10_Dim<-c(Dim_TimeSeries$TAC10, rev(Dim_TimeSeries$TACMed))
yyTAC90_Dim<-c(Dim_TimeSeries$TAC90, rev(Dim_TimeSeries$TACMed))

#Status Qup
yyTAC10_SQ<-c(Status_TimeSeries$TAC10, rev(Status_TimeSeries$TACMed))
yyTAC90_SQ<-c(Status_TimeSeries$TAC90, rev(Status_TimeSeries$TACMed))



### HCR Colors 

HCR1<-rgb(1, 0, 0,0.2)
HCR2_1<-rgb(0, 0, .1,0.2)
HCR2_2<-rgb(0, 0, .2,0.2)
HCR2_3<-rgb(0, 0, .3,0.2)
HCR6_3<-rgb(0, 0, .5,0.2)
HCR4<-rgb(1, 0, 1,0.2)
HCR5<-rgb(1, 0, .6,0.2)
HCR7<-rgb(1,1 ,0 ,0.3)


#Actual Plots 
par(mfrow = c(1, 1))
par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
plot(1,type="n",
     xlim=c(startYear,endYear),ylim=c(0,max(30)),
     yaxs="i", xaxs="i") #no gap between 0 and the axis

mtext(side=3, line=1," Median Total Allowable Catch" )
mtext(side=1, line=2, "Years")
mtext(side=2, line=2, "Thousands of Tons")


#Dimmer
polygon(xx,yyTAC10_Dim, col=HCR4, border = NA)
polygon(xx,yyTAC90_Dim,col=HCR4, border = NA)
lines(years, Dim_TimeSeries$TACMed, type = "l", col = "grey70")

#Dimmer
polygon(xx,yyTAC10_Dim, col=HCR4, border = NA)
polygon(xx,yyTAC90_Dim,col=HCR4, border = NA)
lines(years, Dim_TimeSeries$TACMed, type = "l", col = "grey70")

#Block
polygon(xx,yyTAC10_Block,col=HCR5, border = NA)
polygon(xx,yyTAC90_Block,col=HCR5, border = NA)
lines(years, Block_TimeSeries$TACMed, type = "l", col = "grey70")

#MaleR2
polygon(xx,yyTAC10_MaleR2,col=HCR2_2, border = NA)
polygon(xx,yyTAC90_MaleR2,col=HCR2_2, border = NA)
lines(years, MaleR2_TimeSeries$TACMed, type = "l",col = "grey70")

#StatusQuo
polygon(xx,yyTAC10_SQ,col=HCR7, border = NA)
polygon(xx,yyTAC90_SQ,col=HCR7, border = NA)
lines(years, Status_TimeSeries$TACMed, type = "l",col = "grey20",lwd=2)

legend("topright",c("HCR2_R2","HCR4","HCR5", "HCR7"),lty=c(1,1,1,1),col=c(HCR2_2,HCR4,HCR5, HCR7),lwd=c(6,6,6,6),cex=.8)


##########
#OFL
##########
  
  #Fem
  yyOFL10_fem<-c(Fem_TimeSeries$OFL10, rev(Fem_TimeSeries$OFLMed))
  yyOFL90_fem<-c(Fem_TimeSeries$OFL90, rev(Fem_TimeSeries$OFLMed))
  
  #MaleR1
  yyOFL10_MaleR1<-c(MaleR1_TimeSeries$OFL10, rev(MaleR1_TimeSeries$OFLMed))
  yyOFL90_MaleR1<-c(MaleR1_TimeSeries$OFL90, rev(MaleR1_TimeSeries$OFLMed))
  
  #MaleR2
  yyOFL10_MaleR2<-c(MaleR2_TimeSeries$OFL10, rev(MaleR2_TimeSeries$OFLMed))
  yyOFL90_MaleR2<-c(MaleR2_TimeSeries$OFL90, rev(MaleR2_TimeSeries$OFLMed))

  #MaleR3
  yyOFL10_MaleR3<-c(MaleR3_TimeSeries$OFL10, rev(MaleR3_TimeSeries$OFLMed))
  yyOFL90_MaleR3<-c(MaleR3_TimeSeries$OFL90, rev(MaleR3_TimeSeries$OFLMed))

  #Dimmer
  yyOFL10_Dim<-c(Dim_TimeSeries$OFL10, rev(Dim_TimeSeries$OFLMed))
  yyOFL90_Dim<-c(Dim_TimeSeries$OFL90, rev(Dim_TimeSeries$OFLMed))
  
  #Block
  yyOFL10_Block<-c(Block_TimeSeries$OFL10, rev(Block_TimeSeries$OFLMed))
  yyOFL90_Block<-c(Block_TimeSeries$OFL90, rev(Block_TimeSeries$OFLMed))
  
  #ELM50
  yyOFL10_ELM50<-c(ELM50_TimeSeries$OFL10, rev(ELM50_TimeSeries$OFLMed))
  yyOFL90_ELM50<-c(ELM50_TimeSeries$OFL90, rev(ELM50_TimeSeries$OFLMed))
  
  #Status Qup
  yyOFL10_SQ<-c(Status_TimeSeries$OFL10, rev(Status_TimeSeries$OFLMed))
  yyOFL90_SQ<-c(Status_TimeSeries$OFL90, rev(Status_TimeSeries$OFLMed))
  
  ### HCR Colors 
  
  HCR1<-rgb(.8, 0, 0,0.2)
  HCR2_1<-rgb(0, 0, 1,0.2)
  HCR2_2<-rgb(0, 0, 1,0.2)
  HCR2_3<-rgb(0, 0, 1,0.2)
  HCR6_3<-rgb(0, 0, .3,0.2)
  HCR4<-rgb(.5, 0, .5,0.2)
  HCR5<-rgb(.5, 0, .5,0.2)
  HCR7<-rgb(1,0 ,0 ,0.2)
  
  
  
  
  #Actual Plots 
  par(mfrow = c(1, 1))
  par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
  plot(1,type="n",
       xlim=c(startYear,endYear),ylim=c(0,max(75)),
       yaxs="i", xaxs="i") #no gap between 0 and the axis
  
  mtext(side=3, line=1," Median OFL" )
  mtext(side=1, line=2, "Years")
  mtext(side=2, line=2, "Thousands of Tons")
  
  #female
  polygon(xx,yyOFL10_fem,col=HCR1, border = NA)
  polygon(xx,yyOFL90_fem,col=HCR1, border = NA)
  lines(years, Fem_TimeSeries$OFLMed, type = "l",col = "grey70")
  
  #MaleR1
  polygon(xx,yyOFL10_MaleR1,col=HCR2_1, border = NA)
  polygon(xx,yyOFL90_MaleR1,col=HCR2_1, border = NA)
  lines(years, MaleR1_TimeSeries$OFLMed, type = "l",col = "grey70")
 
   #MaleR2
  polygon(xx,yyOFL10_MaleR2,col=HCR2_2, border = NA)
  polygon(xx,yyOFL90_MaleR2,col=HCR2_2, border = NA)
  lines(years, MaleR2_TimeSeries$OFLMed, type = "l",col = "grey70")
  
  #MaleR2
  polygon(xx,yyOFL10_MaleR3,col=HCR2_3, border = NA)
  polygon(xx,yyOFL90_MaleR3,col=HCR2_3, border = NA)
  lines(years, MaleR3_TimeSeries$OFLMed, type = "l",col = "grey70")
  
  
  #Dimmer
  polygon(xx,yyOFL10_Dim, col=HCR4, border = NA)
  polygon(xx,yyOFL90_Dim,col=HCR4, border = NA)
  lines(years, Dim_TimeSeries$OFLMed, type = "l", col = "grey70")
  
  #Block
  polygon(xx,yyOFL10_Block,col=HCR5, border = NA)
  polygon(xx,yyOFL90_Block,col=HCR5, border = NA)
  lines(years, Block_TimeSeries$OFLMed, type = "l", col = "grey70")
  
  #ELM
  polygon(xx,yyOFL10_ELM50,col=HCR6_3, border = NA)
  polygon(xx,yyOFL90_ELM50,col=HCR6_3, border = NA)
  lines(years, ELM50_TimeSeries$OFLMed, type = "l",col = "grey70")
  
  #StatusQuo
  polygon(xx,yyOFL10_SQ,col=HCR7, border = NA)
  polygon(xx,yyOFL90_SQ,col=HCR7, border = NA)
  lines(years, Status_TimeSeries$OFLMed, type = "l",col = "grey20", lwd=3)
  
  legend("topright",c("HCR1", "HCR2_R1","HCR2_R2","HCR2_R3", "HCR4","HCR5", "HCR6_50", "HCR7"),lty=c(1,1,1,1,1,1,1,1),col=c(HCR1, HCR2_1, HCR2_2, HCR2_3, HCR4,HCR5, HCR6_3, HCR7),lwd=c(6,6,6,6,6,6,6,6),cex=.8)
  

#$$$$$$$$$$$$$
#  BIO DATA 
#$$$$$$$$$$$$$
  
  ### Recruitment ####
  
  #Model Runs to be Plotted 
  Dim_TimeSeries2<-getBioData(Dimmer)
  Block_TimeSeries2<-getBioData(Block)
  MaleR2_TimeSeries2<-getBioData(MaleR2)
  Status_TimeSeries2<-getBioData(StatusQuo)
  
  
  
  
  
  #Dimmer
  yyRec10_Dim<-c(Dim_TimeSeries2$Rec10, rev(Dim_TimeSeries2$RecMed))
  yyRec90_Dim<-c(Dim_TimeSeries2$Rec90, rev(Dim_TimeSeries2$RecMed))
  
  #Block
  yyRec10_Block<-c(Block_TimeSeries2$Rec10, rev(Block_TimeSeries2$RecMed))
  yyRec90_Block<-c(Block_TimeSeries2$Rec90, rev(Block_TimeSeries2$RecMed))
  
  #MaleR2
  yyRec10_MaleR2<-c(MaleR2_TimeSeries2$Rec10, rev(MaleR2_TimeSeries2$RecMed))
  yyRec90_MaleR2<-c(MaleR2_TimeSeries2$Rec90, rev(MaleR2_TimeSeries2$RecMed))
  
  #Status Qup
  yyRec10_SQ<-c(Status_TimeSeries2$Rec10, rev(Status_TimeSeries2$RecMed))
  yyRec90_SQ<-c(Status_TimeSeries2$Rec90, rev(Status_TimeSeries2$RecMed))
  
  
  
  ### HCR Colors 
  
  HCR2_2<-rgb(0, 0, 1,0.2)
  HCR4<-rgb(.5, 0, .5,0.2)
  HCR5<-rgb(.6, 0, .6,0.2)
  HCR7<-rgb(1,0 ,0 ,0.3)
  
  
  #Actual Plots 
  par(mfrow = c(1, 1))
  par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
  plot(1,type="n",
       xlim=c(startYear,endYear),ylim=c(0,max(2000)),
       yaxs="i", xaxs="i") #no gap between 0 and the axis
  
  mtext(side=3, line=1," Median Rec" )
  mtext(side=1, line=2, "Years")
  mtext(side=2, line=2, "Thousands of Tons")
  
  
  #Dimmer
  polygon(xx,yyRec10_Dim, col=HCR4, border = NA)
  polygon(xx,yyRec90_Dim,col=HCR4, border = NA)
  lines(years, Dim_TimeSeries2$RecMed, type = "l", col = "grey70")
  
  #Block
  polygon(xx,yyRec10_Block,col=HCR5, border = NA)
  polygon(xx,yyRec90_Block,col=HCR5, border = NA)
  lines(years, Block_TimeSeries$RecMed, type = "l", col = "grey70")
  
  #MaleR2
  polygon(xx,yyRec10_MaleR2,col=HCR2_2, border = NA)
  polygon(xx,yyRec90_MaleR2,col=HCR2_2, border = NA)
  lines(years, MaleR2_TimeSeries$RecMed, type = "l",col = "grey70")
  
  #StatusQuo
  polygon(xx,yyRec10_SQ,col=HCR7, border = NA)
  polygon(xx,yyRec90_SQ,col=HCR7, border = NA)
  lines(years, Status_TimeSeries$RecMed, type = "l",col = "grey20", lwd=2)
  
  legend("topright",c("HCR2_R2","HCR4","HCR5", "HCR7"),lty=c(1,1,1,1),col=c(HCR2_2,HCR4,HCR5, HCR7),lwd=c(6,6,6,6),cex=.8)
  
  
  
 
  


#############################################
##### MMB > MMBave by mean TAC
############################################## 
#5millions lbs = 2.2 thousand metrics tons 
#10 = 4.5
#15 = 6.8
#20 = 9

#Ave MMB
#Ave MFB

aveMMB <-5.1764925267e+001 
closeMMB<-0.25*aveMMB

aveMFB <-2.3514569603e+001


#Fem 
MMB_fem<-FemOnly$MMB
counter<-0

for(i in 1:sims){
  for(j in 1:length(years)){
    print(MMB_dim[i,j])
    if(MMB_fem[i,j] > (aveMMB)){
      counter<-counter+1
      probMMB_fem<-counter/length(MMB_fem)
    }
  }
}

#MaleR1
MMB_MR1<-MaleR1$MMB
counter<-0

for(i in 1:sims){
  for(j in 1:length(years)){
    print(MMB_MR1[i,j])
    if(MMB_MR1[i,j] > aveMMB){
      counter<-counter+1
      probMMB_MR1<-counter/length(MMB_MR1)
    }
  }
}

#MaleR2
MMB_MR2<-MaleR2$MMB
counter<-0

for(i in 1:sims){
  for(j in 1:length(years)){
    print(MMB_MR2[i,j])
    if(MMB_MR2[i,j] > aveMMB){
      counter<-counter+1
      probMMB_MR2<-counter/length(MMB_MR2)
    }
  }
}

#MaleR3
MMB_MR3<-MaleR3$MMB
counter<-0

for(i in 1:sims){
  for(j in 1:length(years)){
    print(MMB_MR3[i,j])
    if(MMB_MR3[i,j] > aveMMB){
      counter<-counter+1
      probMMB_MR3<-counter/length(MMB_MR3)
    }
  }
}

#Dimmer 
MMB_dim<-Dimmer$MMB
counter<-0

for(i in 1:sims){
  for(j in 1:length(years)){
    print(MMB_dim[i,j])
    if(MMB_dim[i,j] > (aveMMB)){
      counter<-counter+1
      probMMB_dim<-counter/length(MMB_dim)
    }
  }
}

#Block
MMB_blk<-Block$MMB
counter<-0

for(i in 1:sims){
  for(j in 1:length(years)){
    print(MMB_blk[i,j])
    if(MMB_blk[i,j] > aveMMB){
      counter<-counter+1
      probMMB_blk<-counter/length(MMB_blk)
    }
  }
}

#ELM50
MMB_ELM50<-ELM50$MMB
counter<-0

for(i in 1:sims){
  for(j in 1:length(years)){
    print(MMB_ELM50[i,j])
    if(MMB_ELM50[i,j] > aveMMB){
      counter<-counter+1
      probMMB_ELM50<-counter/length(MMB_ELM50)
    }
  }
}



#Status Quo
MMB_Stat<-StatusQuo$MMB
counter<-0

for(i in 1:sims){
  for(j in 1:length(years)){
    print(MMB_Stat[i,j])
    if(MMB_Stat[i,j] > aveMMB){
      counter<-counter+1
      probMMB_Stat<-counter/length(MMB_Stat)
    }
  }
}

  

# GET HELP on Prob logic 

# Get Values for Mean TAC

TACave_fem<-getLimitVals(FemOnly)[3]
TACmin_fem<-getLimitVals(FemOnly)[4]
TACmax_fem<-getLimitVals(FemOnly)[5]

TACave_MR1<-getLimitVals(MaleR1)[3]
TACmin_MR1<-getLimitVals(MaleR1)[4]
TACmax_MR1<-getLimitVals(MaleR1)[5]

TACave_MR2<-getLimitVals(MaleR2)[3]
TACmin_MR2<-getLimitVals(MaleR2)[4]
TACmax_MR2<-getLimitVals(MaleR2)[5]

TACave_MR3<-getLimitVals(MaleR3)[3]
TACmin_MR3<-getLimitVals(MaleR3)[4]
TACmax_MR3<-getLimitVals(MaleR3)[5]

TACave_dim<-getLimitVals(Dimmer)[3]
TACmin_dim<-getLimitVals(Dimmer)[4]
TACmax_dim<-getLimitVals(Dimmer)[5]

TACave_blk<-getLimitVals(Block)[3]
TACmin_blk<-getLimitVals(Block)[4]
TACmax_blk<-getLimitVals(Block)[5]

TACave_ELM50<-getLimitVals(ELM50)[3]
TACmin_ELM50<-getLimitVals(ELM50)[4]
TACmax_ELM50<-getLimitVals(ELM50)[5]

TACave_Stat<-getLimitVals(StatusQuo)[3]
TACmin_Stat<-getLimitVals(StatusQuo)[4]
TACmax_Stat<-getLimitVals(StatusQuo)[5]


# Colors again


### HCR Colors 

HCR1<-rgb(.8, 0, 0,0.4)
HCR2_1<-rgb(0, 0, 1,0.3)
HCR2_2<-rgb(0, 0, 1,0.6)
HCR2_3<-rgb(0, 0, 1,0.9)
HCR6_3<-rgb(0, 0, .3,0.8)
HCR4<-rgb(.5, 0, .5,0.4)
HCR5<-rgb(.5, 0, .5,0.8)
HCR7<-rgb(1,0 ,0 ,1)

#Actual Plots 
par(mfrow = c(1, 1))
par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
plot(1,type="n",
     xlim=c(0,1),ylim=c(0,max(30)),
     yaxs="i", xaxs="i") #no gap between 0 and the axis

mtext(side=3, line=1, "TAC Trade Example" )
mtext(side=1, line=2, "Proportion(MMB > aveMMB)")
mtext(side=2, line=2, "Average TAC")

abline(v =0.5, lty = 2, col = "grey70")

#Female
points(probMMB_fem,TACave_fem, pch = 19, col=HCR1, cex = 2 )
arrows(probMMB_fem, TACmin_fem, probMMB_fem, TACmax_fem, length=0.05, angle=90, code=3, col = HCR1)

#MaleR1
points(probMMB_MR1,TACave_MR1, pch = 19, col=HCR2_1,cex = 2 )
arrows(probMMB_MR1, TACmin_MR1, probMMB_MR1, TACmax_MR1, length=0.05, angle=90, code=3, col = HCR2_1)

#MaleR2
points(probMMB_MR2,TACave_MR2, pch = 19, col=HCR2_2,cex = 2 )
arrows(probMMB_MR2, TACmin_MR2, probMMB_MR2, TACmax_MR2, length=0.05, angle=90, code=3, col = HCR2_2)

#MaleR3
points(probMMB_MR3,TACave_MR3, pch = 19, col=HCR2_3,cex = 2 )
arrows(probMMB_MR3, TACmin_MR3, probMMB_MR3, TACmax_MR3, length=0.05, angle=90, code=3, col = HCR2_3)

#Dimmer
points(probMMB_dim,TACave_dim, pch = 19, col=HCR4,cex = 2 )
arrows(probMMB_dim, TACmin_dim, probMMB_dim, TACmax_dim, length=0.05, angle=90, code=3, col = HCR4)


#Block
points(probMMB_blk,TACave_blk, pch = 19, col=HCR5,cex = 2 )
arrows(probMMB_blk, TACmin_blk, probMMB_blk, TACmax_blk, length=0.05, angle=90, code=3, col = HCR5)


#ELM50
points(probMMB_ELM50,TACave_ELM50, pch = 19, col=HCR6_3,cex = 2 )
arrows(probMMB_ELM50, TACmin_ELM50, probMMB_ELM50, TACmax_ELM50, length=0.05, angle=90, code=3, col = HCR6_3)


#StatusQuo
points(probMMB_Stat,TACave_Stat, pch = 13, col=HCR7,cex = 2 )
arrows(probMMB_Stat, TACmin_Stat, probMMB_Stat, TACmax_Stat, length=0.05, angle=90, code=3, col = HCR7)


legend("topright",c("HCR1","HCR2_R1", "HCR2_R2" ,"HCR2_R3","HCR4","HCR5", "HCR6_50", "HCR7"),pch=c(19,19,19,19,19,19,19,13),col=c(HCR1, HCR2_1,HCR2_2,HCR2_3,HCR4,HCR5, HCR6_3,HCR7),cex=1)


########################################
#####   [3]    Make the Tables      ####
########################################

#Get Limit Vals 

Fem_ltab<-getLimitVals(FemOnly)
MaleR1_ltab<-getLimitVals(MaleR1)
MaleR2_ltab<-getLimitVals(MaleR2)
MaleR3_ltab<-getLimitVals(MaleR3)

Dim_ltab<-getLimitVals(Dimmer)
Block_ltab<-getLimitVals(Block)

ELM50_ltab<-getLimitVals(ELM50)
Status_ltab<-getLimitVals(StatusQuo)


LimitTab<-matrix(data= c(Fem_ltab,
                         MaleR1_ltab,
                         MaleR2_ltab,
                         MaleR3_ltab,
                         Dim_ltab,
                         Block_ltab,
                         ELM50_ltab,
                         Status_ltab), nrow = 8, ncol = 8, byrow = TRUE)

colnames(LimitTab)<-c("aveNoTac", "cv", "aveTAC", "minTAC", "maxTAC", "cv", "aveOFL", "cv")
row.names(LimitTab)<-c("Female", "10% MMB", "15% MMB","20% MMB","Dimmer","Block","ELM50","Status Quo")

write.csv(LimitTab, "LimitTable.csv")


#Get Bio Vals 

Fem_btab<-getBioVals(FemOnly)
MaleR1_btab<-getBioVals(MaleR1)
MaleR2_btab<-getBioVals(MaleR2)
MaleR3_btab<-getBioVals(MaleR3)

Dim_btab<-getBioVals(Dimmer)
Block_btab<-getBioVals(Block)

ELM50_btab<-getBioVals(ELM50)
Status_btab<-getBioVals(StatusQuo)

BioTab<-matrix(data= c(Fem_btab,
                         MaleR1_btab,
                         MaleR2_btab,
                         MaleR3_btab,
                         Dim_btab,
                         Block_btab,
                         ELM50_btab,
                         Status_btab), nrow = 8, ncol = 14, byrow = TRUE)

colnames(BioTab)<-c("aveB0", "cv", "aveBmsy", "cv", "aveMMB","cv", "aveMFB", "cv", "aveELMB", "cv", "aveRecAve", "cv", "aveRec", "cv")
row.names(BioTab)<-c("Female", "10% MMB", "15% MMB","20% MMB","Dimmer","Block","ELM50","Status Quo")

write.csv(BioTab, "BioTable.csv")


########################################
#Extra example Figs 



par(mfrow = c(1, 1))
par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
plot(1,type="n",
     xlim=c(startYear,endYear),ylim=c(0,max(OFL+2)),
     yaxs="i", xaxs="i") #no gap between 0 and the axis

mtext(side=3, line=2,"Total Allowable Catch" )
mtext(side=3, line=1, "ELM 50%")
mtext(side=1, line=2, "Years")
mtext(side=2, line=2, "1000's of Tons")

#TAC Error and Line 
polygon(xx,yy1_tac,col="grey70", border = NA)
polygon(xx,yy2_tac,col = "grey70", border = NA)

lines(startYear:endYear, aveTAC, type = "l",lwd = 1.5)

#OFL Line 
lines(startYear:endYear, aveOFL, type = "l", lty=5, col= "dark red")
#polygon(xx,yy1_ofl,col="grey95")
#polygon(xx,yy2_ofl,col = "grey95")

lines(startYear:endYear, (aveOFL*.80), type = "l", lty=3, col = "navy", lwd=1.5)

legend("topright",c("TAC","OFL", "ABC"),lty=c(1,5,3),col=c("black","dark red","navy"),lwd=c(1.5,1,1.5),cex=1.2)

