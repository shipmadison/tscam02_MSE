#test file 

setwd("C:/MSE_Runs")

### Extract Data from runs 
FemOnly<-getData(Type = "EstMod", 2018, 2028, SimNum = 3, folderName = "FemOnly_CPT_example", FileOutName= "test")
MaleOnlyModSurvBM<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "MaleOnlyModSurvBM_CPT_example", FileOutName= "test")
StatusQuo<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "StatusQuo_CPT_example", FileOutName= "test")
MaleOnlyR1<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "MaleOnlyR1_CPT_example", FileOutName= "test")
MaleOnlyR2<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "MaleOnlyR2_CPT_example", FileOutName= "test")
MaleOnlyR3<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "MaleOnlyR3_CPT_example", FileOutName= "test")
ABCdat<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "ABC20_CPT_example", FileOutName= "test")
Dimmerdat<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "Dimmer_CPT_example", FileOutName= "test")
Blockdat<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "Block_CPT_example", FileOutName= "test")
ELM30dat<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "ELM30_CPT_example", FileOutName= "test")
ELM40dat<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "ELM40_CPT_example", FileOutName= "test")
ELM50dat<-getData(Type = "EstMod", 2018, 2038, SimNum= 3, folderName = "ELM50_CPT_example", FileOutName= "test")
MaleSurveyBMdat<-getData(Type = "EstMod", 2018, 2028, SimNum= 3, folderName = "MaleOnlySurvBiomass_CPT_example", FileOutName= "test")



### New Runs Updated 

data<-MaleR2_Plots


# Plot results for 1 HCR
plotMSEresults<-function(modelType= "modelType", HCRname = "HCRname", startYear, endYear, data = data){

 if(modelType=="OpMod"){ 
   
   startYear<-startYear
   endYear<-endYear-1
 }
if(modelType=="EstMod"){
   startYear<-startYear+1
   endYear<-endYear
 }
  
  print(startYear)
  print(endYear)
  data = data
  
  years<-seq(startYear, endYear)
  
  HCRname<-HCRname
#Extract all the pertinent Data
#data<-data
  TACset<-data$`TAC Set`
  TAC<-data$TAC
  OFL<-data$OFL
  B0<-data$B0
  Bmsy<-data$Bmsy
  MMB<-data$MMB
  MFB<-data$MFB
  ELMB<-data$ELMB
  ELMB_State<-data$ELMB_State
  RecAve<-data$AveRec
  Rec<-data$Rec
  ELMC<-data$ELMC
  MMCB<-data$MMCB
  IMCB<-data$IMCB
  MFDB<-data$MFDB
  IFDB<-data$IFDB
  ELMD<-data$ELMD
  MMDB<-data$MMDB
  IMDB<-data$IMDB


#Extract Averages and standard deviations 
  
#TAC_Set
  ClosedYears<-apply(TACset,1,sum)
  aveClosedYears<-mean(ClosedYears)
  simNOTAC<-rowMeans(TACset)
  aveNOTAC<-mean(simNOTAC)
  sdNOTAC<-sd(simNOTAC)
  cvNOTAC<-sdNOTAC/aveNOTAC
  #medians and 90% and 10% CI
  TACSetMed<-apply(TACset,2,median)
  TAC10<-colQuantiles(TACset, probs = .10)
  TAC90<-colQuantiles(TACset, probs = .90)
  

#TAC
  aveTAC<-colMeans(TAC)
  sdTAC<-apply(TAC,2,sd)
  cvTAC<-sdTAC/aveTAC
  upperTAC<-aveTAC+sdTAC
  lowerTAC<-aveTAC-sdTAC
  print(cvTAC)
  #medians and 90% and 10% CI
  TACMed<-apply(TAC,2,median)
  TAC10<-colQuantiles(TAC, probs = .10)
  TAC90<-colQuantiles(TAC, probs = .90)
  

#OFL 
  aveOFL<-colMeans(OFL)
  sdOFL<-apply(OFL,2,sd)
  cvOFL<-sdOFL/aveOFL
  upperOFL<-aveOFL+sdOFL
  lowerOFL<-aveOFL-sdOFL
  #medians and 90% and 10% CI
  OFLMed<-apply(OFL,2,median)
  OFL10<-colQuantiles(OFL, probs = .10)
  OFL90<-colQuantiles(OFL, probs = .90)
  
#B0 
  aveB0<-colMeans(B0)
  sdB0<-apply(B0,2,sd)
  cvB0<-sdB0/aveB0
  upperB0<-aveB0+sdB0
  lowerB0<-aveB0-sdB0
  #medians and 90% and 10% CI
  B0Med<-apply(B0,2,median)
  B010<-colQuantiles(B0, probs = .10)
  B090<-colQuantiles(B0, probs = .90)
  
#Bmsy 
  aveBmsy<-colMeans(Bmsy)
  sdBmsy<-apply(Bmsy,2,sd)
  cvBmsy<-sdBmsy/aveBmsy
  upperBmsy<-aveBmsy+sdBmsy
  lowerBmsy<-aveBmsy-sdBmsy
  #medians and 90% and 10% CI
  BmsyMed<-apply(Bmsy,2,median)
  Bmsy10<-colQuantiles(Bmsy, probs = .10)
  Bmsy90<-colQuantiles(Bmsy, probs = .90)
  


#MMB

  aveMMB<-colMeans(MMB)
  sdMMB<-apply(MMB,2,sd)
  cvMMB<-sdMMB/aveMMB
  upperMMB<-aveMMB+sdMMB
  lowerMMB<-aveMMB-sdMMB
  #medians and 90% and 10% CI
  MMBMed<-apply(MMB,2,median)
  MMB10<-colQuantiles(MMB, probs = .10)
  MMB90<-colQuantiles(MMB, probs = .90)
  

#MFB

  aveMFB<-colMeans(MFB)
  sdMFB<-apply(MFB,2,sd)
  cvMFB<-sdMFB/aveMFB
  upperMFB<-aveMFB+sdMFB
  lowerMFB<-aveMFB-sdMFB
  #medians and 90% and 10% CI
  MFBMed<-apply(MFB,2,median)
  MFB10<-colQuantiles(MFB, probs = .10)
  MFB90<-colQuantiles(MFB, probs = .90)
  

#ELMB
  aveELMB<-colMeans(ELMB_State)
  sdELMB<-apply(ELMB_State,2,sd)
  cvELMB<-sdELMB/aveELMB
  upperELMB<-aveELMB+sdELMB
  lowerELMB<-aveELMB-sdELMB
  #medians and 90% and 10% CI
  ELMBMed<-apply(ELMB,2,median)
  ELMB10<-colQuantiles(ELMB, probs = .10)
  ELMB90<-colQuantiles(ELMB, probs = .90)

#Recruitment Average
  aveRecAve<-colMeans(RecAve)
  sdRecAve<-apply(RecAve,2,sd)
  cvRecAve<-sdRecAve/aveRecAve
  upperRecAve<-aveRecAve+sdRecAve
  lowerRecAve<-aveRecAve-sdRecAve
  print(cvRecAve)
  #medians and 90% and 10% CI
  RecAveMed<-apply(RecAve,2,median)
  RecAve10<-colQuantiles(RecAve, probs = .10)
  RecAve90<-colQuantiles(RecAve, probs = .90)
  
#Recruitment 
  aveRec<-colMeans(Rec)
  sdRec<-apply(Rec,2,sd)
  cvRec<-sdRec/aveRec
  upperRec<-aveRec+sdRec
  lowerRec<-aveRec-sdRec
  print(cvRec)
  #medians and 90% and 10% CI
  RecMed<-apply(Rec,2,median)
  Rec10<-colQuantiles(Rec, probs = .10)
  Rec90<-colQuantiles(Rec, probs = .90)
  
#Catches 
  #ELMB
  aveELMC<-colMeans(ELMC)
  sdELMC<-apply(ELMC,2,sd)
  cvELMC<-sdELMC/aveELMC
  upperELMC<-aveELMC+sdELMC
  lowerELMC<-aveELMC-sdELMC
  #medians and 90% and 10% CI
  ELMCMed<-apply(ELMC,2,median)
  ELMC10<-colQuantiles(ELMC, probs = .10)
  ELMC90<-colQuantiles(ELMC, probs = .90)
  
  #MMBC
  aveMMCB<-colMeans(MMCB)
  sdMMCB<-apply(MMCB,2,sd)
  cvMMCB<-sdMMCB/aveMMCB
  upperMMCB<-aveMMCB+sdMMCB
  lowerMMCB<-aveMMCB-sdMMCB
  #medians and 90% and 10% CI
  MMCBMed<-apply(MMCB,2,median)
  MMCB10<-colQuantiles(MMCB, probs = .10)
  MMCB90<-colQuantiles(MMCB, probs = .90)


  years<-seq(startYear, endYear)
  
  xx<-c(years, rev(years))


#TAC
  yy1_tac<-c(upperTAC, rev(aveTAC))
  yy2_tac<-c(lowerTAC, rev(aveTAC))

#OFL
  yy1_ofl<-c(upperOFL, rev(aveOFL))
  yy2_ofl<-c(lowerOFL, rev(aveOFL))
  
#B0
  yy1_B0<-c(upperB0, rev(aveB0))
  yy2_B0<-c(lowerB0, rev(aveB0))
  
#Bmsy
  yy1_Bmsy<-c(upperBmsy, rev(aveBmsy))
  yy2_Bmsy<-c(lowerBmsy, rev(aveBmsy))

#MMB
  yy1_mmb<-c(upperMMB, rev(aveMMB))
  yy2_mmb<-c(lowerMMB, rev(aveMMB))

#MFB
  yy1_mfb<-c(upperMFB, rev(aveMFB))
  yy2_mfb<-c(lowerMFB, rev(aveMFB))

#ELMB
  yy1_elmb<-c(upperELMB, rev(aveELMB))
  yy2_elmb<-c(lowerELMB, rev(aveELMB))

#RecruitmentAve
  yy1_RecAve<-c(upperRecAve, rev(aveRecAve))
  yy2_RecAve<-c(lowerRecAve, rev(aveRecAve))
  
  #Recruitment
  yy1_Rec<-c(upperRec, rev(aveRec))
  yy2_Rec<-c(lowerRec, rev(aveRec))
  
#Catches 
  #ELMC
  yy1_ELMC<-c(upperELMC, rev(aveELMC))
  yy2_ELMC<-c(lowerELMC, rev(aveELMC))
  
  #MMCB
  yy1_MMCB<-c(upperMMCB, rev(aveMMCB))
  yy2_MMCB<-c(lowerMMCB, rev(aveMMCB))
  
  
  ###########################################################
  #     Actual Figures Below
  ###########################################################
  

  graphLables<-c("TAC", "ELMB", "MMB", "MFB")
  par(mfrow = c(2, 2))
  par(cex = 0.6)
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))

  for (i in 1:4) {
   plot(1, axes = FALSE, type = "n",xlim=c(startYear,endYear),ylim=c(0,max(MMB+5)),yaxs="i")
    
  
  mtext(graphLables[i], side = 3, line = -2, adj = 0.1, cex = 1.2,col = "grey10")
           
    if (i %in% c(3, 4))
     axis(1, col = "grey40", col.axis = "grey20", at = seq(startYear,endYear, 5))
                                                             
     if (i %in% c(1,3))
       axis(2, col = "grey40", col.axis = "grey20", at = seq(0,max(MMB+5), 25))
                                                             
   box(col = "grey60")
   
   if( i ==1){ ### PLOT TAC 
     
     #standard dev
     polygon(xx,yy1_tac,col="grey70",border = NA)
     polygon(xx,yy2_tac,col = "grey70",border = NA)
     
     #TAC
     lines(startYear:endYear, aveTAC, type = "l", lwd= 1.5)
     
     #OFL
     lines(startYear:endYear, aveOFL, type = "l", lty=2, col= "dark green")
     
     legend("topright",c("TAC","OFL"),lty=c(1,2),col=c("black","dark green"),lwd=c(1.5,1), cex=1.2)
     #polygon(xx,yy1_ofl,col="grey95")
     #polygon(xx,yy2_ofl,col = "grey95")
     
     #lines(2019:2023, (aveOFL*.80), type = "l", lty=2, col = "blue")
     
   }
   
   if(i ==2){ #PLOT ELMB
     
     #Standard Dev
     polygon(xx,yy1_elmb,col="grey70",border = NA)
     polygon(xx,yy2_elmb,col = "grey70",border = NA)
     
     #ELMB
     lines(startYear:endYear, aveELMB, type = "l", lwd=1.5)
     
     #MMB Long term average 
     MMB_LTA<-51.7
     abline(h=MMB_LTA, lty=3, col= "blue")
     legend("topright",c("ELBM","MMBave"),lty=c(1,3),col=c("black","blue"),lwd=c(1.5,1),cex=1.2)
     #polygon(xx,yy1_ofl,col="grey95")
     #polygon(xx,yy2_ofl,col = "grey95")
     
     
   }
   if(i ==3){ #PLOT MMB
    
     #Standard Dev
     polygon(xx,yy1_mmb,col="grey70",border = NA)
     polygon(xx,yy2_mmb,col = "grey70",border = NA)
     
     #MMB
     lines(startYear:endYear, aveMMB, type = "l", lwd=1.5)
     
     
     MMB_LTA<-51.7
  
     abline(h=MMB_LTA, lty=3, col= "blue")
     
     legend("bottomleft",c("MMB","MMBave"),lty=c(1,3),col=c("black","blue"),lwd=c(1.5,1),cex=1.2)
     
     
   }
   if(i ==4){ #PLOT MFB
     
     #Standard Dev 
     polygon(xx,yy1_mfb,col="grey70", border = NA)
     polygon(xx,yy2_mfb,col = "grey70", border = NA)
     
     #MFB
     lines(startYear:endYear, aveMFB, type = "l", lwd=1.5)
     
     #MFBThreshold 
     MFB_LTA<-23.5
     abline(h=MFB_LTA, lty=3, col= "red")
     legend("topright",c("MFB","MFBave"),lty=c(1,3),col=c("black","red"),lwd=c(1.5,1),cex=1.2)
     
     
     ### CONSIDER PLOTTING THRESHOLD AVE INSTEAD OF OFL 
     
   }
  }
 mtext("years", side = 1, outer = TRUE, cex = 1, line = 2.2, col = "grey20")
        
 mtext("Millions of Lbs", side = 2, outer = TRUE, cex = 1, line = 2.2,col = "grey20")
        

#======================================================

#Single Plots
 
#=======================================================
 
par(mfrow = c(1, 1))
par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
plot(1,type="n",
     xlim=c(startYear,endYear),ylim=c(0,max(OFL+2)),
     yaxs="i", xaxs="i") #no gap between 0 and the axis

mtext(side=3, line=2,"Total Allowable Catch" )
mtext(side=3, line=1, HCRname)
mtext(side=1, line=2, "Years")
mtext(side=2, line=2, "Millions of Lbs")

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
####################
#      ELMB       #
####################
plot(1,type="n",
     xlim=c(startYear,endYear),ylim=c(0,max(ELMB+5)),
     yaxs="i", xaxs="i") #no gap between 0 and the axis

mtext(side=3, line=2,"Exploitable Legal Male Biomass" )
mtext(side=3, line=1, HCRname)
mtext(side=1, line=2, "Years")
mtext(side=2, line=2, "Millions of Lbs")

#TAC Error and Line 
polygon(xx,yy1_elmb,col="grey70", border = NA)
polygon(xx,yy2_elmb,col = "grey70", border = NA)

lines(startYear:endYear, aveELMB, type = "l", lwd = 1.5)

#OFL Line 
abline(h=MMB_LTA, lty=5, col= "dark red")
#polygon(xx,yy1_ofl,col="grey95")
#polygon(xx,yy2_ofl,col = "grey95")

legend("topright",c("ELMB", "MMBave"),lty=c(1,5),col=c("black","dark red"),lwd=c(1.5,1),cex=1.2)


###### MMB

plot(1,type="n",
     xlim=c(startYear,endYear),ylim=c(0,max(MMB+5)),
     yaxs="i", xaxs="i") #no gap between 0 and the axis

mtext(side=3, line=2,"Mature Male Biomass" )
mtext(side=3, line=1, HCRname)
mtext(side=1, line=2, "Years")
mtext(side=2, line=2, "Millions of Lbs")

#TAC Error and Line 
polygon(xx,yy1_mmb,col="grey70", border = NA)
polygon(xx,yy2_mmb,col = "grey70", border = NA)

lines(startYear:endYear, aveMMB, type = "l", lwd = 1.5)

#OFL Line 
abline(h=MMB_LTA, lty=5, col= " dark red")

legend("bottomleft",c("MMB","MMBave"),lty=c(1,5),col=c("black","dark red"),lwd=c(1.5,1),cex=1.2)



###### MFB

plot(1,type="n",
     xlim=c(startYear,endYear),ylim=c(0,max(MFB+5)),
     yaxs="i", xaxs="i") #no gap between 0 and the axis

mtext(side=3, line=2,"Mature Female Biomass" )
mtext(side=3, line=1, HCRname)
mtext(side=1, line=2, "Years")
mtext(side=2, line=2, "Millions of Lbs")

#TAC Error and Line 
polygon(xx,yy1_mfb,col="grey70", border = NA)
polygon(xx,yy2_mfb,col = "grey70", border = NA)

lines(startYear:endYear, aveMFB, type = "l", lwd = 1.5)

#OFL Line 
abline(h=MFB_LTA, lty=5, col= "dark red")

legend("bottomleft",c("MFB","MFBave"),lty=c(1,5),col=c("black","dark red"),lwd=c(1.5,1),cex=1.2)


#### Recruitment 
plot(1,type="n",
     xlim=c(startYear,endYear),ylim=c(80,max(Rec+5)),
     yaxs="i", xaxs="i", xlab = "Years", ylab = "Rec") #no gap between 0 and the axis

mtext(side=3, line=2,"Recruitment Averages" )

polygon(xx,yy1_Rec,col="grey70", border = NA)
polygon(xx,yy2_Rec,col = "grey70", border = NA)

lines(years, aveRec, lwd= 2.5)

###################
#####Catch Plots #####

plot(1,type="n",
     xlim=c(startYear,endYear),ylim=c(0,max(ELMC+5)),
     yaxs="i", xaxs="i") #no gap between 0 and the axis

mtext(side=3, line=2,"Catch of Legal Male Biomass" )
mtext(side=3, line=1, HCRname)
mtext(side=1, line=2, "Years")
mtext(side=2, line=2, "1000's of Tonnes")

#TAC Error and Line 
polygon(xx,yy1_ELMC,col="grey70", border = NA)
polygon(xx,yy2_ELMC,col = "grey70", border = NA)

lines(startYear:endYear, aveELMC, type = "l", lwd = 1.5)

#OFL Line 
lines(startYear:endYear, aveMMCB, lty=5, col= "dark red")
#polygon(xx,yy1_ofl,col="grey95")
#polygon(xx,yy2_ofl,col = "grey95")

legend("topright",c("ELMC", "MMCB"),lty=c(1,5),col=c("black","dark red"),lwd=c(1.5,1),cex=1.2)




}
# Add File creation for ease of access 
# Change errors to CVs 


# Plot Results 
FemPlots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:Female Only", startYear=2018, endYear=2028, data = FemOnly)
MaleModSurvBMPlots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:Males 15% ModSurvBiomass", startYear=2018, endYear=2028, data = MaleOnlyModSurvBM)
StatusQuoPlots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:Status Quo", startYear=2018, endYear=2028, data = StatusQuo)
MaleOnlyR1Plots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:Male Only 10% exp Cap", startYear=2018, endYear=2028, data = MaleOnlyR1)
MaleOnlyR2Plots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:Male Only 15% exp Cap", startYear=2018, endYear=2028, data = MaleOnlyR2)
MaleOnlyR3Plots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:Male Only 20% exp Cap", startYear=2018, endYear=2028, data = MaleOnlyR3)
ABC<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:TAC=ABC", startYear=2018, endYear=2028, data = ABCdat)
DimmerPlots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:TAC=ABC", startYear=2018, endYear=2028, data = Dimmerdat)
BlockPlots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:TAC=ABC", startYear=2018, endYear=2028, data = Blockdat)
ELM30Plots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:TAC=ABC", startYear=2018, endYear=2028, data = ELM30dat)
ELM40Plots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:TAC=ABC", startYear=2018, endYear=2028, data = ELM40dat)
ELM50Plots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:ELM50%", startYear=2018, endYear=2038, data = ELM50)
SurvBMPlots<-plotMSEresults(modelType = "EstMod", HCRname = "Scenario:TAC=ABC", startYear=2018, endYear=2028, data = MaleSurveyBMdat)

test2plots<-plotMSEresults(modelType = "EstMod", HCRname = "MaleR2EX", startYear=2018, endYear=2023, data = test)


GetTableVals<-function(data){
  
  data<-data
  
  TACset<-data$`TAC Set`
  TAC<-data$TAC
  OFL<-data$OFL
  MMB<-data$MMB
  MFB<-data$MFB
  ELMB<-data$ELMB
  ELMB_State<-data$ELMB_State
  Rec<-data$AveRec
  ELMC<-data$ELMC
  MMCB<-data$MMCB
  IMCB<-data$IMCB
  MFDB<-data$MFDB
  IFDB<-data$IFDB
  ELMD<-data$ELMD
  MMDB<-data$MMDB
  IMDB<-data$IMDB
  
  #Extract Averages and standard deviations
  
  #TAC_Set
  ClosedYears<-apply(TACset,1,sum)
  aveClosedYears<-mean(ClosedYears)
  simNOTAC<-rowMeans(TACset)
  aveNOTAC<-mean(simNOTAC)
  sdNOTAC<-sd(simNOTAC)
  cvNOTAC<-sdNOTAC/aveNOTAC
  
  #TAC for all Years 
  simTacAve<-rowMeans(TAC)
  aveTAC<-mean(simTacAve)
  sdTAC<-sd(simTacAve)
  cvTAC<-sdTAC/aveTAC
  
  
  #TAC by years 
  #aveTAC<-colMeans(TAC)
  #sdTAC<-apply(TAC,2,sd)
  #cvTAC<-sdTAC/aveTAC
  maxTAC<-max(TAC)
  minTAC<-min(TAC)
  
  # OFL for all Years 
  simOFLAve<-rowMeans(OFL)
  aveOFL<-mean(simOFLAve)
  sdOFL<-sd(simOFLAve)
  cvOFL<-sdOFL/aveOFL
  
  #OFL by years 
  #aveOFL<-colMeans(OFL)
  #sdOFL<-apply(OFL,2,sd)
  #cvOFL<-sdOFL/aveOFL
  
  #MMB
  simMMBAve<-rowMeans(MMB)
  aveMMB<-mean(simMMBAve)
  sdMMB<-sd(simMMBAve)
  cvMMB<-sdMMB/aveMMB
  
  #aveMMB<-colMeans(MMB)
  #sdMMB<-apply(MMB,2,sd)
  #cvMMB<-sdMMB/aveMMB
  
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
  
  #Recruitment
  simRecAve<-rowMeans(Rec)
  aveRec<-mean(simRecAve)
  sdRec<-sd(simRecAve)
  cvRec<-sdRec/aveRec
  
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

  #OutPut
  LimitData<-c("aveOFL"=aveOFL, "cvOFL"= cvOFL, aveClosedYears, cvNOTAC, aveTAC, minTAC, maxTAC, cvTAC)
  BioData<-c(aveMMB, cvMMB, aveMFB, cvMFB, aveELMB, cvELMB, aveRec, cvRec)
  CatchData<-c(aveELMC,cvELMC, aveMMCB, cvMMCB, aveIMCB, cvIMCB)
  DiscardData<-c(aveELMD, cvELMD, aveMMDB, cvMMDB, aveIMDB, cvIMDB, aveMFDB, cvMFDB, aveIFDB, cvIMDB)
  
  Report<-list("LimitData"=LimitData,"BioData"=BioData,"CatchData"=CatchData,"DiscardData"=DiscardData)
  return(Report)
}


#Viable Rules 
Fem<-GetTableVals(FemOnly)
malR1<-GetTableVals(MaleOnlyR1)
malR2<-GetTableVals(MaleOnlyR2)
malR3<-GetTableVals(MaleOnlyR3)
Dimmer<-GetTableVals(Dimmerdat)
Block<-GetTableVals(Blockdat)
ELM30<-GetTableVals(ELM30dat)
ELM40<-GetTableVals(ELM40dat)
ELM50<-GetTableVals(ELM50dat)
#Test Rules 
ModSurv<-GetTableVals(MaleOnlyModSurvBM)
StatQ<-GetTableVals(StatusQuo)
ABC<-GetTableVals(ABCdat)
Surv<-GetTableVals(MaleSurveyBMdat)

test<-GetTableVals(test)

### Tables ####
#doTable<-function(data){




ViableHCRs<-c("Female Only", "Male Only 10%", "Male Only 15%", "Male Only 20%","Dimmer", "Block", "ELM 30%", "ELM 40%", "ELM 50%")
SenHCRs<-c("Model Survey", "Survey", "ABC", "Status Quo")
Limits<-matrix(data= NA, nrow=7, ncol = 9)
row.names(Limits)<-c("OFL","OFLsd","Closure Years", "TACave", "TAClow", "TAChigh","TACsd")
col.names(Limits)<-ViableHCRs

Limits[1,1]<-  # Female OFL Ave 

BioDat<-matrix


##### Much Later, show sparklines of average TACS, MMB, MFB, ELMB for each candidate HCR
plot(1, type="n",axes=F,ann=F, 
     xlim=c(startYear,endYear),ylim=c(0,9), #ylim is 0 to the number of sparklines you want to create
     yaxs="i", xaxs="i")  

for(i in 2:12){
  lines(x=startYear:endYear, y=(i-2) +
          HS[,i]/max(HS[,i])) #standardizes each column so they're relative sizes 
  par(xpd=NA)
  text(x=startYear, y=i-1.5, names(HS)[i], pos =4)
  par(xpd=F)
  abline(h=i-2, lty = 2, col = "grey")
} 
