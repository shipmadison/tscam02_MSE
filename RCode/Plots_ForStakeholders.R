# Data Plan

require(matrixStats)
require(Rgb)
library("expss")


#Directory Set Up 
###################
dirHCR1<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR1_PROBLEM"
dirHCR2_1<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR2_R1"
dirHCR2_2<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR2_R2"
dirHCR2_3<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR2_R3"
dirHCR2_MS<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR2_MS"
dirHCR2_Surv<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR2_Surv"
dirHCR3<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR3"
dirHCR4<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR4"
dirHCR5<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR5"
dirHCR6_30<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR6_30"
dirHCR6_40<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR6_40"
dirHCR6_50<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR6_50"
dirHCR7<-"C:\\MSE_Runs\\EC2_results\\Run1_randomRecruits\\HCR7"


#Get Data 
#############################################################
#HCR1
setwd(dirHCR1)
HCR1_EstMod<-readRDS("HCR1_EstMod")

#HCR2
setwd(dirHCR2_1)
HCR2_1_EstMod<-readRDS("HCR2_R1_EstMod")
HCR2_1_OpMod<-readRDS("HCR2_R1_OpMod")

setwd(dirHCR2_2)
HCR2_2_EstMod<-readRDS("HCR2_R2_EstMod")
HCR2_2_OpMod<-readRDS("HCR2_R2_OpMod")

setwd(dirHCR2_3)
HCR2_3_EstMod<-readRDS("HCR2_R3_EstMod")
HCR2_3_OpMod<-readRDS("HCR2_R3_OpMod")

#Model Surv
setwd(dirHCR2_MS)
HCR2_MS_EstMod<-readRDS("HCR2_ms_EstMod")
HCR2_MS_OpMod<-readRDS("HCR2_ms_OpMod")

#SurvBio
setwd(dirHCR2_Surv)
HCR2_Surv_EstMod<-readRDS("HCR2_Surv_EstMod")
HCR2_Surv_OpMod<-readRDS("HCR2_Surv_OpMod")

#ABC
setwd(dirHCR3)
HCR3_EstMod<-readRDS("HCR3_EstMod")
HCR3_OpMod<-readRDS("HCR3_OpMod")

#Dimmer
setwd(dirHCR4)
HCR4_EstMod<-readRDS("HCR4_EstMod")
HCR4_OpMod<-readRDS("HCR4_OpMod")

#Block
setwd(dirHCR5)
HCR5_EstMod<-readRDS("HCR5_EstMod")
HCR5_OpMod<-readRDS("HCR5_OpMod")

#ELM30
setwd(dirHCR6_30)
HCR6_30_EstMod<-readRDS("HCR6_30_EstMod")
HCR6_30_OpMod<-readRDS("HCR6_30_OpMod")

#EML40
setwd(dirHCR6_40)
HCR6_40_EstMod<-readRDS("HCR6_40_EstMod")
HCR6_40_OpMod<-readRDS("HCR6_40_OpMod")

#ELM50
setwd(dirHCR6_50)
HCR6_50_EstMod<-readRDS("HCR6_50_EstMod")
HCR6_50_OpMod<-readRDS("HCR6_50_OpMod")

#ELM50
setwd(dirHCR6_50)
HCR6_50_EstMod<-readRDS("HCR6_50_EstMod")
HCR6_50_OpMod<-readRDS("HCR6_50_OpMod")

#StatusQuo
setwd(dirHCR7)
HCR7_EstMod<-readRDS("HCR7_EstMod")
HCR7_OpMod<-readRDS("HCR7_OpMod")

# Data est 
data1e<-HCR1_EstMod
data21e<-HCR2_1_EstMod
data22e<-HCR2_2_EstMod
data23e<-HCR2_3_EstMod
data24e<-HCR2_MS_EstMod
data25e<-HCR2_Surv_EstMod
data3e<-HCR3_EstMod
data4e<-HCR4_EstMod
data5e<-HCR5_EstMod
data61e<-HCR6_30_EstMod
data62e<-HCR6_40_EstMod
data63e<-HCR6_50_EstMod
data7e<-HCR7_EstMod


# Data OP 
data21o<-HCR2_1_OpMod
data22o<-HCR2_2_OpMod
data23o<-HCR2_3_OpMod
data24o<-HCR2_MS_OpMod
data25o<-HCR2_Surv_OpMod
data3o<-HCR3_OpMod
data4o<-HCR4_OpMod
data5o<-HCR5_OpMod
data61o<-HCR6_30_OpMod
data62o<-HCR6_40_OpMod
data63o<-HCR6_50_OpMod
data7o<-HCR7_OpMod


# ######## RE DO THE SPAGHTTII Plot Function


doSpaghettiPlots<-function(data, startYear, endYear, modTyp, sims, perfMet,  title = "title"){
  
  
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
    lines(years, OFL[i,])
  
  
  
}
# 
# # Spaghetti Plots TAC
# #########################
# #Group 1 --> Single Sex Females 
# par(mfrow = c(3, 1))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR1_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Female Only -- TAC")
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC")
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Status Quo -- TAC")
# 
# #Group 2 --> Single Sex Males 
# par(mfrow = c(5, 1))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR2_1_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 10% -- TAC")
# doSpaghettiPlots(HCR2_2_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 15% -- TAC")
# doSpaghettiPlots(HCR2_3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 20% -- TAC")
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC" )
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Status Quo -- TAC")
# 
# #Group 2 -2 --> Sub Group
# par(mfrow = c(3, 1))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR2_2_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 15% -- TAC")
# doSpaghettiPlots(HCR2_MS_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Model Survey Estimates -- TAC")
# doSpaghettiPlots(HCR2_Surv_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Survey Estimates -- TAC")
# 
# 
# # Group 3 --> Single Sex ELM
# par(mfrow = c(5, 1))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR6_30_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 30% -- TAC")
# doSpaghettiPlots(HCR6_40_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 40% -- TAC")
# doSpaghettiPlots(HCR6_50_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 50% -- TAC")
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC" )
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Status Quo -- TAC")
# 
# #################################
# 
# # Spaghetti Plots TAC vs Catches 
# #########################
# #Group 1 --> Single Sex Females 
# par(mfrow = c(3, 2))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 2, 0))
# 
# doSpaghettiPlots(HCR1_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Female Only -- TAC")
# doSpaghettiPlots(HCR1_EstMod, 2018, 2118, modTyp =1, 100, 8, title = "Female Only -- Catches")
# 
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC")
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "TAC = ABC -- Catches")
# 
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7,title = "Status Quo -- TAC")
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Status Quo -- Catches")
# 
# #Group 2 --> Single Sex Males 
# par(mfrow = c(5, 2))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 2, 0))
# 
# doSpaghettiPlots(HCR2_1_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 10% -- TAC")
# doSpaghettiPlots(HCR2_1_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Male 10% -- Catches")
# 
# doSpaghettiPlots(HCR2_2_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 15% -- TAC")
# doSpaghettiPlots(HCR2_2_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Male 15% -- Catches")
# 
# doSpaghettiPlots(HCR2_3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 20% -- TAC")
# doSpaghettiPlots(HCR2_3_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Male 20% -- Catches")
# 
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC" )
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "TAC = ABC -- Catches" )
# 
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Status Quo -- TAC")
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Status Quo -- Catches")
# 
# 
# #Group 2 -2 --> Sub Group
# par(mfrow = c(3, 2))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 2, 0))
# 
# doSpaghettiPlots(HCR2_2_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 15% -- TAC")
# doSpaghettiPlots(HCR2_2_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Male 15% -- Catches")
# 
# doSpaghettiPlots(HCR2_MS_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Model Survey Estimates -- TAC")
# doSpaghettiPlots(HCR2_MS_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Model Survey Estimates -- Catches")
# 
# 
# doSpaghettiPlots(HCR2_Surv_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Survey Estimates -- TAC")
# doSpaghettiPlots(HCR2_Surv_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Survey Estimates -- Catches")
# 
# 
# # Group 3 --> Single Sex ELM
# par(mfrow = c(5, 2))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR6_30_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 30% -- TAC")
# doSpaghettiPlots(HCR6_30_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "ELM 30% -- Catches")
# 
# doSpaghettiPlots(HCR6_40_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 40% -- TAC")
# doSpaghettiPlots(HCR6_40_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "ELM 30% -- Catches")
# 
# doSpaghettiPlots(HCR6_50_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 50% -- TAC")
# doSpaghettiPlots(HCR6_30_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "ELM 30% -- Catches")
# 
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC" )
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "TAC = ABC -- Catches" )
# 
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Status Quo -- TAC")
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 8, title = "Status Quo -- Catches")
# 
# #################################
# 
# # Spaghetti Plots TAC vs Discards 
# #########################
# #Group 1 --> Single Sex Females 
# par(mfrow = c(3, 2))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 2, 0))
# 
# doSpaghettiPlots(HCR1_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Female Only -- TAC")
# doSpaghettiPlots(HCR1_EstMod, 2018, 2118, modTyp =1, 100, 9, title = "Female Only -- Discards")
# 
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC")
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "TAC = ABC -- Discards")
# 
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7,title = "Status Quo -- TAC")
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Status Quo -- Discards")
# 
# #Group 2 --> Single Sex Males 
# par(mfrow = c(5, 2))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 2, 0))
# 
# doSpaghettiPlots(HCR2_1_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 10% -- TAC")
# doSpaghettiPlots(HCR2_1_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Male 10% -- Discards")
# 
# doSpaghettiPlots(HCR2_2_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 15% -- TAC")
# doSpaghettiPlots(HCR2_2_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Male 15% -- Discards")
# 
# doSpaghettiPlots(HCR2_3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 20% -- TAC")
# doSpaghettiPlots(HCR2_3_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Male 20% -- Discards")
# 
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC" )
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "TAC = ABC -- Discards" )
# 
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Status Quo -- TAC")
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Status Quo -- Discards")
# 
# 
# #Group 2 -2 --> Sub Group
# par(mfrow = c(3, 2))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 2, 0))
# 
# doSpaghettiPlots(HCR2_2_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Male 15% -- TAC")
# doSpaghettiPlots(HCR2_2_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Male 15% -- Discards")
# 
# doSpaghettiPlots(HCR2_MS_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Model Survey Estimates -- TAC")
# doSpaghettiPlots(HCR2_MS_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Model Survey Estimates -- Discards")
# 
# 
# doSpaghettiPlots(HCR2_Surv_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Survey Estimates -- TAC")
# doSpaghettiPlots(HCR2_Surv_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Survey Estimates -- Discards")
# 
# 
# # Group 3 --> Single Sex ELM
# par(mfrow = c(5, 2))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR6_30_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 30% -- TAC")
# doSpaghettiPlots(HCR6_30_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "ELM 30% -- Discards")
# 
# doSpaghettiPlots(HCR6_40_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 40% -- TAC")
# doSpaghettiPlots(HCR6_40_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "ELM 30% -- Discards")
# 
# doSpaghettiPlots(HCR6_50_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "ELM 50% -- TAC")
# doSpaghettiPlots(HCR6_30_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "ELM 30% -- Discards")
# 
# doSpaghettiPlots(HCR3_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "TAC = ABC -- TAC" )
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "TAC = ABC -- Discards" )
# 
# doSpaghettiPlots(HCR7_EstMod, 2018, 2118, modTyp =1, 100, 7, title = "Status Quo -- TAC")
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 9, title = "Status Quo -- Discards")
# 
# #################################
# 
# # Spaghetti Plots REC
# #####################
# #Group 1 --> Single Sex Females 
# par(mfrow = c(3, 1))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR1_EstMod, 2018, 2118, modTyp =1, 100, 3, title = "Female Only -- Rec")
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "TAC = ABC -- Rec")
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Status Quo -- Rec")
# 
# #Group 2 --> Single Sex Males 
# par(mfrow = c(5, 1))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR2_1_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Male 10% -- Rec")
# doSpaghettiPlots(HCR2_2_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Male 15% -- Rec")
# doSpaghettiPlots(HCR2_3_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Male 20% -- Rec")
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "TAC = ABC -- Rec" )
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Status Quo -- Rec")
# 
# #Group 2 -2 --> Sub Group
# par(mfrow = c(3, 1))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR2_2_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Male 15% -- Rec")
# doSpaghettiPlots(HCR2_MS_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Model Survey Estimates -- Rec")
# doSpaghettiPlots(HCR2_Surv_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Survey Estimates -- Rec")
# 
# 
# # Group 3 --> Single Sex ELM
# par(mfrow = c(5, 1))
# par(mar = c(2, 4, 2, 1), oma = c(2, 2, 1, 0))
# 
# doSpaghettiPlots(HCR6_30_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "ELM 30% -- Rec")
# doSpaghettiPlots(HCR6_40_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "ELM 40% -- Rec")
# doSpaghettiPlots(HCR6_50_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "ELM 50% -- Rec")
# doSpaghettiPlots(HCR3_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "TAC = ABC -- Rec" )
# doSpaghettiPlots(HCR7_OpMod, 2018, 2118, modTyp =2, 100, 3, title = "Status Quo -- Rec")
# 
# ###########################################################

  
# Trade Off Plots 
################################################
#trade Off
# Trades<-function(data1, data2, data3, data4, data5, sims, modTyp, startYear, endYear, 
#                  lab1 = "lab1", 
#                  lab2 = "lab2", 
#                  lab3 = "lab3", 
#                  lab4 = "lab4", 
#                  lab5= "lab5",
#                  PlotType){
  
##### Years 
 yearsOp<-seq(2018, 2117)
 yearsEst<-seq(2019, 2118)
 
 aveMMB <-5.1764925267e+001 
 closeMMB<-0.25*aveMMB
 aveMFB <-2.3514569603e+001
 
  
#   #Data Pull
# #####################
#  
#   
#   TACSET1<-data1e$`TAC Set`
#   TACSET21<-data21e$`TAC Set`
#   TACSET22<-data22e$`TAC Set`
#   TACSET23<-data23e$`TAC Set`
#   TACSET24<-data24e$`TAC Set`
#   TACSET25<-data25e$`TAC Set`
#   TACSET3<-data3e$`TAC Set`
#   TACSET4<-data4e$`TAC Set`
#   TACSET5<-data5e$`TAC Set`
#   TACSET61<-data61e$`TAC Set`
#   TACSET62<-data62e$`TAC Set`
#   TACSET63<-data63e$`TAC Set`
#   TACSET7<-data7e$`TAC Set`
#   
#   TAC1<-data1e$TAC
#   TAC21<-data21e$TAC
#   TAC22<-data22e$TAC
#   TAC23<-data23e$TAC
#   TAC24<-data24e$TAC
#   TAC25<-data25e$TAC
#   TAC3<-data3e$TAC
#   TAC4<-data4e$TAC
#   TAC5<-data5e$TAC
#   TAC61<-data61e$TAC
#   TAC62<-data62e$TAC
#   TAC63<-data63e$TAC
#   TAC7<-data7e$TAC
#   
#   
#   OFL1<-data1e$OFL
#   OFL21<-data21o$OFL
#   OFL22<-data22o$OFL
#   OFL23<-data23o$OFL
#   OFL24<-data24o$OFL
#   OFL25<-data25o$OFL
#   OFL3<-data3o$OFL
#   OFL4<-data4o$OFL
#   OFL5<-data5o$OFL
#   OFL61<-data61o$OFL
#   OFL62<-data62o$OFL
#   OFL63<-data63o$OFL
#   OFL7<-data7o$OFL 
#  
#   
#   B01<-data1e$B0
#   B021<-data21o$B0
#   B022<-data22o$B0
#   B023<-data23o$B0
#   B024<-data24o$B0
#   B025<-data25o$B0
#   B03<-data3o$B0
#   B04<-data4o$B0
#   B05<-data5o$B0
#   B061<-data61o$B0
#   B062<-data62o$B0
#   B063<-data63o$B0
#   B07<-data7o$B0 
#   
#   Bmsy1<-data1e$Bmsy
#   Bmsy21<-data21o$Bmsy
#   Bmsy22<-data22o$Bmsy
#   Bmsy23<-data23o$Bmsy
#   Bmsy24<-data24o$Bmsy
#   Bmsy25<-data25o$Bmsy
#   Bmsy3<-data3o$Bmsy
#   Bmsy4<-data4o$Bmsy
#   Bmsy5<-data5o$Bmsy
#   Bmsy61<-data61o$Bmsy
#   Bmsy62<-data62o$Bmsy
#   Bmsy63<-data63o$Bmsy
#   Bmsy7<-data7o$Bmsy 
#   
#   MMB1<-data1e$MMB
#   MMB21<-data21o$MMB
#   MMB22<-data22o$MMB
#   MMB23<-data23o$MMB
#   MMB24<-data24o$MMB
#   MMB25<-data25o$MMB
#   MMB3<-data3o$MMB
#   MMB4<-data4o$MMB
#   MMB5<-data5o$MMB
#   MMB61<-data61o$MMB
#   MMB62<-data62o$MMB
#   MMB63<-data63o$MMB
#   MMB7<-data7o$MMB 
#   
#   MFB1<-data1e$MFB
#   MFB21<-data21o$MFB
#   MFB22<-data22o$MFB
#   MFB23<-data23o$MFB
#   MFB24<-data24o$MFB
#   MFB25<-data25o$MFB
#   MFB3<-data3o$MFB
#   MFB4<-data4o$MFB
#   MFB5<-data5o$MFB
#   MFB61<-data61o$MFB
#   MFB62<-data62o$MFB
#   MFB63<-data63o$MFB
#   MFB7<-data7o$MFB 
#   
#   ELMB1<-data1e$ELMB
#   ELMB21<-data21o$ELMB
#   ELMB22<-data22o$ELMB
#   ELMB23<-data23o$ELMB
#   ELMB24<-data24o$ELMB
#   ELMB25<-data25o$ELMB
#   ELMB4<-data4o$ELMB
#   ELMB5<-data5o$ELMB
#   ELMB61<-data61o$ELMB
#   ELMB62<-data62o$ELMB
#   ELMB63<-data63o$ELMB
#   ELMB7<-data7o$ELMB 
#   
#   ELMBs1<-data1e$ELMB_State
#   ELMBs21<-data21o$ELMB_State
#   ELMBs22<-data22o$ELMB_State
#   ELMBs23<-data23o$ELMB_State
#   ELMBs24<-data24o$ELMB_State
#   ELMBs25<-data25o$ELMB_State
#   ELMBs4<-data4o$ELMB_State
#   ELMBs5<-data5o$ELMB_State
#   ELMBs61<-data61o$ELMB_State
#   ELMBs62<-data62o$ELMB_State
#   ELMBs63<-data63o$ELMB_State
#   ELMBs7<-data7o$ELMB_State 
#   
#   Rec1<-data1e$Rec
#   Rec21<-data21o$Rec
#   Rec22<-data22o$Rec
#   Rec23<-data23o$Rec
#   Rec24<-data24o$Rec
#   Rec25<-data25o$Rec
#   Rec3<-data3o$Rec
#   Rec4<-data4o$Rec
#   Rec5<-data5o$Rec
#   Rec61<-data61o$Rec
#   Rec62<-data62o$Rec
#   Rec63<-data63o$Rec
#   RecB7<-data7o$Rec 
#   
# ##################################
  
# More breakdown 
#Extract Averages and standard deviations
#########################  
  
PlotNums<-function(dataOp, dataEst, perfMets){
    
    
    data = dataOp
    
    dataE = dataEst
    
    
    #Closed Years 
    
    if(perfMets == 1){
    
    ClosedYears<-apply(dataE$'TAC Set',1,sum)
    aveNOTAC<-mean(dataE$'TAC Set')
    simNOTAC<-rowMeans(dataE$'TAC Set')
    NoTACMed<-median(dataE$'TAC Set')
    NoTACquants<-quantile(dataE$'TAC Set', prob = c(.05, .25, .75, .95))
    sdNOTAC<-sd(simNOTAC)
    cvNOTAC<-sdNOTAC/aveNOTAC
    
    info <- list("ClosedYears" = ClosedYears,
                 "aveNoTAC"=aveNOTAC, 
                 "simNOTAC"=simNOTAC,
                 "NoTACMed"=NoTACMed,
                 "NoTACquants" = NoTACquants,
                 "sdNOTAC"=sdNOTAC,
                 "cvNOTAC"=cvNOTAC)
   
    }
   
     #TAC
    if( perfMets==2){
    
    SumTAC<-apply(dataE$TAC, 1, sum)
    aveTAC<-mean(dataE$TAC)
    simTacAve<-rowMeans(dataE$TAC)
    TACMed<-median(dataE$TAC)
    TACQuant<-quantile(dataE$TAC, probs =c(.05, .25, .5, .75, .95))
    sdTAC<-sd(simTacAve)
    cvTAC<-sdTAC/aveTAC
    
    maxTAC<-max(dataE$TAC)
    minTAC<-min(dataE$TAC)
    
    info <- list("SumTAC" = SumTAC,
                 "aveTAC"=aveTAC,
                 "simTacAve" = simTacAve,
                 "TACMed" = TACMed,
                 "TACQuant"=TACQuant,
                 "sdTAC"=sdTAC,
                 "cvTAC"=cvTAC,
                 "maxTAC" = maxTAC,
                 "minTAC" = minTAC)
    
    }
    
    #OFL 
    if(perfMets ==3){
   
    simOFLAve<-rowMeans(data$OFL)
    aveOFL<-mean(data$OFL)
    OFLMed<-median(data$OFL)
    OFLQuant<-quantile(data$OFL, probs =c(.05, .25, .5, .75, .95))
    sdOFL<-sd(simOFLAve)
    cvOFL<-sdOFL/aveOFL
    
    info <- list("simOFLAve" = simOFLAve,
                 "aveOFL"=aveOFL,
                 "OFLMed" = OFLMed,
                 "OFLQuant"=OFLQuant,
                 "sdOFL"=sdOFL,
                 "cvOFL"=cvOFL)
    
    }
    
    #B0
    if(perfMets ==4){
    simB0Ave<-rowMeans(data$B0)
    aveB0<-mean(data$B0)
    B0Med<-median(data$B0)
    B0Quant<-quantile(data$B0, probs =c(.05, .25, .5, .75, .95))
    sdB0<-sd(simB0Ave1)
    cvB0<-sdB0/aveB0
    
    info <- list("simB0Ave" = simB0Ave,
                 "aveB0"=aveB0,
                 "B0Med" =B0Med,
                 "B0Quant"=B0Quant,
                 "sdB0"=sdB0,
                 "cvB0"=cvB0)
    }
    
    #Bmsy
    if(perfMets ==4){
    simBmsyAve<-rowMeans(data$Bmsy)
    aveBmsy<-mean(data$Bmsy)
    BmsyMed<-median(data$Bmsy)
    BmsyQuant<-quantile(data$Bmsy, probs =c(.05, .25, .5, .75, .95))
    sdBmsy<-sd(simBmsyAve)
    cvBmsy<-sdBmsy/aveBmsy
    
    info <- list("simBmsyAve" = simBmsyAve,
                 "aveBmsy"=aveBmsy,
                 "BmsyMed" =BmsyMed,
                 "BmsyQuant"=BmsyQuant,
                 "sdBmsy"=sdBmsy,
                 "cvBmsy"=cvBmsy)
    
    }
    
    
    #MMB
    if(perfMets ==5){
    simMMBAve<-rowMeans(data$MMB)
    aveMMB<-mean(data$MMB)
    MMBMed<-median(data$MMB)
    MMBQuant<-quantile(data$MMB, probs =c(.05, .25, .5, .75, .95))
    sdMMB<-sd(simMMBAve)
    cvMMB<-sdMMB/aveMMB
    
    info <- list("simMMBAve" = simMMBAve,
                 "aveMMB"=aveMMB,
                 "MMBMed" =MMBMed,
                 "MMBQuant"=MMBQuant,
                 "sdMMB"=sdMMB,
                 "cvMMB"=cvMMB)
    
    }
    
    
    #MFB
    if(perfMets ==6){
    simMFBAve<-rowMeans(data$MFB)
    aveMFB<-mean(data$MFB)
    MFBMed<-median(data$MFB)
    MFBQuant<-quantile(data$MFB, probs =c(.05, .25, .5, .75, .95))
    sdMFB<-sd(simMFBAve)
    cvMFB<-sdMFB/aveMFB
    
    info <- list("simMFBAve" = simMFBAve,
                 "aveMFB"=aveMFB,
                 "MFBMed" =MFBMed,
                 "MFBQuant"=MFBQuant,
                 "sdMFB"=sdMFB,
                 "cvMMB"=cvMMB)
    }
    
    #ELMB_State
    if(perfMets ==7){
    simELMBAve<-rowMeans(data$ELMB_State)
    aveELMB<-mean(data$ELMB_State)
    ELMBMed<-median(data$ELMB_State)
    ELMBQuant<-quantile(data$ELMB_State, probs =c(.05, .25, .5, .75, .95))
    sdELMB<-sd(simELMBAve)
    cvELMB<-sdELMB/aveELMB
    
    info <- list("simELMBAve" = simELMBAve,
                 "aveELMB"=aveELMB,
                 "ELMBMed" =ELMBMed,
                 "ELMBQuant"=ELMBQuant,
                 "sdELMB"=sdELMB,
                 "cvELMB"=cvELMB)
    }
    
    #Recruitment
    if(perfMets == 8){
    simRecAve<-rowMeans(data$Rec)
    aveRec<-mean(data$Rec)
    RecMed<-median(data$Rec)
    RecQuant<-quantile(data$Rec, probs =c(.05, .25, .5, .75, .95))
    sdRec<-sd(simRecAve)
    cvRec<-sdRec/aveRec
    
    info <- list("simRecAve" = simRecAve,
                 "aveRec"=aveRec,
                 "RecMed" =RecMed,
                 "RecQuant"=RecQuant,
                 "sdRec"=sdRec,
                 "cvRec"=cvRec)
    
    }
    return(info)
  }
  
  
 TACHCR1<-PlotNums(dataOp = data1e, dataEst = data1e, 2)
 TACHCR21<-PlotNums(dataOp = data21o, dataEst = data21e, 2)
 TACHCR22<-PlotNums(dataOp = data22o, dataEst = data22e, 2)
 TACHCR23<-PlotNums(dataOp = data23o, dataEst = data23e, 2)
 TACHCR24<-PlotNums(dataOp = data24o, dataEst = data24e, 2)
 TACHCR25<-PlotNums(dataOp = data25o, dataEst = data25e, 2)
 TACHCR3<-PlotNums(dataOp = data3o, dataEst = data3e, 2)
 TACHCR4<-PlotNums(dataOp = data4o, dataEst = data4e, 2)
 TACHCR5<-PlotNums(dataOp = data5o, dataEst = data5e, 2)
 TACHCR61<-PlotNums(dataOp = data61o, dataEst = data61e, 2)
 TACHCR62<-PlotNums(dataOp = data62o, dataEst = data62e, 2)
 TACHCR63<-PlotNums(dataOp = data63o, dataEst = data63e, 2)
 TACHCR7<-PlotNums(dataOp = data7o, dataEst = data7e, 2)
 
 #OFL 
 
 MMBHCR1<-PlotNums(dataOp = data1e, dataEst = data1e, 5)
 MMBHCR21<-PlotNums(dataOp = data21o, dataEst = data21e, 5)
 MMBHCR22<-PlotNums(dataOp = data22o, dataEst = data22e, 5)
 MMBHCR23<-PlotNums(dataOp = data23o, dataEst = data23e, 5)
 MMBHCR24<-PlotNums(dataOp = data24o, dataEst = data24e, 5)
 MMBHCR25<-PlotNums(dataOp = data25o, dataEst = data25e, 5)
 MMBHCR3<-PlotNums(dataOp = data3o, dataEst = data3e, 5)
 MMBHCR4<-PlotNums(dataOp = data4o, dataEst = data4e, 5)
 MMBHCR5<-PlotNums(dataOp = data5o, dataEst = data5e, 5)
 MMBHCR61<-PlotNums(dataOp = data61o, dataEst = data61e, 5)
 MMBHCR62<-PlotNums(dataOp = data62o, dataEst = data62e, 5)
 MMBHCR63<-PlotNums(dataOp = data63o, dataEst = data63e, 5)
 MMBHCR7<-PlotNums(dataOp = data7o, dataEst = data7e, 5)
 
 
 
  
  
### Plot of TAC mins and Maxes and averages 
 
 library(RColorBrewer)
 
 plot(1,type="n",
      xlim=c(0,14),ylim=c(0,200),
      yaxs="i", xaxs="i", xaxt='n') #no gap between 0 and the axis
 
 axis(side=1, at= c(1,2,3,4,5,6,7,8,9,10,11,12,13), labels= c("Fem",
                                                              "Male10", 
                                                              "Male15", 
                                                              "Male20",
                                                              "FemDim",
                                                              "FemBlk",
                                                              "ELM30",
                                                              "ELM40",
                                                              "ELM50",
                                                              "TAC=ABC",
                                                              "Status Quo",
                                                              "ModSurv",
                                                              "Surv"))
 
 mtext(side=3, line=1, "TAC Ranges" )
 mtext(side=1, line=2, "HCR")
 mtext(side=2, line=2, "1000's of Tons")
 
 #abline(v =0.5, lty = 2, col = "grey70")
 #abline(v =0, lty = 2, col = "grey70")
 #data1
 #HCRcol<-colorRampPalette(brewer.pal(8,"Dark2"))
 
 HCRcol<-viridis(26)
 
 #data1
 points(1,TACHCR1$aveTAC, pch = 19,  cex = 2, col = HCRcol[1] )
 arrows(1, TACHCR1$minTAC, 1, TACHCR1$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[1])
 
 points(2,TACHCR21$aveTAC, pch = 19,  cex = 2,col = HCRcol[3] )
 arrows(2, TACHCR21$minTAC, 2, TACHCR21$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[3])
 
 points(3,TACHCR22$aveTAC, pch = 19,  cex = 2,col = HCRcol[5] )
 arrows(3, TACHCR22$minTAC, 3, TACHCR22$maxTAC, length=0.05, angle=90, code=3, col = HCRcol[5])
 
 points(4,TACHCR23$aveTAC, pch = 19,  cex = 2,col = HCRcol[7] )
 arrows(4, TACHCR23$minTAC, 4, TACHCR22$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[7])
 
 points(5,TACHCR4$aveTAC, pch = 19,  cex = 2,col = HCRcol[9] )
 arrows(5, TACHCR4$minTAC, 5, TACHCR4$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[9])

 points(6,TACHCR5$aveTAC, pch = 19,  cex = 2,col = HCRcol[11] )
 arrows(6, TACHCR5$minTAC, 6, TACHCR5$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[11])
 
 points(7,TACHCR61$aveTAC, pch = 19,  cex = 2,col = HCRcol[13])
 arrows(7, TACHCR61$minTAC, 7, TACHCR61$maxTAC, length=0.05, angle=90, code=3, col = HCRcol[13])
 
 points(8,TACHCR62$aveTAC, pch = 19,  cex = 2, col = HCRcol[15] )
 arrows(8, TACHCR62$minTAC, 8, TACHCR62$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[15])
 
 points(9,TACHCR63$aveTAC, pch = 19,  cex = 2,col = HCRcol[17] )
 arrows(9, TACHCR63$minTAC, 9, TACHCR63$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[17])
 
 points(10,TACHCR3$aveTAC, pch = 13,  cex = 2,col = HCRcol[19] )
 arrows(10, TACHCR3$minTAC, 10, TACHCR3$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[19])
 
 points(11,TACHCR7$aveTAC, pch = 13,  cex = 2 )
 arrows(11, TACHCR7$minTAC, 11, TACHCR7$maxTAC, length=0.05, angle=90, code=3)
 
 points(12,TACHCR24$aveTAC, pch = 13,  cex = 2,col = HCRcol[23] )
 arrows(12, TACHCR24$minTAC, 12, TACHCR24$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[23])
 
 points(13,TACHCR25$aveTAC, pch = 13,  cex = 2,col = HCRcol[25] )
 arrows(13, TACHCR25$minTAC, 13, TACHCR25$maxTAC, length=0.05, angle=90, code=3,col = HCRcol[25])
 
 
 
 legend("topleft",c("Fem",
                    "Male10", 
                    "Male15", 
                    "Male20",
                    "FemDim",
                    "FemBlk",
                    "ELM30",
                    "ELM40",
                    "ELM50",
                    "TAC=ABC",
                    "Status Quo",
                    "ModSurv",
                    "Surv"),pch=c(19,19,19,19,19,19,19,19,19,13,13,13,13), col=c(HCRcol[1],
                                                                                 HCRcol[3],
                                                                                 HCRcol[5],
                                                                                 HCRcol[7],
                                                                                 HCRcol[9],
                                                                                 HCRcol[11],
                                                                                 HCRcol[13],
                                                                                 HCRcol[15],
                                                                                 HCRcol[17],
                                                                                 HCRcol[19],
                                                                                 "black",
                                                                                 HCRcol[23],
                                                                                 HCRcol[25]),cex=.5)
 
  
 # TAC Cv vs MMBprob
 ####################### 
 par(mfrow = c(1, 1))
 par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
 plot(1,type="n",
      xlim=c(0,1),ylim=c(0,20),
      yaxs="i", xaxs="i") #no gap between 0 and the axis
 
 mtext(side=3, line=1, "TAC Median by Variability" )
 mtext(side=2, line=2, "TAC Median 1000's of Tons")
 mtext(side=1, line=2, "TAC cv")
 
 
# abline(v =0.5, lty = 2, col = "grey70")
 
 #data1
 points(TACHCR1$cvTAC,TACHCR1$TACMed, pch = 19, col=HCRcol[1], cex = 3 )
 arrows(TACHCR1$cvTAC, TACHCR1$TACQuant[1], TACHCR1$cvTAC, TACHCR1$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[1])
 
 points(TACHCR21$cvTAC,TACHCR21$TACMed, pch = 19, col=HCRcol[3], cex = 3 )
 arrows(TACHCR21$cvTAC, TACHCR21$TACQuant[1], TACHCR21$cvTAC, TACHCR21$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[3])
 
 points(TACHCR22$cvTAC, TACHCR22$TACMed, pch = 19, col=HCRcol[5], cex = 3 )
 arrows(TACHCR22$cvTAC, TACHCR22$TACQuant[1], TACHCR22$cvTAC, TACHCR22$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[5])
 
 points(TACHCR23$cvTAC,TACHCR23$TACMed, pch = 19, col=HCRcol[7], cex = 3 )
 arrows(TACHCR23$cvTAC, TACHCR23$TACQuant[1], TACHCR23$cvTAC, TACHCR23$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[7])
 
 points(TACHCR4$cvTAC,TACHCR4$TACMed, pch = 19, col=HCRcol[9], cex = 3 )
 arrows(TACHCR4$cvTAC, TACHCR4$TACQuant[1], TACHCR4$cvTAC, TACHCR4$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[9])
 
 points(TACHCR5$cvTAC,TACHCR5$TACMed, pch = 19, col=HCRcol[11], cex = 3 )
 arrows(TACHCR5$cvTAC, TACHCR5$TACQuant[1], TACHCR5$cvTAC, TACHCR5$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[11])
 
 points(TACHCR61$cvTAC,TACHCR61$TACMed, pch = 19, col=HCRcol[13], cex = 3 )
 arrows(TACHCR61$cvTAC, TACHCR61$TACQuant[1], TACHCR61$cvTAC, TACHCR61$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[13])
 
 points(TACHCR62$cvTAC,TACHCR62$TACMed, pch = 19, col=HCRcol[15], cex = 3 )
 arrows(TACHCR62$cvTAC, TACHCR62$TACQuant[1], TACHCR62$cvTAC, TACHCR62$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[15])
 
 points(TACHCR63$cvTAC,TACHCR63$TACMed, pch = 19, col=HCRcol[17], cex = 3 )
 arrows(TACHCR63$cvTAC, TACHCR63$TACQuant[1], TACHCR63$cvTAC, TACHCR63$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[17])
 
 points(TACHCR3$cvTAC,TACHCR3$TACMed, pch = 19, col=HCRcol[19], cex = 3 )
 arrows(TACHCR3$cvTAC, TACHCR3$TACQuant[1], TACHCR3$cvTAC, TACHCR3$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[19])
 
 points(TACHCR7$cvTAC,TACHCR7$TACMed, pch = 19, col=HCRcol[21], cex = 3 )
 arrows(TACHCR7$cvTAC, TACHCR7$TACQuant[1], TACHCR7$cvTAC, TACHCR7$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[21])
 
 points(TACHCR24$cvTAC,TACHCR24$TACMed, pch = 19, col=HCRcol[23], cex = 3 )
 arrows(TACHCR24$cvTAC, TACHCR24$TACQuant[1], TACHCR24$cvTAC, TACHCR24$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[23])
 
 points(TACHCR25$cvTAC,TACHCR25$TACMed, pch = 19, col=HCRcol[25], cex = 3 )
 arrows(TACHCR25$cvTAC, TACHCR25$TACQuant[1], TACHCR25$cvTAC, TACHCR25$TACQuant[5], length=0.05, angle=90, code=3,col = HCRcol[25])
 
 
 legend("topleft",c("Fem",
                    "Male10", 
                    "Male15", 
                    "Male20",
                    "FemDim",
                    "FemBlk",
                    "ELM30",
                    "ELM40",
                    "ELM50",
                    "TAC=ABC",
                    "Status Quo",
                    "ModSurv",
                    "Surv"),pch=c(19,19,19,19,19,19,19,19,19,13,13,13,13), col=c(HCRcol[1],
                                                                                 HCRcol[3],
                                                                                 HCRcol[5],
                                                                                 HCRcol[7],
                                                                                 HCRcol[9],
                                                                                 HCRcol[11],
                                                                                 HCRcol[13],
                                                                                 HCRcol[15],
                                                                                 HCRcol[17],
                                                                                 HCRcol[19],
                                                                                 HCRcol[21],
                                                                                 HCRcol[23],
                                                                                 HCRcol[25]),cex=.5)
  
 
 ##############
 
 pastMMB<-read.csv("MMB_past.csv") 
 
 par(mfrow = c(1, 1))
 par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
 plot(1,type="n",
      xlim=c(1983,2118),ylim=c(0,150),
      yaxs="i", xaxs="i") #no gap between 0 and the axis


  lines(pastMMB$years, pastMMB$MMB, lty = 2)
  lines(yearsOp, MMBHCR1$simMMBAve, col = HCRcol[1])
  lines(yearsOp, MMBHCR21$simMMBAve, col = HCRcol[3])
  lines(yearsOp, MMBHCR22$simMMBAve, col = HCRcol[5])
  lines(yearsOp, MMBHCR23$simMMBAve, col = HCRcol[7])
  lines(yearsOp, MMBHCR4$simMMBAve, col = HCRcol[9])
  lines(yearsOp, MMBHCR5$simMMBAve, col = HCRcol[11])
  lines(yearsOp, MMBHCR61$simMMBAve, col = HCRcol[13])
  lines(yearsOp, MMBHCR62$simMMBAve, col = HCRcol[15])
  lines(yearsOp, MMBHCR63$simMMBAve, col = HCRcol[17])
  lines(yearsOp, MMBHCR3$simMMBAve, col = HCRcol[19])
  lines(yearsOp, MMBHCR24$simMMBAve, col = HCRcol[23])
  lines(yearsOp, MMBHCR25$simMMBAve, col = HCRcol[25])
  lines(yearsOp, MMBHCR7$simMMBAve,  lwd = 3)
  
 
 
  
##############################  
#TAC_Set
  #data1
  ClosedYears1<-apply(TACSET1,1,sum)
  aveClosedYears1<-mean(ClosedYears1)
  NoTACMed1<-median(TACSET1)
  NoTACquants1<-quantile(TACSET1, prob = c(.05, .25, .75, .95))
  simNOTAC1<-rowMeans(TACSET1)
  aveNOTAC1<-mean(simNOTAC1)
  sdNOTAC1<-sd(simNOTAC1)
  cvNOTAC1<-sdNOTAC1/aveNOTAC1
  
  #data21
  ClosedYears21<-apply(TACSET21,1,sum)
  aveClosedYears21<-mean(ClosedYears21)
  NoTACMed21<-median(TACSET21)
  NoTACquants21<-quantile(TACSET21, prob = c(.05, .25, .75, .95))
  simNOTAC21<-rowMeans(TACSET21)
  aveNOTAC21<-mean(simNOTAC21)
  sdNOTAC21<-sd(simNOTAC21)
  cvNOTAC21<-sdNOTAC21/aveNOTAC21
  
  #data22
  ClosedYears22<-apply(TACSET22,1,sum)
  aveClosedYears22<-mean(ClosedYears22)
  NoTACMed22<-median(TACSET22)
  NoTACquants22<-quantile(TACSET22, prob = c(.05, .25, .75, .95))
  simNOTAC22<-rowMeans(TACSET22)
  aveNOTAC22<-mean(simNOTAC22)
  sdNOTAC22<-sd(simNOTAC22)
  cvNOTAC22<-sdNOTAC22/aveNOTAC22
  
  #data23
  ClosedYears23<-apply(TACSET23,1,sum)
  aveClosedYears23<-mean(ClosedYears23)
  NoTACMed23<-median(TACSET23)
  NoTACquants23<-quantile(TACSET23, prob = c(.05, .25, .75, .95))
  simNOTAC23<-rowMeans(TACSET23)
  aveNOTAC23<-mean(simNOTAC23)
  sdNOTAC23<-sd(simNOTAC23)
  cvNOTAC23<-sdNOTAC23/aveNOTAC23
  
  #data24
  ClosedYears24<-apply(TACSET24,1,sum)
  aveClosedYears24<-mean(ClosedYears24)
  NoTACMed24<-median(TACSET24)
  NoTACquants24<-quantile(TACSET24, prob = c(.05, .25, .75, .95))
  simNOTAC24<-rowMeans(TACSET24)
  aveNOTAC24<-mean(simNOTAC24)
  sdNOTAC24<-sd(simNOTAC24)
  cvNOTAC24<-sdNOTAC24/aveNOTAC24
  
  #data25
  ClosedYears25<-apply(TACSET25,1,sum)
  aveClosedYears25<-mean(ClosedYears25)
  NoTACMed25<-median(TACSET25)
  NoTACquants25<-quantile(TACSET25, prob = c(.05, .25, .75, .95))
  simNOTAC25<-rowMeans(TACSET25)
  aveNOTAC25<-mean(simNOTAC25)
  sdNOTAC25<-sd(simNOTAC25)
  cvNOTAC25<-sdNOTAC25/aveNOTAC25
  
  #data3
  ClosedYears3<-apply(TACSET3,1,sum)
  aveClosedYears3<-mean(ClosedYears3)
  NoTACMed3<-median(TACSET3)
  NoTACQuants3<-colQuantiles(TACSET3, probs = c(.05, .25, .5, .75, .95))
  simNOTAC3<-rowMeans(TACSET3)
  aveNOTAC3<-mean(simNOTAC3)
  sdNOTAC3<-sd(simNOTAC3)
  cvNOTAC3<-sdNOTAC3/aveNOTAC3
  
  #data4
  ClosedYears4<-apply(TACSET4,1,sum)
  aveClosedYears4<-mean(ClosedYears4)
  NoTACMed4<-median(TACSET4)
  NoTACQuants4<-colQuantiles(TACSET4, probs = c(.05, .25, .5, .75, .95))
  simNOTAC4<-rowMeans(TACSET4)
  aveNOTAC4<-mean(simNOTAC4)
  sdNOTAC4<-sd(simNOTAC4)
  cvNOTAC4<-sdNOTAC4/aveNOTAC4
  
  #data5
  ClosedYears5<-apply(TACSET5,1,sum)
  aveClosedYears5<-mean(ClosedYears5)
  NoTACMed5<-median(TACSET5)
  NoTACQuants5<-colQuantiles(TACSET5, probs = c(.05, .25, .5, .75, .95))
  simNOTAC5<-rowMeans(TACSET5)
  aveNOTAC5<-mean(simNOTAC5)
  sdNOTAC5<-sd(simNOTAC5)
  cvNOTAC5<-sdNOTAC5/aveNOTAC5
  
  #data61
  ClosedYears61<-apply(TACSET61,1,sum)
  aveClosedYears61<-mean(ClosedYears61)
  NoTACMed61<-median(TACSET61)
  NoTACQuants5<-colQuantiles(TACSET61, probs = c(.05, .25, .5, .75, .95))
  simNOTAC61<-rowMeans(TACSET61)
  aveNOTAC61<-mean(simNOTAC61)
  sdNOTAC61<-sd(simNOTAC61)
  cvNOTAC61<-sdNOTAC61/aveNOTAC61
  
  #data62
  ClosedYears62<-apply(TACSET62,1,sum)
  aveClosedYears62<-mean(ClosedYears62)
  NoTACMed62<-median(TACSET62)
  NoTACQuants5<-colQuantiles(TACSET62, probs = c(.05, .25, .5, .75, .95))
  simNOTAC62<-rowMeans(TACSET62)
  aveNOTAC62<-mean(simNOTAC62)
  sdNOTAC62<-sd(simNOTAC62)
  cvNOTAC62<-sdNOTAC62/aveNOTAC62
  
  #data63
  ClosedYears63<-apply(TACSET63,1,sum)
  aveClosedYears63<-mean(ClosedYears63)
  NoTACMed63<-median(TACSET63)
  NoTACQuants5<-colQuantiles(TACSET63, probs = c(.05, .25, .5, .75, .95))
  simNOTAC63<-rowMeans(TACSET63)
  aveNOTAC63<-mean(simNOTAC63)
  sdNOTAC63<-sd(simNOTAC63)
  cvNOTAC63<-sdNOTAC63/aveNOTAC63
  
  #data7
  ClosedYears7<-apply(TACSET7,1,sum)
  aveClosedYears7<-mean(ClosedYears7)
  NoTACMed7<-median(TACSET7)
  NoTACQuants5<-colQuantiles(TACSET7, probs = c(.05, .25, .5, .75, .95))
  simNOTAC7<-rowMeans(TACSET7)
  aveNOTAC7<-mean(simNOTAC7)
  sdNOTAC7<-sd(simNOTAC7)
  cvNOTAC7<-sdNOTAC7/aveNOTAC7
  
  
#TAC 
  
  #data1
  SumTAC1<-apply(TAC1, 1, sum)
  simTacAve1<-rowMeans(TAC1)
  aveTAC1<-mean(simTacAve1)
  TACMed1<-median(TAC1)
  TACQuant1<-quantile(TAC1, probs =c(.05, .25, .5, .75, .95))
  sdTAC1<-sd(simTacAve1)
  cvTAC1<-sdTAC1/aveTAC1
  
  maxTAC1<-max(TAC1)
  minTAC1<-min(TAC1)
  
  #data2_1
  SumTAC21<-apply(TAC21, 1, sum)
  simTacAve21<-rowMeans(TAC21)
  aveTAC21<-mean(simTacAve21)
  TACMed21<-median(TAC21)
  TACQuant21<-quantile(TAC21, probs =c(.05, .25, .5, .75, .95))
  sdTAC21<-sd(simTacAve21)
  cvTAC21<-sdTAC21/aveTAC21
  
  maxTAC21<-max(TAC21)
  minTAC21<-min(TAC21)
  
  #data2_2
  SumTAC22<-apply(TAC22, 1, sum)
  simTacAve22<-rowMeans(TAC22)
  aveTAC22<-mean(simTacAve22)
  TACMed22<-median(TAC22)
  TACQuant22<-quantile(TAC22, probs =c(.05, .25, .5, .75, .95))
  sdTAC22<-sd(simTacAve22)
  cvTAC22<-sdTAC22/aveTAC22
  
  maxTAC22<-max(TAC22)
  minTAC22<-min(TAC22)
  
  #data2_3
  SumTAC23<-apply(TAC23, 1, sum)
  simTacAve23<-rowMeans(TAC23)
  aveTAC23<-mean(simTacAve23)
  TACMed23<-median(TAC23)
  TACQuant23<-quantile(TAC23, probs =c(.05, .25, .5, .75, .95))
  sdTAC23<-sd(simTacAve23)
  cvTAC23<-sdTAC23/aveTAC23
  
  maxTAC23<-max(TAC23)
  minTAC23<-min(TAC23)
  
  #data2_4 (ModSurv)
  SumTAC24<-apply(TAC24, 1, sum)
  simTacAve24<-rowMeans(TAC24)
  aveTAC24<-mean(simTacAve24)
  TACMed24<-median(TAC24)
  TACQuant24<-quantile(TAC24, probs =c(.05, .25, .5, .75, .95))
  sdTAC24<-sd(simTacAve24)
  cvTAC24<-sdTAC24/aveTAC24
  
  maxTAC24<-max(TAC24)
  minTAC24<-min(TAC24)
  
  #data2_5 (Surv)
  SumTAC25<-apply(TAC25, 1, sum)
  simTacAve25<-rowMeans(TAC25)
  aveTAC25<-mean(simTacAve25)
  TACMed25<-median(TAC25)
  TACQuant25<-quantile(TAC25, probs =c(.05, .25, .5, .75, .95))
  sdTAC25<-sd(simTacAve25)
  cvTAC25<-sdTAC25/aveTAC25
  
  maxTAC25<-max(TAC25)
  minTAC25<-min(TAC25)
  
  #data3
  SumTAC3<-apply(TAC3, 1, sum)
  simTacAve3<-rowMeans(TAC3)
  aveTAC3<-mean(simTacAve3)
  TACMed3<-median(TAC3)
  TACQuant3<-quantile(TAC3, probs =c(.05, .25, .5, .75, .95))
  sdTAC3<-sd(simTacAve3)
  cvTAC3<-sdTAC3/aveTAC3
  
  maxTAC3<-max(TAC3)
  minTAC3<-min(TAC3)
  
  #data4
  simTacAve4<-rowMeans(TAC4)
  aveTAC4<-mean(simTacAve4)
  TACMed4<-median(TAC4)
  TACQuant4<-quantile(TAC4, probs=c(.05,.95))
  sdTAC4<-sd(simTacAve4)
  cvTAC4<-sdTAC4/aveTAC4
  
  maxTAC4<-max(TAC4)
  minTAC4<-min(TAC4)
  
  #data5
  simTacAve5<-rowMeans(TAC5)
  aveTAC5<-mean(simTacAve5)
  TACMed5<-median(TAC5)
  TACQuant5<-quantile(TAC5, probs=c(.05,.95))
  sdTAC5<-sd(simTacAve5)
  cvTAC5<-sdTAC5/aveTAC5
  
  maxTAC5<-max(TAC5)
  minTAC5<-min(TAC5)
  
  #data61
  simTacAve5<-rowMeans(TAC5)
  aveTAC5<-mean(simTacAve5)
  TACMed5<-median(TAC5)
  TACQuant5<-quantile(TAC5, probs=c(.05,.95))
  sdTAC5<-sd(simTacAve5)
  cvTAC5<-sdTAC5/aveTAC5
  
  maxTAC5<-max(TAC5)
  minTAC5<-min(TAC5)
  
  #data62
  simTacAve5<-rowMeans(TAC5)
  aveTAC5<-mean(simTacAve5)
  TACMed5<-median(TAC5)
  TACQuant5<-quantile(TAC5, probs=c(.05,.95))
  sdTAC5<-sd(simTacAve5)
  cvTAC5<-sdTAC5/aveTAC5
  
  maxTAC5<-max(TAC5)
  minTAC5<-min(TAC5)
  
  
#OFL 
  #Data1
  simOFLAve1<-rowMeans(OFL1)
  aveOFL1<-mean(simOFLAve1)
  OFLMed1<-median(OFL1)
  OFLQuant1<-quantile(OFL1, probs=c(.05,.95))
  sdOFL1<-sd(simOFLAve1)
  cvOFL1<-sdOFL1/aveOFL1
  
  #Data2
  simOFLAve2<-rowMeans(OFL2)
  aveOFL2<-mean(simOFLAve2)
  OFLMed2<-median(OFL2)
  OFLQuant2<-quantile(OFL2, probs=c(.05,.95))
  sdOFL2<-sd(simOFLAve2)
  cvOFL2<-sdOFL2/aveOFL2
  
  #Data3
  simOFLAve3<-rowMeans(OFL3)
  aveOFL3<-mean(simOFLAve3)
  OFLMed3<-median(OFL3)
  OFLQuant3<-quantile(OFL3, probs=c(.05,.95))
  sdOFL3<-sd(simOFLAve3)
  cvOFL3<-sdOFL3/aveOFL3
  
  #Data4
  simOFLAve4<-rowMeans(OFL4)
  aveOFL4<-mean(simOFLAve4)
  OFLMed4<-median(OFL4)
  OFLQuant4<-quantile(OFL4, probs=c(.05,.95))
  sdOFL4<-sd(simOFLAve4)
  cvOFL4<-sdOFL4/aveOFL4
  
  #Data5
  simOFLAve5<-rowMeans(OFL5)
  aveOFL5<-mean(simOFLAve5)
  OFLMed5<-median(OFL5)
  OFLQuant5<-quantile(OFL5, probs=c(.05,.95))
  sdOFL5<-sd(simOFLAve5)
  cvOFL5<-sdOFL5/aveOFL5
  
#B0
  #Data1
  simB0Ave1<-rowMeans(B01)
  aveB01<-mean(simB0Ave1)
  B0Med1<-median(B01)
  B0Quant1<-quantile(B01, probs=c(.05,.95))
  sdB01<-sd(simB0Ave1)
  cvB01<-sdB01/aveB01
  
  #Data2
  simB0Ave2<-rowMeans(B02)
  aveB02<-mean(simB0Ave2)
  B0Med2<-median(B02)
  B0Quant2<-quantile(B02, probs=c(.05,.95))
  sdB02<-sd(simB0Ave2)
  cvB02<-sdB02/aveB02
  
  #Data3
  simB0Ave3<-rowMeans(B03)
  aveB03<-mean(simB0Ave3)
  B0Med3<-median(B03)
  B0Quant3<-quantile(B03, probs=c(.05,.95))
  sdB03<-sd(simB0Ave3)
  cvB03<-sdB03/aveB03
  
  #Data4
  simB0Ave4<-rowMeans(B04)
  aveB04<-mean(simB0Ave4)
  B0Med4<-median(B04)
  B0Quant4<-quantile(B04, probs=c(.05,.95))
  sdB04<-sd(simB0Ave4)
  cvB04<-sdB04/aveB04
  
  #Data5
  simB0Ave5<-rowMeans(B05)
  aveB05<-mean(simB0Ave5)
  B0Med5<-median(B05)
  B0Quant5<-quantile(B05, probs=c(.05,.95))
  sdB05<-sd(simB0Ave5)
  cvB05<-sdB05/aveB05
  
#BMSY
  #Data1
  simBmsyAve1<-rowMeans(Bmsy1)
  aveBmsy1<-mean(simBmsyAve1)
  BmsyMed1<-median(Bmsy1)
  BmsyQuant1<-quantile(Bmsy1, probs=c(.05,.95))
  sdBmsy1<-sd(simBmsyAve1)
  cvBmsy1<-sdBmsy1/aveBmsy1
  
  #Data2
  simBmsyAve2<-rowMeans(Bmsy2)
  aveBmsy2<-mean(simBmsyAve2)
  BmsyMed2<-median(Bmsy2)
  BmsyQuant2<-quantile(Bmsy2, probs=c(.05,.95))
  sdBmsy2<-sd(simBmsyAve2)
  cvBmsy2<-sdBmsy2/aveBmsy2
  
  #Data3
  simBmsyAve3<-rowMeans(Bmsy3)
  aveBmsy3<-mean(simBmsyAve3)
  BmsyMed3<-median(Bmsy3)
  BmsyQuant3<-quantile(Bmsy3, probs=c(.05,.95))
  sdBmsy3<-sd(simBmsyAve3)
  cvBmsy3<-sdBmsy3/aveBmsy3
  
  #Data4
  simBmsyAve4<-rowMeans(Bmsy4)
  aveBmsy4<-mean(simBmsyAve4)
  BmsyMed4<-median(Bmsy4)
  BmsyQuant4<-quantile(Bmsy4, probs=c(.05,.95))
  sdBmsy4<-sd(simBmsyAve4)
  cvBmsy4<-sdBmsy4/aveBmsy4
  
  #Data5
  simBmsyAve5<-rowMeans(Bmsy5)
  aveBmsy5<-mean(simBmsyAve5)
  BmsyMed5<-median(Bmsy5)
  BmsyQuant5<-quantile(Bmsy5, probs=c(.05,.95))
  sdBmsy5<-sd(simBmsyAve5)
  cvBmsy5<-sdBmsy5/aveBmsy5
  
  
#MMB
  #Data1
  simMMBAve1<-rowMeans(MMB1)
  aveMMB1<-mean(simMMBAve1)
  MMBMed1<-median(MMB1)
  MMBQuant1<-quantile(MMB1, probs=c(.05,.95))
  sdMMB1<-sd(simMMBAve1)
  cvMMB1<-sdMMB1/aveMMB1
  
  #Data2
  simMMBAve2<-rowMeans(MMB2)
  aveMMB2<-mean(simMMBAve2)
  MMBMed2<-median(MMB2)
  MMBQuant2<-quantile(MMB2, probs=c(.05,.95))
  sdMMB2<-sd(simMMBAve2)
  cvMMB2<-sdMMB2/aveMMB2
  
  #Data3
  simMMBAve3<-rowMeans(MMB3)
  aveMMB3<-mean(simMMBAve3)
  MMBMed3<-median(MMB3)
  MMBQuant3<-quantile(MMB3, probs=c(.05,.95))
  sdMMB3<-sd(simMMBAve3)
  cvMMB3<-sdMMB3/aveMMB3
  
  #Data4
  simMMBAve4<-rowMeans(MMB4)
  aveMMB4<-mean(simMMBAve4)
  MMBMed4<-median(MMB4)
  MMBQuant4<-quantile(MMB4, probs=c(.05,.95))
  sdMMB4<-sd(simMMBAve4)
  cvMMB4<-sdMMB4/aveMMB4
  
  #Data5
  simMMBAve5<-rowMeans(MMB5)
  aveMMB5<-mean(simMMBAve5)
  MMBMed5<-median(MMB5)
  MMBQuant5<-quantile(MMB5, probs=c(.05,.95))
  sdMMB5<-sd(simMMBAve5)
  cvMMB5<-sdMMB5/aveMMB5
  
#MFB
  #Data1
  simMFBAve1<-rowMeans(MFB1)
  aveMFB1<-mean(simMFBAve1)
  MFBMed1<-median(MFB1)
  MFBQuant1<-quantile(MFB1, probs=c(.05,.95))
  sdMFB1<-sd(simMFBAve1)
  cvMFB1<-sdMFB1/aveMFB1
  
  #Data2
  simMFBAve2<-rowMeans(MFB2)
  aveMFB2<-mean(simMFBAve2)
  MFBMed2<-median(MFB2)
  MFBQuant2<-quantile(MFB2, probs=c(.05,.95))
  sdMFB2<-sd(simMFBAve2)
  cvMFB2<-sdMFB2/aveMFB2
  
  #Data3
  simMFBAve3<-rowMeans(MFB3)
  aveMFB3<-mean(simMFBAve3)
  MFBMed3<-median(MFB3)
  MFBQuant3<-quantile(MFB3, probs=c(.05,.95))
  sdMFB3<-sd(simMFBAve3)
  cvMFB3<-sdMFB3/aveMFB3
  
  #Data4
  simMFBAve4<-rowMeans(MFB4)
  aveMFB4<-mean(simMFBAve4)
  MFBMed4<-median(MFB4)
  MFBQuant4<-quantile(MFB4, probs=c(.05,.95))
  sdMFB4<-sd(simMFBAve4)
  cvMFB4<-sdMFB4/aveMFB4
  
  #Data5
  simMFBAve5<-rowMeans(MFB5)
  aveMFB5<-mean(simMFBAve5)
  MFBMed5<-median(MFB5)
  MFBQuant5<-quantile(MFB5, probs=c(.05,.95))
  sdMFB5<-sd(simMFBAve5)
  cvMFB5<-sdMFB5/aveMFB5
  
  #ELMB_State
  #data1
  simELMBAve1<-rowMeans(ELMB_S1)
  aveELMB1<-mean(simELMBAve1)
  ELMBMed1<-median(ELMB_S1)
  ELMBQuant1<-quantile(ELMB_S1, probs=c(.05,.95))
  sdELMB1<-sd(simELMBAve1)
  cvELMB1<-sdELMB1/aveELMB1
  
  #data2
  simELMBAve2<-rowMeans(ELMB_S2)
  aveELMB2<-mean(simELMBAve2)
  ELMBMed2<-median(ELMB_S2)
  ELMBQuant2<-quantile(ELMB_S2, probs=c(.05,.95))
  sdELMB2<-sd(simELMBAve2)
  cvELMB2<-sdELMB2/aveELMB2
  
  #data3
  simELMBAve3<-rowMeans(ELMB_S3)
  aveELMB3<-mean(simELMBAve3)
  ELMBMed3<-median(ELMB_S3)
  ELMBQuant3<-quantile(ELMB_S3, probs=c(.05,.95))
  sdELMB3<-sd(simELMBAve3)
  cvELMB3<-sdELMB3/aveELMB3
  
  #data4
  simELMBAve4<-rowMeans(ELMB_S4)
  aveELMB4<-mean(simELMBAve4)
  ELMBMed4<-median(ELMB_S4)
  ELMBQuant4<-quantile(ELMB_S4, probs=c(.05,.95))
  sdELMB4<-sd(simELMBAve4)
  cvELMB4<-sdELMB4/aveELMB4
  
  #data5
  simELMBAve5<-rowMeans(ELMB_S5)
  aveELMB5<-mean(simELMBAve5)
  ELMBMed5<-median(ELMB_S5)
  ELMBQuant5<-quantile(ELMB_S5, probs=c(.05,.95))
  sdELMB5<-sd(simELMBAve5)
  cvELMB5<-sdELMB5/aveELMB5
  
  #Recruitment
  #data1
  simRecAve1<-rowMeans(Rec1)
  aveRec1<-mean(simRecAve1)
  RecMed1<-median(Rec1)
  RecQuant1<-quantile(Rec1, probs=c(.05,.95))
  sdRec1<-sd(simRecAve1)
  cvRec1<-sdRec1/aveRec1
  
  #data2
  simRecAve2<-rowMeans(Rec2)
  aveRec2<-mean(simRecAve2)
  RecMed2<-median(Rec2)
  RecQuant2<-quantile(Rec2, probs=c(.05,.95))
  sdRec2<-sd(simRecAve2)
  cvRec2<-sdRec2/aveRec2
  
  #data3
  simRecAve3<-rowMeans(Rec3)
  aveRec3<-mean(simRecAve3)
  RecMed3<-median(Rec3)
  RecQuant3<-quantile(Rec3, probs=c(.05,.95))
  sdRec3<-sd(simRecAve3)
  cvRec3<-sdRec3/aveRec3
  
  #data4
  simRecAve4<-rowMeans(Rec4)
  aveRec4<-mean(simRecAve4)
  RecMed4<-median(Rec4)
  RecQuant4<-quantile(Rec4, probs=c(.05,.95))
  sdRec4<-sd(simRecAve4)
  cvRec4<-sdRec4/aveRec4
  
  #data5
  simRecAve5<-rowMeans(Rec5)
  aveRec5<-mean(simRecAve5)
  RecMed5<-median(Rec5)
  RecQuant5<-quantile(Rec5, probs=c(.05,.95))
  sdRec5<-sd(simRecAve5)
  cvRec5<-sdRec5/aveRec5
  
  ##
  
  years<-seq(startYear, endYear)
  
  aveMMB <-5.1764925267e+001 
  closeMMB<-0.25*aveMMB
  aveMFB <-2.3514569603e+001
  
#############################  
  
# Proportion of MMB < MMBave
########################################  
  counter<-0
  
#Data1 
  for(i in 1:sims){
    for(j in 1:length(years)){
      #print(MMB1[i,j])
      if(MMB1[i,j] < aveMMB){
        counter<-counter+1
        propMMB1<-counter/length(MMB1)
      }
    }
  }
  
  counter<-0
#Data2 
  for(i in 1:sims){
    for(j in 1:length(years)){
      #print(MMB2[i,j])
      if(MMB2[i,j] < aveMMB){
        counter<-counter+1
        propMMB2<-counter/length(MMB2)
      }
    }
  }
  
  counter<-0
#Data3 
  for(i in 1:sims){
    for(j in 1:length(years)){
      #print(MMB3[i,j])
      if(MMB3[i,j] < aveMMB){
        counter<-counter+1
        propMMB3<-counter/length(MMB3)
      }
    }
  }  
  
  counter<-0
#Data4 
  for(i in 1:sims){
    for(j in 1:length(years)){
      #print(MMB4[i,j])
      if(MMB4[i,j] < aveMMB){
        counter<-counter+1
        propMMB4<-counter/length(MMB4)
      }
    }
  }  
  
  counter<-0
#Data5 
  for(i in 1:sims){
    for(j in 1:length(years)){
      #print(MMB5[i,j])
      if(MMB5[i,j] < aveMMB){
        counter<-counter+1
        propMMB5<-counter/length(MMB5)
      }
    }
  }    
  
##################################
 
# TAC CV vs Prob MMB < 0.25
  
# Proportion of MMB < 0.25MMBave
########################################  
  counter<-0
  
  #Data1 
  for(i in 1:sims){
    for(j in 1:length(years)){
     # print(MMB1[i,j])
      if(MMB1[i,j] < closeMMB){
        counter<-counter+1
        }
      }
     }
  if(counter == 0){
    propclose1 <- 0}else{
      propclose1<-counter/length(MMB1)
  }
  
  counter<-0
  
  #Data2 
  for(i in 1:sims){   
    for(j in 1:length(years)){
      #print(MMB2[i,j])
      if(MMB2[i,j] < closeMMB){
        counter<-counter+1
      }
    }
  }
  if(counter == 0){
    propclose2 <- 0}else{
      propclose2<-counter/length(MMB2)
    }
  
  counter<-0
  
  #Data3 
  for(i in 1:sims){
    for(j in 1:length(years)){
      #print(MMB3[i,j])
      if(MMB3[i,j] < closeMMB){
        counter<-counter+1
      }
    }
  } 
  if(counter == 0){
    propclose3 <- 0}else{
      propclose3<-counter/length(MMB3)
    }
  
  
  counter<-0
  #Data4 
  for(i in 1:sims){
    for(j in 1:length(years)){
      #print(MMB4[i,j])
      if(MMB4[i,j] < closeMMB){
        counter<-counter+1
      }
    }
  }  
  if(counter == 0){
    propclose4 <- 0}else{
      propclose4<-counter/length(MMB4)
    }
  
  counter<-0
  #Data5 
  for(i in 1:sims){
    for(j in 1:length(years)){
      #print(MMB5[i,j])
      if(MMB5[i,j] < closeMMB){
        counter<-counter+1
      }
    }
  }    
  if(counter == 0){
    propclose5 <- 0}else{
      propclose5<-counter/length(MMB5)
    }
########################################  
    
    # TAC CV vs Prob MFB < MFBave
    
########################################  
    counter<-0
    
    #Data1 
    for(i in 1:sims){
      for(j in 1:length(years)){
        #print(MFB1[i,j])
        if(MFB1[i,j] < aveMFB){
          counter<-counter+1
          propMFB1<-counter/length(MFB1)
        }
      }
    }
  
    counter<-0 
    #Data2 
    for(i in 1:sims){
      for(j in 1:length(years)){
        #print(MMB2[i,j])
        if(MFB2[i,j] < aveMFB){
          counter<-counter+1
          propMFB2<-counter/length(MFB2)
        }
      }
    }
    
    counter<-0 
    #Data3 
    for(i in 1:sims){
      for(j in 1:length(years)){
        #print(MMB3[i,j])
        if(MFB3[i,j] < aveMFB){
          counter<-counter+1
          propMFB3<-counter/length(MFB3)
        }
      }
    } 
   
    
    counter<-0
    #Data4 
    for(i in 1:sims){
      for(j in 1:length(years)){
        #print(MFB4[i,j])
        if(MFB4[i,j] < aveMFB){
          counter<-counter+1
          propMFB4<-counter/length(MFB4)
        }
      }
    }  
    
    
    counter<-0 
    #Data5 
    for(i in 1:sims){
      for(j in 1:length(years)){
        #print(MMB5[i,j])
        if(MFB5[i,j] < aveMFB){
          counter<-counter+1
          propMFB5<-counter/length(MFB5)
        }
      }
    }    
    
###############################################
    
    # Colors
    ##########################  
    HCR1<-rgb(.2, .2, 1,0.4)
    HCR2<-rgb(.2, .4, .6 ,0.4)
    HCR3<-rgb(.8, 0, .6 ,0.4)
    HCR4<-rgb(0, .2, .4,1)
    HCR5<-rgb(.8, 0, 0.4,1)
    ########################### 
    
    # TAC Cv vs MMBprob
    ####################### 
    par(mfrow = c(1, 1))
    par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
    plot(1,type="n",
         xlim=c(0,1),ylim=c(0,1),
         yaxs="i", xaxs="i") #no gap between 0 and the axis
    
    mtext(side=3, line=1, "TAC Trade Example 1" )
    mtext(side=1, line=2, "Proportion(MMB < aveMMB)")
    mtext(side=2, line=2, "TAC cv")
    
    abline(v =0.5, lty = 2, col = "grey70")
    
    #data1
    points(propMMB1,cvTAC1, pch = 19, col=HCR1, cex = 3 )
    points(propMMB1,cvTAC1, pch = 19, col=HCR1, cex = 3 )
    points(propMMB2,cvTAC2, pch = 19, col=HCR2, cex = 3 )
    points(propMMB2,cvTAC2, pch = 19, col=HCR2, cex = 3 )
    points(propMMB3,cvTAC3, pch = 19, col=HCR3, cex = 3 )
    points(propMMB4,cvTAC4, pch = 13, col=HCR4, cex = 3 )
    points(propMMB5,cvTAC5, pch = 13, col=HCR5, cex = 3 )
    
    legend("topright",c(lab1,lab2, lab3 ,lab4,lab5),pch=c(19,19,19,13,13),col=c(HCR1, HCR2,HCR3,HCR4,HCR5),cex=1)
    
  
    ###############  
    
    # TAC Cv vs MMBclose
    ####################### 
    par(mfrow = c(1, 1))
    par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
    plot(1,type="n",
         xlim=c(-.05,1),ylim=c(0,.5),
         yaxs="i", xaxs="i") #no gap between 0 and the axis
    
    mtext(side=3, line=1, "TAC Trade Example 3" )
    mtext(side=1, line=2, "Proportion(MMB < 0.25aveMMB)")
    mtext(side=2, line=2, "TAC cv")
    
    abline(v =0.5, lty = 2, col = "grey70")
    abline(v =0.0, lty = 2, col = "grey40")
    
    #data1
    points(propclose1,cvTAC1, pch = 19, col=HCR1, cex = 2 )
    points(propclose2,cvTAC2, pch = 19, col=HCR2, cex = 2 )
    points(propclose3,cvTAC3, pch = 19, col=HCR3, cex = 2 )
    points(propclose4,cvTAC4, pch = 13, col=HCR4, cex = 3 )
    points(propclose5,cvTAC5, pch = 13, col=HCR5, cex = 3 )
    
    legend("topright",c(lab1,lab2, lab3 ,lab4,lab5),pch=c(19,19,19,13,13),col=c(HCR1, HCR2,HCR3,HCR4,HCR5),cex=1)
    
    ############################ 
    
    #   Tac Ranges vs MMBprob
    ##############################
    
    plot(1,type="n",
         xlim=c(-.05,1),ylim=c(0,60),
         yaxs="i", xaxs="i") #no gap between 0 and the axis
    
    mtext(side=3, line=1, "TAC Trade Example 4" )
    mtext(side=1, line=2, "Proportion(MMB < aveMMB)")
    mtext(side=2, line=2, "TAC Median")
    
    abline(v =0.5, lty = 2, col = "grey70")
    abline(v =0, lty = 2, col = "grey70")
    #data1
    points(propclose1,TACMed1, pch = 19, col=HCR1, cex = 2 )
    arrows(propclose1, TACQuant1[1], propclose1, TACQuant1[2], length=0.05, angle=90, code=3, col = HCR1)
    
    points(propclose2,TACMed2, pch = 19, col=HCR2, cex = 2 )
    arrows(propclose2, TACQuant2[1], propclose2, TACQuant2[2], length=0.05, angle=90, code=3, col = HCR2)
    
    points(propclose3,TACMed3, pch = 19, col=HCR3, cex = 2 )
    arrows(propclose3, TACQuant3[1], propclose3, TACQuant3[2], length=0.05, angle=90, code=3, col = HCR3)
    
    points(propclose4,TACMed4, pch = 13, col=HCR4, cex = 2 )
    arrows(propclose4, TACQuant4[1], propclose4, TACQuant4[2], length=0.05, angle=90, code=3, col = HCR4)
    
    points(propclose5,TACMed5, pch = 13, col=HCR5, cex = 2 )
    arrows(propclose5, TACQuant5[1], propclose5, TACQuant5[2], length=0.05, angle=90, code=3, col = HCR5)
    
    legend("topright",c(lab1,lab2, lab3 ,lab4,lab5),pch=c(19,19,19,13,13),col=c(HCR1, HCR2,HCR3,HCR4,HCR5),cex=1)
    
    #######################################  
    
    
    
    
    ###############  
    
    # TAC Cv vs MMBclose
    ####################### 
    par(mfrow = c(1, 1))
    par(mar = c(1, 1, 2, 1), oma = c(2.5, 2.5, 1, 0.5))
    plot(1,type="n",
         xlim=c(0,1),ylim=c(0,1),
         yaxs="i", xaxs="i") #no gap between 0 and the axis
    
    mtext(side=3, line=1, "TAC Trade Example" )
    mtext(side=1, line=2, "Proportion(MFB < aveMFB)")
    mtext(side=2, line=2, "TAC cv")
    
    abline(v =0.5, lty = 2, col = "grey70")
    abline(v =0.0, lty = 2, col = "grey40")
    
    #data1
    points(propMFB1,cvTAC1, pch = 19, col=HCR1, cex = 2 )
    points(propMFB2,cvTAC2, pch = 19, col=HCR2, cex = 2 )
    points(propMFB3,cvTAC3, pch = 19, col=HCR3, cex = 2 )
    points(propMFB4,cvTAC4, pch = 13, col=HCR4, cex = 2 )
    points(propMFB5,cvTAC5, pch = 13, col=HCR5, cex = 2 )
    
    legend("topright",c(lab1,lab2, lab3 ,lab4,lab5),pch=c(19,19,19,13,13),col=c(HCR1, HCR2,HCR3,HCR4,HCR5),cex=1)
    
    ############################ 
    
    #   Tac Ranges vs MMBprob
    ##############################
    
    plot(1,type="n",
         xlim=c(0,1),ylim=c(0,50),
         yaxs="i", xaxs="i") #no gap between 0 and the axis
    
    mtext(side=3, line=1, "TAC Trade Example" )
    mtext(side=1, line=2, "Proportion(MMB < aveMMB)")
    mtext(side=2, line=2, "TAC Mean")
    
    abline(v =0.5, lty = 2, col = "grey70")
    
    #data1
    points(propMFB1,TACMed1, pch = 19, col=HCR1, cex = 2 )
    arrows(propMFB1, TACQuant1[1], propMFB1, TACQuant1[2], length=0.05, angle=90, code=3, col = HCR1)
    
    points(propMFB2,TACMed2, pch = 19, col=HCR2, cex = 2 )
    arrows(propMFB2, TACQuant2[1], propMFB2, TACQuant2[2], length=0.05, angle=90, code=3, col = HCR2)
    
    points(propMFB3,TACMed3, pch = 19, col=HCR3, cex = 2 )
    arrows(propMFB3, TACQuant3[1], propMFB3, TACQuant3[2], length=0.05, angle=90, code=3, col = HCR3)
    
    points(propMFB4,TACMed4, pch = 13, col=HCR4, cex = 2 )
    arrows(propMFB4, TACQuant4[1], propMFB4, TACQuant4[2], length=0.05, angle=90, code=3, col = HCR4)
    
    points(propMFB5,TACMed5, pch = 13, col=HCR5, cex = 2 )
    arrows(propMFB5, TACQuant5[1], propMFB5, TACQuant5[2], length=0.05, angle=90, code=3, col = HCR5)
    
    legend("topright",c(lab1,lab2, lab3 ,lab4,lab5),pch=c(19,19,19,13,13),col=c(HCR1, HCR2,HCR3,HCR4,HCR5),cex=1)
    
    #######################################  
   }  
  






Trades(data1=HCR2_1_OpMod, data2=HCR2_2_OpMod, data3=HCR2_3_OpMod, data4=HCR3_OpMod, data5=HCR2_1_OpMod, sims=100, modTyp=2, 2018, 2118, 
  lab1 = "Male 10", 
  lab2 = "Male 15", 
  lab3 = "Male 20", 
  lab4 = "TAC = ABC",
  lab5 = "test",
  PlotType = 2)
