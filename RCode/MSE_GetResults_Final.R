# R script for Reading data Files and putting them together 

#This set of Functions Pulls the data from the runs out for use by the state.
#Ultimately it will be run 13 times, once for each HCR.

########
#For Andre!!!!!
# SET THE WORKING DIRECTORY TO WHERE THIS FILE IS SAVED AND PUT
# THE UNZIPPED ATTACHED DATA FOLDER In THE SAME FOLDER

#setwd()
#setwd("C:/MSE_Runs/")
setwd("/bigVol/model/MSE_Runs")

#####################################################    
## Create the name of the folders for simulations 
#################################################### 

getSimsFolders <- function(SimNum) {
  SimNum <- SimNum
  run <- vector(length = SimNum)
  for (i in 1:length(run)) {
    if (i < 10) {
      run[i] <- paste("run00", i, sep = "")
    }
    if (i >= 10 && i < 100) {
      run[i] <- paste("run0", i, sep = "")
    }
    
    if (i >= 100) {
      run[i] <- paste("run", i, sep = "")
      
    }
  }
  return(run)
}

getSimsFolders(10)


#####################################################
## Create the name of the OpMod and EstMod folders 
#####################################################

getFolderNames<-function(Type ="Type", StartYear, EndYear){
  Type <- Type
  #StartYear<-StartYear
  #EndYear<-EndYear
  #Years<-seq(StartYear,EndYear)  
  mod<-vector()
  
  if(Type == "OpMod"){
    StartYear<-StartYear
    EndYear<-EndYear-1
    Years<-seq(StartYear,EndYear)
  }else{
    StartYear<-StartYear+1
    EndYear<-EndYear
    Years<-seq(StartYear,EndYear)
  }
  
  
  for(i in 1:length(Years)) {
    
    if (Type == "OpMod") {
      # StartYear<-StartYear
      #EndYear<-EndYear-1
      #Years<-seq(StartYear,EndYear)
      year=(StartYear-1)+i
      mod[i] <- paste(year, ".OpMod", sep = "")
    }
    if (Type == "EstMod") {
      #StartYear<-StartYear+1
      #EndYear<-EndYear
      year=(StartYear-1)+i
      mod[i] <- paste(year, ".EstMod", sep = "")
    }
  }
  #EndFunctin
  #print(mod)
  return(mod)
}

getFolderNames(Type="EstMod", 2018, 2028)

#######################################################
# Create the names of the performance metric files 
######################################################

getFileName<-function(Type = "Type", StartYear, EndYear){
  Type <- Type
  
  if( Type == "OpMod"){
    StartYear<-StartYear
    EndYear<-EndYear-1
    Years<-seq(StartYear,EndYear) 
  }else{
    StartYear<-StartYear+1
    EndYear<-EndYear
    Years<-seq(StartYear,EndYear)
  }
  
  fileName<-vector()
  
  for(i in 1:length(Years)) {
    if (Type == "OpMod") {
      year=(StartYear-1)+i
      fileName[i] <-"perfMetricsOPMod.txt"
    }
    if (Type == "EstMod") {
      year=(StartYear-1)+i
      fileName[i] <- paste("perfMetricsEstMod_", year, ".txt", sep = "")
    }
  }
  
  return(fileName)
  
}

getFileName(Type="EstMod", 2018, 2028)
getFileName(Type="OpMod", 2018,2038)

#########################################
#Get Starting Directory
########################################

startingDir<-function(StartingDir, folderName){
  folderName<-folderName
  StartingDir<-StartingDir
  NewDir<-paste(StartingDir, "/", folderName, sep ="")
  setwd(NewDir)
  getwd()
  #print("New Working Directory Set")
}


##################################
# Clean NaNs 
#################################

cleanData<-function(file){
  
  data<-read.csv2(file, header = TRUE)
  
  data<-readLines(file, -1)
  if(data[4] == "TAC=-nan, "){data[4]="TAC=-9999, "}
  if(data[8] == "Fmsy=-inf, "){data[8]="Fmsy= Inf, "}
  if(data[9] == "Fofl=-inf, "){data[9]="Fofl= Inf, "}
  
  writeLines(data,file)
  
}


#currentDir<-(getwd())
#startingDir(currentDir, folderName = "testMSE_19")

#########################################
# EXTRACT THE DATA
########################################

getData<-function(Type = Type, StartYear, EndYear, SimNum, folderName, FileName){
  Type <- Type
  
  FileName<-FileName ## What kind of data file to put out? Csv? Txt? Dat?
  
  StartYear<-StartYear
  EndYear<-EndYear 
  
  print(StartYear)
  print(EndYear)
  
  SimNum <- SimNum
  folderName<-folderName
  
  currentDir<-(getwd())
  startingDir(currentDir, folderName = folderName)
  
  print(run<-getSimsFolders(SimNum)) # Number of simulations
  print(mod<-getFolderNames(Type, StartYear, EndYear)) # number of years 
  print(file<-getFileName(Type, StartYear, EndYear))
  
  # Create Matrices to hold data 
  
  TAC_SET<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(TAC_SET)<-run
  colnames(TAC_SET)<-mod
  
  TAC<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(TAC)<-run
  colnames(TAC)<-mod
  
  B0<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(B0)<-run
  colnames(B0)<-mod
  
  Bmsy<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(Bmsy)<-run
  colnames(Bmsy)<-mod
  
  OFL<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(OFL)<-run
  colnames(OFL)<-mod
  
  Fmsy<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(Fmsy)<-run
  colnames(Fmsy)<-mod
  
  Fofl<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(Fofl)<-run
  colnames(Fofl)<-mod
  
  MMB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(MMB)<-run
  colnames(MMB)<-mod
  
  MFB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(MFB)<-run
  colnames(MFB)<-mod
  
  ELMB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(ELMB)<-run
  colnames(ELMB)<-mod
  
  ELMB_State<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(ELMB_State)<-run
  colnames(ELMB_State)<-mod
  
  AveRec<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(AveRec)<-run
  colnames(AveRec)<-mod
  
  Rec<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(Rec)<-run
  colnames(Rec)<-mod
  
  MFCB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(MFCB)<-run
  colnames(MFCB)<-mod
  
  IFCB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(IFCB)<-run
  colnames(IFCB)<-mod
  
  ELMC<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(ELMC)<-run
  colnames(ELMC)<-mod
  
  MMCB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(MMCB)<-run
  colnames(MMCB)<-mod
  
  IMCB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(IMCB)<-run
  colnames(IMCB)<-mod
  
  MFDB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(MFDB)<-run
  colnames(MFDB)<-mod
  
  IFDB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(IFDB)<-run
  colnames(IFDB)<-mod
  
  ELMD<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(ELMD)<-run
  colnames(ELMD)<-mod
  
  MMDB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(MMDB)<-run
  colnames(MMDB)<-mod
  
  IMDB<-matrix(data= NA,nrow =length(run), ncol = length(mod))
  row.names(IMDB)<-run
  colnames(IMDB)<-mod
  
  
  for(i in 1:length(run)){
    print("FirstPart")
    dir<-getwd()
    startingDir(dir, folderName = run[i])
    for(j in 1:length(mod)){
      print("SecondPart")
      dir<-getwd()
      print(startingDir(dir,folderName = mod[j]))
      print(file[j])
      cleanData(file[j])
      data<-source(file[j])
      print(data)
      TAC_SET[i,j]<-data$value$TACset
      TAC[i,j]<-data$value$TAC
      B0[i,j]<-data$value$B0
      Bmsy[i,j]<-data$value$Bmsy
      OFL[i,j]<-data$value$OFL
      Fmsy[i,j]<-data$value$Fmsy
      Fofl[i,j]<-data$value$Fofl
      MMB[i,j]<-data$value$MMB
      MFB[i,j]<-data$value$MFB
      ELMB[i,j]<-data$value$ELMB
      ELMB_State[i,j]<-data$value$ELMB_State
      Rec[i,j]<-data$value$Rec
      if(Type=="OpMod"){
        AveRec[i,j]<-data$value$RecAve
      }else{
        AveRec[i,j]<-data$value$AveRec
      }
      MFCB[i,j]<-data$value$MFCB
      IFCB[i,j]<-data$value$IFCB
      ELMC[i,j]<-data$value$ELMC
      MMCB[i,j]<-data$value$MMCB
      IMCB[i,j]<-data$value$IMCB
      MFDB[i,j]<-data$value$MFDB
      IFDB[i,j]<-data$value$IFDB
      ELMD[i,j]<-data$value$ELMD
      MMDB[i,j]<-data$value$MMDB
      IMDB[i,j]<-data$value$IMDB
      setwd('..')
    }
    setwd('..')
    print(getwd())
  }
  
  #Return this   
  print("Assembling Data")
  Data <-  list("TAC Set"=TAC_SET,
                "TAC" =TAC, 
                "OFL"=OFL, 
                "B0"=B0,
                "Bmsy"= Bmsy,
                "Fmsy"=Fmsy,
                "Fofl"=Fofl,
                "MMB"=MMB, 
                "MFB"=MFB, 
                "ELMB"=ELMB,
                "ELMB_State"=ELMB_State,
                "AveRec"=AveRec, 
                "Rec"=Rec,
                "MFCB"=MFCB,
                "IFCB"=IFCB,
                "ELMC"=ELMC,
                "MMCB"=MMCB,
                "IMCB"=IMCB,
                "MFDB"=MFDB,
                "IFDB"=IFDB,
                "ELMD"=ELMD,
                "MMDB"=MMDB,
                "IMDB"=IMDB)   
  
  saveRDS(Data, FileName, compress = FALSE)
  setwd("..")
  return(Data)  
}


setwd("/bigVol/model/MSE_Runs")
#setwd("C:/MSE_Runs/")


getData(Type = "OpMod", 2018, 2118, 50, "StatQuo_2", FileName = "HCR7_OpMod2")
getData(Type = "EstMod", 2018, 2118, 50, "StatQuo_2", FileName = "HCR7_EstMod2")

# getData(Type = "OpMod", 2018, 2028, 3, "Dim_10225_30", FileName = "Dim43_op")
# getData(Type = "EstMod", 2018, 2028, 3, "Dim_10225_30", FileName = "Dim43_est")
# 
# getData(Type = "OpMod", 2018, 2028, 2, "MaleR4", FileName = "MaleR4_op")
# getData(Type = "EstMod", 2018, 2028, 2, "MaleR4", FileName = "MaleR4_est")




#look<-getData(Type = "EstMod", 2018, 2118, 100, "StatusQuo", FileName = "HCR7_EstMod")
#getData(Type = "OpMod", 2018, 2118, 100, "StatusQuo", FileName = "HCR7_OpMod")

#Save the information for later use 
#saveRDS(test, "testfile2.RDS", compress = FALSE)
# setwd("C:/MSE_Runs/Cmsy_test")
# test1 <- readRDS("cmsytest")
# test2 <- readRDS("cmsytest2")


