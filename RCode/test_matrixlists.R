

###TEST 

setwd("C:\\MSE_Runs\\testMSE_19\\run001\\2019.EstMod")

test1<-source("perfMetricsEstMod_2019.txt")
typeof(test1)

print(test1$value$DISCARDS)

Discards<-test1$value$DISCARDS
typeof(Discards)

dim(as.array(Discards)) #array fishery, years, sex, maturity, shell, size
Discards[4,,1 , , ,]

Discards[,,1,2,2,]

setwd("C:\\MSE_Runs\\testMSE_19\\run001\\2020.EstMod")

test2<-source("perfMetricsEstMod_2020.txt")

library(data.table)


Discards2<-test2$value$DISCARDS

Disc<-matrix(list(), nrow= 1, ncol=2)

Disc[[1,1]]<-Discards
Disc[[1,2]]<-Discards2
# No issue with vectors or matrices with the single value perf metrics 
TAC<-vector()

TAC[1]<-test1$value$TAC
TAC[2]<-test2$value$TAC


CATCHES<-list()
CATCHES[[1]]<-test1$value$CATCH
CATCHES[[2]]<-test2$value$CATCH

