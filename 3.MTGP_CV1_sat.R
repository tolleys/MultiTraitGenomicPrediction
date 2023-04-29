# FILENAME: TCZm_PHP02_Even.Het_A.BGLR.Genomic.Prediction_Server.R
setwd("~/Documents/Purdue/PhD/Research/Chp2_MTGP/DATA/")
########Loading libraries###########
#install.packages('caret',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
install.packages(c("caret", "AGHmatrix", "readxl", "bWGR", "dplyr"))
library('caret', verbose=FALSE)
library("AGHmatrix", verbose=FALSE)
library("readxl", verbose=FALSE)
library("bWGR", verbose=FALSE)
library("dplyr", verbose=FALSE)
library("remotes")
remotes::install_github("gdlc/BGLR-R")
library(BGLR)





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################[2FACC_Temp]###########################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
Grin1 <- read.csv("Grin1.csv", header=T)

DataOutput <- data.frame()
DataOutput1 <- data.frame()

for(n in 1:10){
    set.seed(24+n)
  
  BLUEs <- as.data.frame(read_excel("1.TCZm17_18_EP.BLUEs.xlsx", sheet=1))
  

  BLUEs_2FACC <- subset(BLUEs, BLUEs$Pedigree2=="2FACC")
  BLUEs_2FACC <- merge(Grin1, BLUEs_2FACC, by="Grin1")
  BLUEs_2FACC$REFREFYLD18 <- as.numeric(BLUEs_2FACC$REFYLD18)
      
    BLUEs_2FACC_IO   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="Iodent")
    BLUEs_2FACC_IO   <- subset(BLUEs_2FACC_IO, BLUEs_2FACC_IO$REFYLD18!="NA")
      flds <- createFolds(BLUEs_2FACC_IO$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
      BLUEs_2FACC_IO <- cbind(BLUEs_2FACC_IO, flds)
      
    BLUEs_2FACC_NS   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="NS")
    BLUEs_2FACC_NS   <- subset(BLUEs_2FACC_NS, BLUEs_2FACC_NS$REFYLD18!="NA")
      flds <- createFolds(BLUEs_2FACC_NS$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
      BLUEs_2FACC_NS <- cbind(BLUEs_2FACC_NS, flds)
      
    BLUEs_2FACC_SS   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="SS")
    BLUEs_2FACC_SS   <- subset(BLUEs_2FACC_SS, BLUEs_2FACC_SS$REFYLD18!="NA")
      flds <- createFolds(BLUEs_2FACC_SS$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
      BLUEs_2FACC_SS <- cbind(BLUEs_2FACC_SS, flds)

      BLUEs_2FACC <- rbind(BLUEs_2FACC_IO, BLUEs_2FACC_NS, BLUEs_2FACC_SS)
      BLUEs_2FACC <- BLUEs_2FACC %>%
        arrange(BLUEs_2FACC$Grin1)
  
      Grin1 <- as.data.frame(BLUEs_2FACC$Grin1)
      Grin1 <- Grin1 %>%
        rename(Grin1 = `BLUEs_2FACC$Grin1`)
      

      FACC2_Add <- read.table('2FACC_Add.txt', header=T, sep=" ")
      Grin1_2FACC_Add <- FACC2_Add$Hybrids
      Grin1_2FACC_Add <- gsub("\t", "", x=Grin1_2FACC_Add, perl = T)
      Grin1_2FACC_Add <- gsub("0_", "", x=Grin1_2FACC_Add, perl=T)

      FACC2_Add <- FACC2_Add[,-c(1)]
      rownames(FACC2_Add) <- Grin1_2FACC_Add
      FACC2_Add <- cbind(Grin1_2FACC_Add, FACC2_Add)
      
      FACC2_Add <- FACC2_Add %>%
        rename(Grin1 = Grin1_2FACC_Add)
      
      FACC2_Add <- merge(Grin1, FACC2_Add, by="Grin1") #Match geno and pheno in this run

      rownames(FACC2_Add) <- FACC2_Add$Grin1
      FACC2_Add <- FACC2_Add[,-c(1)]
FACC2_Add <- as.matrix(FACC2_Add)

FACC2_Add_VanRaden <- Gmatrix(SNPmatrix = FACC2_Add, missingValue=-9,
                              maf=0.05, method="VanRaden")

EVD_2FACC_Add_VanRaden <- eigen(FACC2_Add_VanRaden)
rownames(EVD_2FACC_Add_VanRaden$vectors) <- rownames(FACC2_Add_VanRaden)


      for(i in 1:5){

        pheno <- BLUEs_2FACC #make a copy that we can use for our analysis
        pheno1 <- pheno %>%
          select(REFYLD18, PHTYLD, PHTKPE, EARAREA) #<<------Change
        pheno1$REFYLD18 <- as.numeric(pheno1$REFYLD18)
        pheno1$REFYLD18 <- scale(pheno1$REFYLD18, center=TRUE, scale=TRUE)
        pheno1$PHTYLD <- scale(pheno1$PHTYLD, center=TRUE, scale=TRUE)
        #pheno1$EARLGT <- scale(pheno1$EARLGT, center=TRUE, scale=TRUE)
        pheno1$PHTKPE <- scale(pheno1$PHTKPE, center=TRUE, scale=TRUE)
        pheno1$EARAREA <- scale(pheno1$EARAREA, center=TRUE, scale=TRUE)
        pheno1$flds <- pheno$flds
  
        pheno1$REFYLD18[pheno1$flds == i] <- NA
        pheno1$PHTYLD[pheno1$flds == i] <- NA
        pheno1$PHTKPE[pheno1$flds == i] <- NA
        #pheno1$EARLGT[pheno1$flds == i] <- NA
        pheno1$EARAREA[pheno1$flds == i] <- NA
        pheno1 <- as.matrix(pheno1[,-5])
        rownames(pheno1) <- BLUEs_2FACC$Grin1
        
        fm=Multitrait(y=pheno1,
                ETA=list(list(EVD=EVD_2FACC_Add_VanRaden,
                                model='RKHS')),
                nIter=100000,
                burnIn=25000,
                thin=100,
                verbose=F,
                saveAt='output')


        yhat1 <- fm$ETAHat[,1]
        data1 <- cbind(pheno,yhat1)
        Correlation <- cor(data1$yhat1[data1$flds==i], data1$REFYLD18[data1$flds==i]) ### <<<------ Need to make the change here
        trait <- "REFYLD18"
        Repeated <- n


        data1 <- cbind(Repeated,i,trait,Correlation)
        DataOutput <- rbind(DataOutput, data1)
        

        
        print(n)
        print(i)
      }
      
    }
    



mean(as.numeric(DataOutput$Correlation), na.rm=TRUE)
sd(as.numeric(DataOutput$Correlation))


save.image(file='2FACC_MTM_REFYLD18_NAPHTYLD.PHTKPE.EARAREA_sat.RData')





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##################################[PHP02_Temp]#####################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
DataOutput <- data.frame()
DataOutput1 <- data.frame()

for(n in 1:10){
  set.seed(24+n)
  
  BLUEs <- as.data.frame(read_excel("1.TCZm17_18_EP.BLUEs.xlsx", sheet=1))
  BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")

  
  BLUEs_PHP02_IO   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="Iodent")
  BLUEs_PHP02_IO <- subset(BLUEs_PHP02_IO, BLUEs_PHP02_IO$REFYLD18!="NA")
  flds <- createFolds(BLUEs_PHP02_IO$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
    BLUEs_PHP02_IO <- cbind(BLUEs_PHP02_IO, flds)
  
  BLUEs_PHP02_NS   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="NS")
  BLUEs_PHP02_NS <- subset(BLUEs_PHP02_NS, BLUEs_PHP02_NS$REFYLD18!="NA")
  flds <- createFolds(BLUEs_PHP02_NS$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
    BLUEs_PHP02_NS <- cbind(BLUEs_PHP02_NS, flds)
  
  BLUEs_PHP02_SS   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="SS")
  BLUEs_PHP02_SS <- subset(BLUEs_PHP02_SS, BLUEs_PHP02_SS$REFYLD18!="NA")
  flds <- createFolds(BLUEs_PHP02_SS$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
   BLUEs_PHP02_SS <- cbind(BLUEs_PHP02_SS, flds)
   
  BLUEs_PHP02 <- rbind(BLUEs_PHP02_IO, BLUEs_PHP02_NS, BLUEs_PHP02_SS)
  
  BLUEs_PHP02 <- BLUEs_PHP02 %>%
    arrange(BLUEs_PHP02$Grin1)
  

  Grin1 <- as.data.frame(BLUEs_PHP02$Grin1)
  Grin1 <- Grin1 %>%
    rename(Grin1 = `BLUEs_PHP02$Grin1`)
  

  PHP02_Add <- read.table('PHP02_Add.txt', header=T, sep=" ")
  Grin1_PHP02_Add <- PHP02_Add$Hybrids
  Grin1_PHP02_Add <- gsub("\t", "", x=Grin1_PHP02_Add, perl = T)
  Grin1_PHP02_Add <- gsub("0_", "", x=Grin1_PHP02_Add, perl=T)
  
  PHP02_Add <- PHP02_Add[,-c(1)]
  rownames(PHP02_Add) <- Grin1_PHP02_Add
  PHP02_Add <- cbind(Grin1_PHP02_Add, PHP02_Add)
  
  PHP02_Add <- PHP02_Add %>%
    rename(Grin1 = Grin1_PHP02_Add)
  
  PHP02_Add <- merge(Grin1, PHP02_Add, by="Grin1") #Match geno and pheno in this run
  rownames(PHP02_Add) <- PHP02_Add$Grin1
  PHP02_Add <- PHP02_Add[,-c(1)]
  PHP02_Add <- as.matrix(PHP02_Add)
  
  
  PHP02_Add_VanRaden <- Gmatrix(SNPmatrix = PHP02_Add, missingValue=-9,
                                maf=0.05, method="VanRaden")
  EVD_PHP02_Add_VanRaden <- eigen(PHP02_Add_VanRaden)
  rownames(EVD_PHP02_Add_VanRaden$vectors) <- rownames(PHP02_Add_VanRaden)

  
  for(i in 1:5){
    
    pheno <- BLUEs_PHP02 #make a copy that we can use for our analysis
    pheno$REFYLD18 <- as.numeric(pheno$REFYLD18)
    pheno1 <- pheno %>%
      select(REFYLD18, PHTYLD, PHTKPE, EARAREA) #<<------Change
    pheno1$REFYLD18 <- scale(pheno1$REFYLD18, center=TRUE, scale=TRUE)
    pheno1$PHTYLD <- scale(pheno1$PHTYLD, center=TRUE, scale=TRUE)
    #pheno1$EARLGT <- scale(pheno1$EARLGT, center=TRUE, scale=TRUE)
    pheno1$PHTKPE <- scale(pheno1$PHTKPE, center=TRUE, scale=TRUE)
    pheno1$EARAREA <- scale(pheno1$EARAREA, center=TRUE, scale=TRUE)
    pheno1$flds <- pheno$flds
    
    pheno1$REFYLD18[pheno1$flds == i] <- NA
    #pheno1$PHTYLD[pheno1$flds == i] <- NA
    #pheno1$PHTKPE[pheno1$flds == i] <- NA
    #pheno1$EARLGT[pheno1$flds == i] <- NA
    #pheno1$EARAREA[pheno1$flds == i] <- NA
    pheno1 <- as.matrix(pheno1[,-5])
    rownames(pheno1) <- BLUEs_PHP02$Grin1
    
    fm=Multitrait(y=pheno1,
                  ETA=list(list(EVD=EVD_PHP02_Add_VanRaden,
                                model='RKHS')),
                  nIter=100000,
                  burnIn=25000,
                  thin=100,
                  verbose=F,
                  saveAt='output')
    
    
    yhat1 <- fm$ETAHat[,1]
    data1 <- cbind(pheno,yhat1)
    Correlation <- cor(data1$yhat1[data1$flds==i], data1$REFYLD18[data1$flds==i]) ### <<<------ Need to make the change here
    trait <- "REFYLD18"
    Repeated <- n
    
    
    data1 <- cbind(Repeated,i,trait,Correlation)
    DataOutput <- rbind(DataOutput, data1)
    
    
    
    print(n)
    print(i)
  }
  
}

mean(as.numeric(DataOutput$Correlation), na.rm=TRUE)
sd(as.numeric(DataOutput$Correlation))


save.image(file='PHP02_Temp_MTM_REFYLD18_PHTYLD.PHTKPE.EARAREA_sat.RData')



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##################################[PHP02_Trop]#####################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
DataOutput <- data.frame()
DataOutput1 <- data.frame()

for(n in 1:10){
  set.seed(24+n)
  
  BLUEs <- as.data.frame(read_excel("1.TCZm17_18_EP.BLUEs.xlsx", sheet=1))
  BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")
  
      BLUEs_PHP02_DTMA <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="DTMA")
      BLUEs_PHP02_DTMA <- subset(BLUEs_PHP02_DTMA, BLUEs_PHP02_DTMA$REFYLD18!="NA")
        flds <- createFolds(BLUEs_PHP02_DTMA$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
        BLUEs_PHP02_DTMA <- cbind(BLUEs_PHP02_DTMA, flds)


  BLUEs_PHP02 <- rbind(BLUEs_PHP02_DTMA)
  
  BLUEs_PHP02 <- BLUEs_PHP02 %>%
    arrange(BLUEs_PHP02$Grin1)
  

  Grin1 <- as.data.frame(BLUEs_PHP02$Grin1)
  Grin1 <- Grin1 %>%
    rename(Grin1 = `BLUEs_PHP02$Grin1`)
  

  PHP02_Add <- read.table('PHP02_Add.txt', header=T, sep=" ")
  Grin1_PHP02_Add <- PHP02_Add$Hybrids
  Grin1_PHP02_Add <- gsub("\t", "", x=Grin1_PHP02_Add, perl = T)
  Grin1_PHP02_Add <- gsub("0_", "", x=Grin1_PHP02_Add, perl=T)
  
  PHP02_Add <- PHP02_Add[,-c(1)]
  rownames(PHP02_Add) <- Grin1_PHP02_Add
  PHP02_Add <- cbind(Grin1_PHP02_Add, PHP02_Add)
  
  PHP02_Add <- PHP02_Add %>%
    rename(Grin1 = Grin1_PHP02_Add)
  
  PHP02_Add <- merge(Grin1, PHP02_Add, by="Grin1") #Match geno and pheno in this run
  rownames(PHP02_Add) <- PHP02_Add$Grin1
  PHP02_Add <- PHP02_Add[,-c(1)]
  PHP02_Add <- as.matrix(PHP02_Add)
  
  
  PHP02_Add_VanRaden <- Gmatrix(SNPmatrix = PHP02_Add, missingValue=-9,
                                maf=0.05, method="VanRaden")
  EVD_PHP02_Add_VanRaden <- eigen(PHP02_Add_VanRaden)
  rownames(EVD_PHP02_Add_VanRaden$vectors) <- rownames(PHP02_Add_VanRaden)
  
  
  for(i in 1:5){
    
    pheno <- BLUEs_PHP02 #make a copy that we can use for our analysis
    pheno$REFYLD18 <- as.numeric(pheno$REFYLD18)
    pheno1 <- pheno %>%
      select(PHTYLD, EARLGT) #<<------Change
    #pheno1$REFYLD18 <- as.numeric(pheno1$REFYLD18)
    #pheno1$REFYLD18 <- scale(pheno1$REFYLD18, center=TRUE, scale=TRUE)
    pheno1$PHTYLD <- scale(pheno1$PHTYLD, center=TRUE, scale=TRUE)
    pheno1$EARLGT <- scale(pheno1$EARLGT, center=TRUE, scale=TRUE)
    #pheno1$PHTKPE <- scale(pheno1$PHTKPE, center=TRUE, scale=TRUE)
    #pheno1$EARAREA <- scale(pheno1$EARAREA, center=TRUE, scale=TRUE)
    pheno1$flds <- pheno$flds
    
    #pheno1$REFYLD18[pheno1$flds == i] <- NA
    pheno1$PHTYLD[pheno1$flds == i] <- NA
    #pheno1$PHTKPE[pheno1$flds == i] <- NA
    #pheno1$EARLGT[pheno1$flds == i] <- NA
    #pheno1$EARAREA[pheno1$flds == i] <- NA
    pheno1 <- as.matrix(pheno1[,-3])
    rownames(pheno1) <- BLUEs_PHP02$Grin1
    
    fm=Multitrait(y=pheno1,
                  ETA=list(list(EVD=EVD_PHP02_Add_VanRaden,
                                model='RKHS')),
                  nIter=100000,
                  burnIn=25000,
                  thin=100,
                  verbose=F,
                  saveAt='output')
    
    
    yhat1 <- fm$ETAHat[,1]
    data1 <- cbind(pheno,yhat1)
    Correlation <- cor(data1$yhat1[data1$flds==i], data1$PHTYLD[data1$flds==i]) ### <<<------ Need to make the change here
    trait <- "PHTYLD"
    Repeated <- n
    
    
    data1 <- cbind(Repeated,i,trait,Correlation)
    DataOutput <- rbind(DataOutput, data1)
    
    
    
    print(n)
    print(i)
  }
  
}

mean(as.numeric(DataOutput$Correlation), na.rm=TRUE)
sd(as.numeric(DataOutput$Correlation))

save.image(file='PHP02_DTMA_MTM_PHTYLD_EARLGT_sat.RData')



