# FILENAME: TCZm_PHP02_Even.Het_A.BGLR.Genomic.Prediction_Server.R
setwd("~/Documents/Purdue/PhD/Research/Chp2_MTGP/DATA/")
########Loading libraries###########
#install.packages('caret',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
#install.packages(c("caret", "BGLR", AGHmatrix", "readxl", "dplyr"))
library('caret', verbose=FALSE)
library("BGLR", verbose=FALSE)
library("AGHmatrix", verbose=FALSE)
library("readxl", verbose=FALSE)
library("dplyr", verbose=FALSE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####################################[2FACC_Temp]##########################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
DataOutput <- data.frame()
DataOutput1 <- data.frame()


for(n in 1:10){
    set.seed(24+n)
  
  BLUEs <- as.data.frame(read_excel("1.TCZm17_18_EP.BLUEs.xlsx", sheet=1))
  

  
  BLUEs_2FACC <- subset(BLUEs, BLUEs$Pedigree2=="2FACC")
  colnum=c(9:ncol(BLUEs_2FACC))
  
  BLUEs_2FACC$REFYLD18 <- as.numeric(BLUEs_2FACC$REFYLD18)
      
    BLUEs_2FACC_IO   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="Iodent")
    BLUEs_2FACC_IO <- subset(BLUEs_2FACC_IO, BLUEs_2FACC_IO$REFYLD18!="NA")
      #BLUEs_2FACC_IO <- BLUEs_2FACC_IO[sample(nrow(BLUEs_2FACC_IO),14), ]
      flds <- createFolds(BLUEs_2FACC_IO$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
      BLUEs_2FACC_IO <- cbind(BLUEs_2FACC_IO, flds)
      
    BLUEs_2FACC_NS   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="NS")
    BLUEs_2FACC_NS <- subset(BLUEs_2FACC_NS, BLUEs_2FACC_NS$REFYLD18!="NA")
      #BLUEs_2FACC_NS <- BLUEs_2FACC_NS[sample(nrow(BLUEs_2FACC_NS),14), ]
      flds <- createFolds(BLUEs_2FACC_NS$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
      BLUEs_2FACC_NS <- cbind(BLUEs_2FACC_NS, flds)
      
    BLUEs_2FACC_SS   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="SS")
    BLUEs_2FACC_SS <- subset(BLUEs_2FACC_SS, BLUEs_2FACC_SS$REFYLD18!="NA")
      #BLUEs_2FACC_SS <- BLUEs_2FACC_SS[sample(nrow(BLUEs_2FACC_SS),14), ]
      flds <- createFolds(BLUEs_2FACC_SS$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
      BLUEs_2FACC_SS <- cbind(BLUEs_2FACC_SS, flds)
#Phenotype File Prepared--------------------------------------------------------------   
BLUEs_2FACC <- rbind(BLUEs_2FACC_IO, BLUEs_2FACC_NS, BLUEs_2FACC_SS)
      BLUEs_2FACC <- BLUEs_2FACC %>%
        arrange(BLUEs_2FACC$Grin1)
#-------------------------------------------------------------------------------------   
      
      Grin1 <- as.data.frame(BLUEs_2FACC$Grin1)
      Grin1 <- Grin1 %>%
        rename(Grin1 = `BLUEs_2FACC$Grin1`)
      
#Genotype File Preparation------------------------------------------------------------   
#Add-----
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

      Grin1 <- data.frame(rownames(FACC2_Add))
      colnames(Grin1) = "Grin1"
            FACC2_Add <- FACC2_Add[,-c(1)]
FACC2_Add <- as.matrix(FACC2_Add)

FACC2_Add_VanRaden <- Gmatrix(SNPmatrix = FACC2_Add, missingValue=-9,
                              maf=0.05, method="VanRaden")

EVD_2FACC_Add_VanRaden <- eigen(FACC2_Add_VanRaden)
rownames(EVD_2FACC_Add_VanRaden$vectors) <- rownames(FACC2_Add_VanRaden)

#-------------------------------------------------------------------------------------

    for(j in 1:26){
      x=colnum[j]  #set the current [i] column as "x"
      trait=colnames(BLUEs_2FACC)[x] #sets the current column header as the trait name
      pheno <- merge(Grin1, BLUEs_2FACC, by="Grin1") #make a copy that we can use for our analysis
      colnames(pheno)[x]="y"  #renames the trait variable as "y" for the model analysis below
      
      for(i in 1:5){

        pheno1 <- pheno
        pheno1$y <- scale(pheno$y, center = TRUE, scale=TRUE)
        pheno1$y[pheno1$flds == i] <- NA


        fm=BGLR(y=pheno1$y,
                ETA=list(G=list(V=EVD_2FACC_Add_VanRaden$vectors,
                                d=EVD_2FACC_Add_VanRaden$values,
                                model='RKHS')),
                nIter=100000,
                burnIn=25000,
                thin=100,
                verbose=T)

        yhat1 <- fm$yHat
        data1 <- cbind(pheno,yhat1)
        Correlation <- cor(data1$yhat1[data1$flds==i], data1$y[data1$flds==i])
    
        Repeated <- n

        varE=scan('varE.dat')
        varA=scan('ETA_G_varU.dat')
        h2=(varA)/(varA+varE)
        
        
        pdf(file=paste0("~/Documents/Purdue/PhD/Research/Chp2_MTGP/PLOTS/TracePlots/2FACC_Temp_", trait,"_",i,"_",n,".pdf"),
            width=8,
            height=6)
        par(mfrow=c(2,2))
        plot(varE,type='l', col="black", ylab="Error Variance")
        plot(varA, type='l', col="navy", ylab="Genotypic Variance")
        plot(h2, type="l", col="orange2", ylab="Heritability")
        dev.off()
        
        h2_mean=(mean(varA[251:1000]))/(mean(varA[251:1000])+mean(varE[251:1000]))
        
        data1 <- cbind(Repeated,i,trait,Correlation, h2_mean)
        DataOutput <- rbind(DataOutput, data1)
        
        data2 <- cbind(Repeated, i, trait, varE, varA, h2)
        DataOutput1 <- rbind(DataOutput1, data2)
        
        print(n)
        print(j)
        print(i)
      }
      
    }
    
}
DataOutput


write.csv(DataOutput, "2FACC_Scaled_A_Single.Trait.Prediction_sat.csv")
write.csv(DataOutput1, '2FACC_Scaled_A_Variance.Comp_sat.csv')


save.image(file='2FACC_Scaled_A_STGP_sat.RData')



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#####################################[PHP02_Temp]##########################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
DataOutput <- data.frame()
DataOutput1 <- data.frame()

for(n in 1:10){
  set.seed(24+n)
  
  BLUEs <- as.data.frame(read_excel("1.TCZm17_18_EP.BLUEs.xlsx", sheet=1))
  BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")
  colnum=c(9:ncol(BLUEs_PHP02))
  
  BLUEs_PHP02_IO   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="Iodent")
  BLUEs_PHP02_IO <- subset(BLUEs_PHP02_IO, BLUEs_PHP02_IO$REFYLD18!="-999")
  #BLUEs_PHP02_IO <- BLUEs_PHP02_IO[sample(nrow(BLUEs_PHP02_IO),42), ]
  flds <- createFolds(BLUEs_PHP02_IO$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
  BLUEs_PHP02_IO <- cbind(BLUEs_PHP02_IO, flds)
  
  BLUEs_PHP02_NS   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="NS")
  BLUEs_PHP02_NS <- subset(BLUEs_PHP02_NS, BLUEs_PHP02_NS$REFYLD18!="-999")
  #BLUEs_PHP02_NS <- BLUEs_PHP02_NS[sample(nrow(BLUEs_PHP02_NS),42), ]
  flds <- createFolds(BLUEs_PHP02_NS$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
  BLUEs_PHP02_NS <- cbind(BLUEs_PHP02_NS, flds)
  
  BLUEs_PHP02_SS   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="SS")
  BLUEs_PHP02_SS <- subset(BLUEs_PHP02_SS, BLUEs_PHP02_SS$REFYLD18!="-999")
  #BLUEs_PHP02_SS <- BLUEs_PHP02_SS[sample(nrow(BLUEs_PHP02_SS),42), ]
  flds <- createFolds(BLUEs_PHP02_SS$PHTYLD, k=5, list=FALSE, returnTrain = FALSE)
  BLUEs_PHP02_SS <- cbind(BLUEs_PHP02_SS, flds)
  #Phenotype File Prepared--------------------------------------------------------------   
  BLUEs_PHP02 <- rbind(BLUEs_PHP02_IO, BLUEs_PHP02_NS, BLUEs_PHP02_SS)
  
  BLUEs_PHP02 <- BLUEs_PHP02 %>%
    arrange(BLUEs_PHP02$Grin1)
  
  #-------------------------------------------------------------------------------------   
  
  Grin1 <- as.data.frame(BLUEs_PHP02$Grin1)
  Grin1 <- Grin1 %>%
    rename(Grin1 = `BLUEs_PHP02$Grin1`)
  
  #Genotype File Preparation------------------------------------------------------------   
  #Add-----
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

  
  #-------------------------------------------------------------------------------------
  
  for(j in 1:26){
    x=colnum[j]  #set the current [i] column as "x"
    trait=colnames(BLUEs_PHP02)[x] #sets the current column header as the trait name
    pheno <- BLUEs_PHP02 #make a copy that we can use for our analysis
    colnames(pheno)[x]="y"  #renames the trait variable as "y" for the model analysis below
    
    for(i in 1:5){
      
      pheno1 <- pheno
      pheno1$y <- scale(pheno$y, center = TRUE, scale=TRUE)
      pheno1$y[pheno1$flds == i] <- NA
      
      fm=BGLR(y=pheno1$y,
                ETA=list(G=list(V=EVD_PHP02_Add_VanRaden$vectors, 
                       d=EVD_PHP02_Add_VanRaden$values,
                       model='RKHS')),
              nIter=100000,
              burnIn=25000,
              thin=100,
              saveAt='PHP02_TempScaled_A_',
              verbose=F)
      yhat1 <- fm$yHat
      data1 <- cbind(pheno,yhat1)
      Correlation <- cor(data1$yhat1[data1$flds==i], data1$y[data1$flds==i])
      Repeated <- n
      
      varE=scan( 'PHP02_TempScaled_A_varE.dat')
      varA=scan('PHP02_TempScaled_A_ETA_A_varU.dat')
      h2=(varA)/(varA+varE)
      h2=(varA)/(varA+varE)
      h2_mean=(mean(varA[251:1000]))/(mean(varA[251:1000])+mean(varE[251:1000]))
      
      data1 <- cbind(Repeated,trait,Correlation, h2_mean)
      DataOutput <- rbind(DataOutput, data1)
      
      data2 <- cbind(Repeated, trait, varE, varA, h2)
      DataOutput1 <- rbind(DataOutput1, data2)
      
      print(n)
      print(j)
      print(i)
    }
    
  }
  
}

write.csv(DataOutput, "PHP02_TempScaled_A_Single.Trait.Prediction_sat.csv")
write.csv(DataOutput1, 'PHP02_TempScaled_A_Variance.Comp_sat.csv')


save.image(file='PHP02_TempScaled_A_STGP_sat.RData')






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
  colnum=c(9:ncol(BLUEs_PHP02))
  
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
  
  
  #-------------------------------------------------------------------------------------
  for(j in 1:26){
    x=colnum[j]  #set the current [i] column as "x"
    trait=colnames(BLUEs_PHP02)[x] #sets the current column header as the trait name
    pheno <- BLUEs_PHP02 #make a copy that we can use for our analysis
    colnames(pheno)[x]="y"  #renames the trait variable as "y" for the model analysis below
    
    for(i in 1:5){
      
      pheno1 <- pheno
      pheno1$y <- scale(pheno$y, center = TRUE, scale=TRUE)
      pheno1$y[pheno1$flds == i] <- NA
      
      fm=BGLR(y=pheno1$y,
                ETA=list(G=list(V=EVD_PHP02_Add_VanRaden$vectors, 
                                d=EVD_PHP02_Add_VanRaden$values,
                                model='RKHS')),
              nIter=100000,
              burnIn=25000,
              thin=100,
              saveAt='PHP02_TropScaled_A_',
              verbose=T)
      yhat1 <- fm$yHat
      data1 <- cbind(pheno,yhat1)
      Correlation <- cor(data1$yhat1[data1$flds==i], data1$y[data1$flds==i])
      Repeated <- n
      
      varE=scan( 'PHP02_TropScaled_A_varE.dat')
      varA=scan('PHP02_TropScaled_A_ETA_A_varU.dat')
      #varD=scan('PHP02_FULL.SPLIT_AD_ETA_D_varU.dat')
      h2=(varA)/(varA+varE)
      h2=(varA)/(varA+varE)
      h2_mean=(mean(varA[251:1000]))/(mean(varA[251:1000])+mean(varE[251:1000]))
      
      data1 <- cbind(Repeated,trait,Correlation, h2_mean)
      DataOutput <- rbind(DataOutput, data1)
      
      data2 <- cbind(Repeated, trait, varE, varA, h2)
      DataOutput1 <- rbind(DataOutput1, data2)
      
      print(n)
      print(j)
      print(i)
    }
    
  }
  
}

write.csv(DataOutput, "PHP02_TropScaled_A_Single.Trait.Prediction_sat.csv")
write.csv(DataOutput1, 'PHP02_TropScaled_A_Variance.Comp_sat.csv')


save.image(file='PHP02_TropScaled_A_STGP_sat.RData')



