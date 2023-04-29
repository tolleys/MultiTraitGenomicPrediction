# FILENAME: TCZm_PHP02_Even.Het_A.BGLR.Genomic.Prediction_Server.R
setwd("~/Documents/Purdue/PhD/Research/Chp2_MTGP/DATA/")
########Loading libraries###########
#install.packages(c('caret', 'BGLR', 'AGHmatrix','readxl',
#                   'bWGR', 'dplyr'),repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
library('caret', verbose=FALSE)
library("BGLR", verbose=FALSE)
library("AGHmatrix", verbose=FALSE)
library("readxl", verbose=FALSE)
library("bWGR", verbose=FALSE)
library("dplyr", verbose=FALSE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###############################[PHP02 Heterotic Groups]####################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
DataOutput <- data.frame()
DataOutput1 <- data.frame()

 for(n in 1:10){ 
   set.seed(24+n)

  for(m in 1:4){
    
    for(k in 1:4){
      BLUEs <- as.data.frame(read_excel("1.TCZm17_18_EP.BLUEs.xlsx", sheet=1))
      BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")
      colnum=c(9:ncol(BLUEs_PHP02))
      BLUEs_PHP02$REFYLD18 <- as.numeric(BLUEs_PHP02$REFYLD18)
      
      if( m == k){
        print("Don't run Script")
      }
      if(m!=k){
                
    BLUEs_PHP02_1 <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1Num==m)
    BLUEs_PHP02_1 <- subset(BLUEs_PHP02_1, BLUEs_PHP02_1$REFYLD18!="NA")
    BLUEs_PHP02_1 <- BLUEs_PHP02_1[sample(nrow(BLUEs_PHP02_1),37), ] #<----Can be blocked out for Supplemental Material Table
    
    
    BLUEs_PHP02_2   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1Num==k)
    BLUEs_PHP02_2   <- subset(BLUEs_PHP02_2, BLUEs_PHP02_2$REFYLD18 != "NA")
#Phenotype File Prepared--------------------------------------------------------------   
BLUEs_PHP02 <- rbind(BLUEs_PHP02_1, BLUEs_PHP02_2)
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
      
        
        pheno1 <- pheno
        pheno1$y <- scale(pheno$y, center = TRUE, scale=TRUE)
        pheno1$y[pheno1$Heterotic1Num == k] <- NA


        fm=BGLR(y=pheno1$y,
                ETA=list(A=list(V=EVD_PHP02_Add_VanRaden$vectors, 
                                d=EVD_PHP02_Add_VanRaden$values,
                                model='RKHS')),
                nIter=100000,
                burnIn=25000,
                thin=100,
                verbose=T,
                saveAt='PHP02_SingleHeteroticGroup_A_')
        yhat1 <- fm$yHat
        data1 <- cbind(pheno,yhat1)
        Correlation <- cor(data1$yhat1[data1$Heterotic1Num==k], data1$y[data1$Heterotic1Num==k])
        Repeated <- n
        
        varE=scan( 'PHP02_SingleHeteroticGroup_A_varE.dat')
        varA=scan('PHP02_SingleHeteroticGroup_A_ETA_A_varU.dat')

        h2=(varA)/(varA+varE)
        h2_mean=(mean(varA[251:1000]))/(mean(varA[251:1000])+mean(varE[251:1000]))
        
        data1 <- cbind(Repeated, m, k, trait, Correlation, h2_mean)
        DataOutput <- rbind(DataOutput, data1)
        
        data2 <- cbind(Repeated, m, k, trait, varE, varA, h2)
        DataOutput1 <- rbind(DataOutput1, data2)
        
        print(n)
        print(j)
        print(m)
        print(k)
          
        }
      }
    }
  }
}   
  



write.csv(DataOutput, "PHP02_SingleHeteroticGroup_A_Single.Trait.Prediction_sat.csv")
write.csv(DataOutput1, 'PHP02_SingleHeteroticGroup_A_Variance.Comp_sat.csv')


save.image(file='PHP02_SingleHeteroticGroup_A_STGP_sat.RData')



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###############################[2FACC Heterotic Groups]####################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

Grin1 <- read.csv("Grin1.csv", header=T)

DataOutput <- data.frame()
DataOutput1 <- data.frame()


for(n in 1:10){ 
  set.seed(24+n)
  
  for(m in 2:4){
    
    for(k in 2:4){
      
      BLUEs <- as.data.frame(read_excel("1.TCZm17_18_EP.BLUEs.xlsx", sheet=1))
      BLUEs_2FACC <- subset(BLUEs, BLUEs$Pedigree2=="2FACC")
      BLUEs_2FACC <- merge(Grin1, BLUEs_2FACC, by="Grin1")
      BLUEs_2FACC$REFYLD18 <- as.numeric(BLUEs_2FACC$REFYLD18)
      colnum=c(9:ncol(BLUEs_2FACC))
      
      if( m == k){
        print("Don't run Script")
      }
      if(m!=k){
        
        BLUEs_2FACC_1 <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1Num==m)
        BLUEs_2FACC_1 <- subset(BLUEs_2FACC_1, BLUEs_2FACC_1$REFYLD18!="NA")
        BLUEs_2FACC_1 <- BLUEs_2FACC_1[sample(nrow(BLUEs_2FACC_1),37), ] #<----Can be blocked out for Supplemental Material Table
        
        BLUEs_2FACC_2   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1Num==k)
        BLUEs_2FACC_2   <- subset(BLUEs_2FACC_2, BLUEs_2FACC_2$REFYLD18 != "NA")
        #Phenotype File Prepared--------------------------------------------------------------   
        BLUEs_2FACC <- rbind(BLUEs_2FACC_1, BLUEs_2FACC_2)
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
          pheno <- BLUEs_2FACC #make a copy that we can use for our analysis
          colnames(pheno)[x]="y"  #renames the trait variable as "y" for the model analysis below
          
          
          pheno1 <- pheno
          pheno1$y <- scale(pheno$y, center = TRUE, scale=TRUE)
          pheno1$y[pheno1$Heterotic1Num == k] <- NA
          
          
          fm=BGLR(y=pheno1$y,
                  ETA=list(A=list(V=EVD_2FACC_Add_VanRaden$vectors, 
                                  d=EVD_2FACC_Add_VanRaden$values,
                                  model='RKHS')),
                  nIter=100000,
                  burnIn=25000,
                  thin=100,
                  verbose=T,
                  saveAt='2FACC_SingleHeteroticGroup_A_')
          yhat1 <- fm$yHat
          data1 <- cbind(pheno,yhat1)
          Correlation <- cor(data1$yhat1[data1$Heterotic1Num==k], data1$y[data1$Heterotic1Num==k])
          Repeated <- n
          
          varE=scan( '2FACC_SingleHeteroticGroup_A_varE.dat')
          varA=scan('2FACC_SingleHeteroticGroup_A_ETA_A_varU.dat')
          
          h2=(varA)/(varA+varE)
          h2_mean=(mean(varA[251:1000]))/(mean(varA[251:1000])+mean(varE[251:1000]))
          
          data1 <- cbind(Repeated, m, k, trait, Correlation, h2_mean)
          DataOutput <- rbind(DataOutput, data1)
          
          data2 <- cbind(Repeated, m, k, trait, varE, varA, h2)
          DataOutput1 <- rbind(DataOutput1, data2)
          
          print(n)
          print(j)
          print(m)
          print(k)
          
        }
      }
    }
  }
}   




write.csv(DataOutput, "2FACC_SingleHeteroticGroup_A_Single.Trait.Prediction_sat.csv")
write.csv(DataOutput1, '2FACC_SingleHeteroticGroup_A_Variance.Comp_sat.csv')


save.image(file='2FACC_SingleHeteroticGroup_A_STGP_sat.RData')


