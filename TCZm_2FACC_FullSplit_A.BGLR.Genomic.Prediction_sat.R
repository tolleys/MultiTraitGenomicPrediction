# FILENAME: TCZm_PHP02_Even.Het_A.BGLR.Genomic.Prediction_Server.R
setwd("~/Documents/Purdue/PhD/EarPhotometry/BGLR_Prediction/DATA/")
########Loading libraries###########
install.packages('caret',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
library('caret', verbose=FALSE)
install.packages('BGLR',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
library("BGLR", verbose=FALSE)
install.packages('AGHmatrix',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
library("AGHmatrix", verbose=FALSE)
install.packages('readxl',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
library("readxl", verbose=FALSE)
install.packages('bWGR',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
library("bWGR", verbose=FALSE)
install.packages('dplyr',repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
library("dplyr", verbose=FALSE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
####################################[Server Loop]##########################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

BLUEs <- as.data.frame(read_excel("EP_BLUEs_Full.xlsx", sheet=2))

#BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")

BLUEs_2FACC <- subset(BLUEs, BLUEs$Pedigree2=="2FACC")
BLUEs_2FACC <- merge(Grin1, BLUEs_2FACC, by="Grin1")

#This was to center and scale the data.
#head(BLUEs_2FACC)
#BLUEs_2FACC_258_1.7 <- BLUEs_2FACC_258[,c(1:7)]
#BLUEs_2FACC_258_8.39 <- scale(BLUEs_2FACC_258[,c(8:39)], center=TRUE, scale=TRUE)
#BLUEs_2FACC_258 <- cbind(BLUEs_2FACC_258.1.7, BLUEs_2FACC_258_8.39)



DataOutput <- data.frame()
DataOutput1 <- data.frame()
colnum=c(11:ncol(BLUEs_2FACC))


for(n in 1:10){
    set.seed(24+n)
  
  BLUEs <- as.data.frame(read_excel("EP_BLUEs_Full.xlsx", sheet=2))
  
  #BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")
  
  BLUEs_2FACC <- subset(BLUEs, BLUEs$Pedigree2=="2FACC")
  BLUEs_2FACC <- merge(Grin1, BLUEs_2FACC, by="Grin1")
  BLUEs_2FACC$yld18 <- as.numeric(BLUEs_2FACC$yld18)
      
    BLUEs_2FACC_IO   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="Iodent")
    BLUEs_2FACC_IO <- subset(BLUEs_2FACC_IO, BLUEs_2FACC_IO$yld18!="NA")
      #BLUEs_2FACC_IO <- BLUEs_2FACC_IO[sample(nrow(BLUEs_2FACC_IO),14), ]
      flds <- createFolds(BLUEs_2FACC_IO$PHTYLD_g, k=5, list=FALSE, returnTrain = FALSE)
      BLUEs_2FACC_IO <- cbind(BLUEs_2FACC_IO, flds)
      
    BLUEs_2FACC_NS   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="NS")
    BLUEs_2FACC_NS <- subset(BLUEs_2FACC_NS, BLUEs_2FACC_NS$yld18!="NA")
      #BLUEs_2FACC_NS <- BLUEs_2FACC_NS[sample(nrow(BLUEs_2FACC_NS),14), ]
      flds <- createFolds(BLUEs_2FACC_NS$PHTYLD_g, k=5, list=FALSE, returnTrain = FALSE)
      BLUEs_2FACC_NS <- cbind(BLUEs_2FACC_NS, flds)
      
    BLUEs_2FACC_SS   <- subset(BLUEs_2FACC, BLUEs_2FACC$Heterotic1=="SS")
    BLUEs_2FACC_SS <- subset(BLUEs_2FACC_SS, BLUEs_2FACC_SS$yld18!="NA")
      #BLUEs_2FACC_SS <- BLUEs_2FACC_SS[sample(nrow(BLUEs_2FACC_SS),14), ]
      flds <- createFolds(BLUEs_2FACC_SS$PHTYLD_g, k=5, list=FALSE, returnTrain = FALSE)
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
      
      for(i in 1:5){

        pheno1 <- pheno
        pheno1$y <- scale(pheno$y, center = TRUE, scale=TRUE)
        pheno1$y[pheno1$flds == i] <- NA


        fm=BGLR(y=pheno1$y,
                ETA=list(A=list(V=EVD_2FACC_Add_VanRaden$vectors,
                                d=EVD_2FACC_Add_VanRaden$values,
                                model='RKHS'),
                         D=list(V=EVD_2FACC_Dom$vectors, 
                                d=EVD_2FACC_Dom$values,
                                model='RKHS')),
                nIter=100000,
                burnIn=25000,
                thin=100,
                verbose=F,
                saveAt='2FACC_Scaled_AD_')

        yhat1 <- fm$yHat
        data1 <- cbind(pheno,yhat1)
        Correlation <- cor(data1$yhat1[data1$flds==i], data1$y[data1$flds==i])
    
        Repeated <- n

        varE=scan( '2FACC_Scaled_AD_varE.dat')
        varA=scan('2FACC_Scaled_AD_ETA_A_varU.dat')
        varD=scan('2FACC_Scaled_AD_ETA_D_varU.dat')
        h2=(varA)/(varA+varE)
        h2_mean=(mean(varA[251:1000]))/(mean(varA[251:1000])+mean(varD[251:1000])+mean(varE[251:1000]))
        
        data1 <- cbind(Repeated,i,trait,Correlation, h2_mean)
        DataOutput <- rbind(DataOutput, data1)
        
        data2 <- cbind(Repeated, i, trait, varE, varA, varD, h2)
        DataOutput1 <- rbind(DataOutput1, data2)
        
        print(n)
        print(j)
        print(i)
      }
      
    }
    
}


write.csv(DataOutput, "2FACC_Scaled_A_Single.Trait.Prediction_sat.csv")
write.csv(DataOutput1, '2FACC_Scaled_A_Variance.Comp_sat.csv')


save.image(file='2FACC_Scaled_A_STGP_sat.RData')



