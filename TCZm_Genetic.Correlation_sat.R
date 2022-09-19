# FILENAME: TCZm_PHP02_Even.Het_A.BGLR.Genomic.Prediction_Server.R
setwd("~/Documents/Purdue/PhD/EarPhotometry/BGLR_Prediction/DATA/")
########Loading libraries###########
install.packages('caret', repos = "http://ftp.ussg.iu.edu/CRAN/", lib="/scratch/bell/tolleys/R_LIBS/")
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
library("corrplot")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###########################################[PHP02_ALL]#################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
BLUEs <- as.data.frame(read_excel("EP_BLUEs_Full.xlsx", sheet=2))

BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")


BLUEs_PHP02_DTMA   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="DTMA")
BLUEs_PHP02_DTMA <- subset(BLUEs_PHP02_DTMA, BLUEs_PHP02_DTMA$yld18!="NA")

BLUEs_PHP02_IO   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="Iodent")
BLUEs_PHP02_IO <- subset(BLUEs_PHP02_IO, BLUEs_PHP02_IO$yld18!="NA")

BLUEs_PHP02_NS   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="NS")
BLUEs_PHP02_NS <- subset(BLUEs_PHP02_NS, BLUEs_PHP02_NS$yld18!="NA")

BLUEs_PHP02_SS   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="SS")
BLUEs_PHP02_SS <- subset(BLUEs_PHP02_SS, BLUEs_PHP02_SS$yld18!="NA")


BLUEs_PHP02 <- rbind(BLUEs_PHP02_DTMA, BLUEs_PHP02_IO, BLUEs_PHP02_NS, BLUEs_PHP02_SS)
BLUEs_PHP02$yld18 <- as.numeric(BLUEs_PHP02$yld18)

BLUEs_PHP021 <- BLUEs_PHP02 %>%
  arrange(BLUEs_PHP02$Grin1)
head(BLUEs_PHP021)  
BLUEs_PHP021 <- BLUEs_PHP02[,-c(1:10)]
BLUEs_PHP021$yld18 <- as.numeric(BLUEs_PHP021$yld18)

BLUEs_PHP021 <- as.matrix(BLUEs_PHP021)



dim(BLUEs_PHP02)

rownames(BLUEs_PHP021) <- BLUEs_PHP02$Grin1
str(BLUEs_PHP021)
head(BLUEs_PHP021)

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

#E = eigen(PHP02_Add_VanRaden,T)
#w = 1:(which((cumsum(E$values)/sum(E$values))>=0.98)[1])
#E$values = E$values[w]
#E$vectors = E$vectors[,w]
#X = E$vectors %*% diag(sqrt(E$values))



BLUEs_PHP02_stand <- scale(BLUEs_PHP021, center=TRUE, scale=TRUE)

head(BLUEs_PHP02_stand)
#model1 = mrr(BLUEs_PHP02_stand, X)


model1 <- mkr(BLUEs_PHP02_stand, PHP02_Add_VanRaden)
gcor_PHP02_All <- model1$GC  

rownames(gcor_PHP02_All) <- c("yld18", "PHTYLD", "PHTKPE", "PHTKPR", "PHTKR", "KERFIL", "TKERAB",
                               "SCTTER", "EARLGT", "EARWTH", "EARCW", "EARAREA", "EARVOL", "EARPER",     
                               "EARBOX", "ETB", "EARTR", "KERCC", "KERARE", "KERWGT", "KERLEN",     
                               "KERWTH", "KERMEAND", "KERMIND", "KERMAXD", "KERPER")
colnames(gcor_PHP02_All) <- c("yld18", "PHTYLD", "PHTKPE", "PHTKPR", "PHTKR", "KERFIL", "TKERAB",
                               "SCTTER", "EARLGT", "EARWTH", "EARCW", "EARAREA", "EARVOL", "EARPER",     
                               "EARBOX", "ETB", "EARTR", "KERCC", "KERARE", "KERWGT", "KERLEN",     
                               "KERWTH", "KERMEAND", "KERMIND", "KERMAXD", "KERPER")

gcor_PHP02_All[,1:2]



corrplot(gcor_PHP02_All, type="lower", method="square", tl.col="black", order="original",
         bg="white", diag=F, tl.srt=45, addCoef.col = "black", addCoefasPercent = TRUE)








### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###########################################[PHP02_Temp]################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
BLUEs <- as.data.frame(read_excel("EP_BLUEs_Full.xlsx", sheet=2))

BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")


BLUEs_PHP02_IO   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="Iodent")
BLUEs_PHP02_IO <- subset(BLUEs_PHP02_IO, BLUEs_PHP02_IO$yld18!="NA")

BLUEs_PHP02_NS   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="NS")
BLUEs_PHP02_NS <- subset(BLUEs_PHP02_NS, BLUEs_PHP02_NS$yld18!="NA")

BLUEs_PHP02_SS   <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="SS")
BLUEs_PHP02_SS <- subset(BLUEs_PHP02_SS, BLUEs_PHP02_SS$yld18!="NA")

  
BLUEs_PHP02 <- rbind(BLUEs_PHP02_IO, BLUEs_PHP02_NS, BLUEs_PHP02_SS)
BLUEs_PHP02$yld18 <- as.numeric(BLUEs_PHP02$yld18)

BLUEs_PHP021 <- BLUEs_PHP02 %>%
  arrange(BLUEs_PHP02$Grin1)
head(BLUEs_PHP021)  
BLUEs_PHP021 <- BLUEs_PHP02[,-c(1:10)]
BLUEs_PHP021$yld18 <- as.numeric(BLUEs_PHP021$yld18)

BLUEs_PHP021 <- as.matrix(BLUEs_PHP021)



dim(BLUEs_PHP02)

rownames(BLUEs_PHP021) <- BLUEs_PHP02$Grin1
str(BLUEs_PHP021)
head(BLUEs_PHP021)

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

#E = eigen(PHP02_Add_VanRaden,T)
#w = 1:(which((cumsum(E$values)/sum(E$values))>=0.98)[1])
#E$values = E$values[w]
#E$vectors = E$vectors[,w]
#X = E$vectors %*% diag(sqrt(E$values))



BLUEs_PHP02_stand <- scale(BLUEs_PHP021, center=TRUE, scale=TRUE)

head(BLUEs_PHP02_stand)
#model1 = mrr(BLUEs_PHP02_stand, X)


model1 <- mkr(BLUEs_PHP02_stand, PHP02_Add_VanRaden)
gcor_PHP02_Temp <- model1$GC  

rownames(gcor_PHP02_Temp) <- c("yld18", "PHTYLD", "PHTKPE", "PHTKPR", "PHTKR", "KERFIL", "TKERAB",
                     "SCTTER", "EARLGT", "EARWTH", "EARCW", "EARAREA", "EARVOL", "EARPER",     
                     "EARBOX", "ETB", "EARTR", "KERCC", "KERARE", "KERWGT", "KERLEN",     
                     "KERWTH", "KERMEAND", "KERMIND", "KERMAXD", "KERPER")
colnames(gcor_PHP02_Temp) <- c("yld18", "PHTYLD", "PHTKPE", "PHTKPR", "PHTKR", "KERFIL", "TKERAB",
                     "SCTTER", "EARLGT", "EARWTH", "EARCW", "EARAREA", "EARVOL", "EARPER",     
                     "EARBOX", "ETB", "EARTR", "KERCC", "KERARE", "KERWGT", "KERLEN",     
                     "KERWTH", "KERMEAND", "KERMIND", "KERMAXD", "KERPER")

gcor_PHP02_Temp[,1:2]



corrplot(gcor_PHP02_Temp, type="lower", method="square", tl.col="black", order="original",
         bg="white", diag=F, tl.srt=45, addCoef.col = "black", addCoefasPercent = TRUE)










### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###########################################[PHP02_DTMA]################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
BLUEs <- as.data.frame(read_excel("EP_BLUEs_Full.xlsx", sheet=2))

BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")

BLUEs_PHP02_DTMA <- subset(BLUEs_PHP02, BLUEs_PHP02$Heterotic1=="DTMA")
BLUEs_PHP02_DTMA <- subset(BLUEs_PHP02_DTMA, BLUEs_PHP02_DTMA$yld18!="NA")

BLUEs_PHP02 <- rbind(BLUEs_PHP02_DTMA)
BLUEs_PHP02$yld18 <- as.numeric(BLUEs_PHP02$yld18)


BLUEs_PHP021 <- BLUEs_PHP02 %>%
  arrange(BLUEs_PHP02$Grin1)
head(BLUEs_PHP021)  
BLUEs_PHP021 <- BLUEs_PHP02[,-c(1:10)]
BLUEs_PHP021$yld18 <- as.numeric(BLUEs_PHP021$yld18)

BLUEs_PHP021 <- as.matrix(BLUEs_PHP021)



dim(BLUEs_PHP02)

rownames(BLUEs_PHP021) <- BLUEs_PHP02$Grin1
str(BLUEs_PHP021)
head(BLUEs_PHP021)


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

E = eigen(PHP02_Add_VanRaden,T)
w = 1:(which((cumsum(E$values)/sum(E$values))>=0.98)[1])
E$values = E$values[w]
E$vectors = E$vectors[,w]
X = E$vectors %*% diag(sqrt(E$values))


BLUEs_PHP02_stand <- scale(BLUEs_PHP021, center=TRUE, scale=TRUE)
head(BLUEs_PHP02_stand)
model1 = mrr(BLUEs_PHP02_stand, X)

gcor_PHP02_DTMA <- model1$GC  

rownames(gcor_PHP02_DTMA) <- c("yld18", "PHTYLD", "PHTKPE", "PHTKPR", "PHTKR", "KERFIL", "TKERAB",
                     "SCTTER", "EARLGT", "EARWTH", "EARCW", "EARAREA", "EARVOL", "EARPER",     
                     "EARBOX", "ETB", "EARTR", "KERCC", "KERARE", "KERWGT", "KERLEN",     
                     "KERWTH", "KERMEAND", "KERMIND", "KERMAXD", "KERPER")
colnames(gcor_PHP02_DTMA) <- c("yld18", "PHTYLD", "PHTKPE", "PHTKPR", "PHTKR", "KERFIL", "TKERAB",
                     "SCTTER", "EARLGT", "EARWTH", "EARCW", "EARAREA", "EARVOL", "EARPER",     
                     "EARBOX", "ETB", "EARTR", "KERCC", "KERARE", "KERWGT", "KERLEN",     
                     "KERWTH", "KERMEAND", "KERMIND", "KERMAXD", "KERPER")

gcor_PHP02_DTMA[,1:2]


corrplot(gcor_PHP02_DTMA, type="lower", method="square", tl.col="black", order="original",
         bg="white", diag=F, tl.srt=45, addCoef.col = "black", addCoefasPercent = TRUE)








### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###########################################[2FACC_Temp]################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
FACC2_Add <- read.table('2FACC_Add.txt', header=T, sep=" ")
Grin1_2FACC_Add <- FACC2_Add$Hybrids
Grin1_2FACC_Add <- gsub("\t", "", x=Grin1_2FACC_Add, perl = T)
Grin1_2FACC_Add <- gsub("0_", "", x=Grin1_2FACC_Add, perl=T)
FACC2_Add[1:5,1:5]
FACC2_Add <- FACC2_Add[,-c(1)]
rownames(FACC2_Add) <- Grin1_2FACC_Add
FACC2_Add$GRIN1 <- Grin1_2FACC_Add
FACC2_Add <- FACC2_Add[,c(52960,1:52959)]
FACC2_Add <- FACC2_Add[!(FACC2_Add$GRIN1=='Ames23447' | FACC2_Add$GRIN1=='Ames27067' |
                           FACC2_Add$GRIN1=='PI559919' | FACC2_Add$GRIN1=='PI601003' |
                           FACC2_Add$GRIN1=='PI601269' | FACC2_Add$GRIN1=='PI601558' ),]
Grin1 <- as.data.frame(FACC2_Add$GRIN1)
Grin1$Grin1 <- Grin1$`FACC2_Add$GRIN1`
Grin1 <- Grin1[,-1]
Grin1 <- as.data.frame(Grin1)
FACC2_Add <- FACC2_Add[,-c(1)]
FACC2_Add <- as.matrix(FACC2_Add)

dim(FACC2_Add)
FACC2_Add_VanRaden <- Gmatrix(SNPmatrix = FACC2_Add, missingValue=-9,
                              maf=0.05, method="VanRaden")

EVD_2FACC_Add_VanRaden <- eigen(FACC2_Add_VanRaden)
rownames(EVD_2FACC_Add_VanRaden$vectors) <- rownames(FACC2_Add_VanRaden)



E = eigen(FACC2_Add_VanRaden,T)
w = 1:(which((cumsum(E$values)/sum(E$values))>=0.98)[1])
E$values = E$values[w]
E$vectors = E$vectors[,w]
X = E$vectors %*% diag(sqrt(E$values))

BLUEs <- as.data.frame(read_excel("EP_BLUEs_Full.xlsx", sheet=2))

#BLUEs_PHP02 <- subset(BLUEs, BLUEs$Pedigree2=="PHP02")

BLUEs_2FACC <- subset(BLUEs, BLUEs$Pedigree2=="2FACC")
BLUEs_2FACC <- merge(Grin1, BLUEs_2FACC, by="Grin1")

head(BLUEs_2FACC)
BLUEs_2FACC1 <- BLUEs_2FACC[,-c(1:10,37)]
BLUEs_2FACC1$yld18 <- as.numeric(BLUEs_2FACC1$yld18)
BLUEs_2FACC1 <- as.matrix(BLUEs_2FACC1)
rownames(BLUEs_2FACC1) <- BLUEs_2FACC$Grin1


BLUEs_2FACC_stand <- scale(BLUEs_2FACC1, center=TRUE, scale=TRUE)
head(BLUEs_2FACC_stand)
model1 = mrr(BLUEs_2FACC_stand, X)



model1 <- mkr(BLUEs_2FACC_stand, FACC2_Add_VanRaden)
gcor_2FACC_Temp <- model1$GC  

#model1 <- mkr(BLUEs_2FACC1, FACC2_Add_VanRaden)


rownames(gcor_2FACC_Temp) <- c("yld18", "PHTYLD", "PHTKPE", "PHTKPR", "PHTKR", "KERFIL", "TKERAB",
                    "SCTTER", "EARLGT", "EARWTH", "EARCW", "EARAREA", "EARVOL", "EARPER",     
                    "EARBOX", "ETB", "EARTR", "KERCC", "KERARE", "KERWGT", "KERLEN",     
                    "KERWTH", "KERMEAND", "KERMIND", "KERMAXD", "KERPER")
colnames(gcor_2FACC_Temp) <- c("yld18", "PHTYLD", "PHTKPE", "PHTKPR", "PHTKR", "KERFIL", "TKERAB",
                    "SCTTER", "EARLGT", "EARWTH", "EARCW", "EARAREA", "EARVOL", "EARPER",     
                    "EARBOX", "ETB", "EARTR", "KERCC", "KERARE", "KERWGT", "KERLEN",     
                    "KERWTH", "KERMEAND", "KERMIND", "KERMAXD", "KERPER")


gcor_2FACC_Temp[,1:2]

corrplot(gcor_2FACC_Temp, type="lower", method="square", tl.col="black", order="original",
         bg="white", diag=F, tl.srt=45, addCoef.col = "black", addCoefasPercent = TRUE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#########################################[Yield Cor Heatmap]###########################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###




Gen_Cor <- read_excel("~/Documents/Purdue/PhD/EarPhotometry/BGLR_Prediction/BGLR_TraitPrediction_sat.xlsx", sheet=2)

head(Gen_Cor)
Gen_Cor1 <- subset(Gen_Cor, Gen_Cor$Population!="PHP02_All")

Gen_Cor1$Gcor_BLUPF90 <- as.numeric(Gen_Cor1$Gcor_BLUPF90)
Gen_Cor1$Gcor_BLUPF90_Rounded <- as.numeric(Gen_Cor1$Gcor_BLUPF90_Rounded)

ggplot(data=Gen_Cor1, aes(y=as.factor(TraitNum_Reverse), x=as.factor(Trait1Num_Correct), fill=Gcor_BLUPF90)) + 
  geom_tile() +
  geom_text(aes(label = Gcor_BLUPF90_Rounded), color = "white", size = 4) +
  theme_bw() +
  scale_fill_gradient2(low="firebrick",mid="cornflowerblue", high="darkblue",
                       limits=c(-1, 1)) +
  theme(axis.text.x = element_text(angle=0, vjust=.5)) +
  scale_y_discrete(labels=c("TKERAB",  "SCTTER", "PHTKR",
                            "PHTKPR",  "PHTKPE",   "KERWTH",  "KERWGT", "KERPER",
                            "KERMIND", "KERMEAND", "KERMAXD", "KERLEN", "KERFIL",
                            "KERCC",  "KERARE",    "ETB",     "EARWTH", "EARVOL", 
                            "EARTR",  "EARPER",    "EARLGT",  "EARCW",  "EARBOX", 
                            "EARAREA", "PHTYLD")) + 
  scale_x_discrete(labels=c("REFYLD18", "PHTYLD", "REFYLD18", "PHTYLD",
                            "REFYLD18", "PHTYLD", "REFYLD18", "PHTYLD")) + 
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text = element_text(size=14)) +
  theme(plot.title = element_text(size=16, face="bold", hjust=.5)) +
  theme(legend.position = "right") +
  labs(x = "",
       y = "Trait",
       fill="Correlation")







