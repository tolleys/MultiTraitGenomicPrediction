library("readxl")
library("ggplot2")

setwd("~/Documents/Purdue/PhD/EarPhotometry/BGLR_Prediction/DATA/")
results <- read.csv("2FACC_SRT_A_Single.Trait.Prediction_sat.csv", header=TRUE, sep=',')



cors_each.repeat <- results %>%
  group_by(trait, Repeated, Het_training, Het_testing) %>%
  summarize(cor_mean = mean(Correlation, na.rm=TRUE))



cors <- cors_each.repeat %>%
  group_by(trait, Het_training, Het_testing) %>%
  summarize(mean = mean(cor_mean, na.rm=TRUE),
            SE = sd(cor_mean)/sqrt(10))


#write.csv(cors, "~/Desktop/2FACC_SRT.Correlation_sat.csv")




ggplot(data=cors, aes(x=as.factor(trait), y=mean)) +
  geom_bar(stat="identity", fill="palegreen3",col="black", alpha=.6)+
  ylim(c(0,1)) +
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=0.2, size=1, color="black") +
  theme_bw() +
  coord_flip() +
  #theme(axis.text.x= element_text(angle=90)) +
  labs(x="Traits",
       y="Cross Validation Correlation")


ggplot(data=cors, aes(x=as.factor(trait), y=mean)) +
  geom_bar(stat="identity", fill="dodgerblue3",col="black", alpha=.6)+
  ylim(c(0,1)) +
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=0.2, size=1, color="black") +
  theme_bw() +
  coord_flip() +
  #theme(axis.text.x= element_text(angle=90)) +
  labs(x="Traits",
       y="Cross Validation Correlation")








Single_Het <- read_excel("~/Desktop/BGLR_TraitPrediction_sat.xlsx", sheet=4)

head(Single_Het)


Single_Het_PHP02 <- subset(Single_Het, Single_Het$Pedigree1=="PHP02")
Single_Het_PHP02 <- subset(Single_Het_PHP02, Single_Het_PHP02$logical37=="1")


ggplot(data=Single_Het_PHP02, aes(y=as.factor(TraitNum_Reverse), x=train_test, fill=Mean)) + 
  geom_tile() +
  geom_text(aes(label = Mean_round), color = "white", size = 4) +
  theme_bw() +
  scale_fill_gradient2(low="firebrick",mid="cornflowerblue", high="darkblue", midpoint=.1,
                       limits=c(-.4, 0.61)) +
  theme(axis.text.x = element_text(angle=20, vjust=.5)) +
  scale_y_discrete(labels=c("Average",   "TKERAB",  "SCTTER", "PHTKR",
                            "PHTKPR",  "PHTKPE",   "KERWTH",  "KERWGT", "KERPER",
                            "KERMIND", "KERMEAND", "KERMAXD", "KERLEN", "KERFIL",
                            "KERCC",  "KERARE",    "ETB",     "EARWTH", "EARVOL", 
                            "EARTR",  "EARPER",    "EARLGT",  "EARCW",  "EARBOX", 
                            "EARAREA", "PHYTYLD",   "REFYLD18")) + 
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text = element_text(size=14)) +
  theme(plot.title = element_text(size=16, face="bold", hjust=.5)) +
  theme(legend.position = "bottom") +
  labs(x = "Training_Testing Heterotic Group",
       y = "Trait",
       fill="Correlation")






Single_Het_2FACC <- subset(Single_Het, Single_Het$Pedigree1=="2FACC")
Single_Het_2FACC <- subset(Single_Het_2FACC, Single_Het_2FACC$logical37=="1")

ggplot(data=Single_Het_2FACC, aes(y=as.factor(TraitNum_Reverse), x=train_test, fill=Mean)) + 
  geom_tile() +
  geom_text(aes(label = Mean_round), color = "white", size = 4) +
  theme_bw() +
  scale_fill_gradient2(low="firebrick",mid="cornflowerblue", high="darkblue", midpoint=.1,
                       limits=c(-.4, 0.61)) +
  theme(axis.text.x = element_text(angle=20, vjust=.5)) +
  scale_y_discrete(position="right", labels=c("Average", "TKERAB",  "SCTTER", "PHTKR",
                                            "PHTKPR",  "PHTKPE",   "KERWTH",  "KERWGT", "KERPER",
                                            "KERMIND", "KERMEAND", "KERMAXD", "KERLEN", "KERFIL",
                                            "KERCC",  "KERARE",    "ETB",     "EARWTH", "EARVOL", 
                                            "EARTR",  "EARPER",    "EARLGT",  "EARCW",  "EARBOX", 
                                            "EARAREA", "PHTYLD", "REFYLD18")) + 
  theme(axis.title = element_text(size=16, face="bold")) +
  theme(axis.text = element_text(size=14)) +
  theme(plot.title = element_text(size=16, face="bold", hjust=.5)) +
  theme(legend.position = "bottom") +
  labs(x = "Training_Testing Heterotic Group",
       y = "Trait",
       fill="Correlation")


