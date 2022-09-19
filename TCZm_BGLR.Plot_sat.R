library("dplyr")
library("readxl")
library("ggplot2")

setwd("~/Documents/Purdue/PhD/EarPhotometry/BGLR_Prediction/DATA/")
results <- read.csv("2FACC_Scaled_A_Single.Trait.Prediction_sat.csv", header=TRUE, sep=',')



cors_each.repeat <- results %>%
  group_by(trait) %>%
  summarize(cor_mean = mean(Correlation, na.rm=TRUE),
            cor_SD = sd(Correlation))



#write.csv(cors_each.repeat, "~/Desktop/2FACC_TempScaled_A.Correlation_sat.csv")
    

        

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
