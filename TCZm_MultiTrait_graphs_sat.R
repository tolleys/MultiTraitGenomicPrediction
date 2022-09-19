library("readxl")
library('ggplot2')

mtgp <- read_excel("~/Documents/Purdue/PhD/EarPhotometry/BGLR_Prediction/BGLR_TraitPrediction_sat.xlsx", sheet=6)[,c(1:9)]
head(mtgp)



mtgp_2FACC_Temp <- subset(mtgp, mtgp$Population=="2FACC_Temp")

  mtgp_2FACC_Temp_yld18 <- subset(mtgp_2FACC_Temp, mtgp_2FACC_Temp$Trait=="yld18")
  mtgp_2FACC_Temp_PHTYLD <- subset(mtgp_2FACC_Temp, mtgp_2FACC_Temp$Trait=="PHTYLD")

  
mtgp_PHP02_Temp <- subset(mtgp, mtgp$Population=="PHP02_Temp")
  
  mtgp_PHP02_Temp_yld18 <- subset(mtgp_PHP02_Temp, mtgp_PHP02_Temp$Trait=="yld18")
  mtgp_PHP02_Temp_PHTYLD <- subset(mtgp_PHP02_Temp, mtgp_PHP02_Temp$Trait=="PHTYLD")

  
mtgp_PHP02_Trop <- subset(mtgp, mtgp$Population=="PHP02_Trop")
  
  mtgp_PHP02_Trop_yld18 <- subset(mtgp_PHP02_Trop, mtgp_PHP02_Trop$Trait=="yld18")
  mtgp_PHP02_Trop_PHTYLD <- subset(mtgp_PHP02_Trop, mtgp_PHP02_Trop$Trait=="PHTYLD")

  
  
  
  
ggplot(data=mtgp_2FACC_Temp_PHTYLD, aes(x=MultiTrait_Num2, y=Mean, fill=Validation)) +
    geom_bar(stat='identity', width=.7, position=position_dodge(width=0.8), colour="black", alpha=.8) +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD),  
                  position=position_dodge(width=0.8),
                  width=0.2, size=.4, col='black') +
    #scale_fill_manual(values=c('firebrick', 'darkgreen', 'dodgerblue3')) +
    scale_fill_manual(values=c('skyblue', 'dodgerblue4', 'dodgerblue3'), 
                      labels=c('Unknown', 'Known', "STGP")) +
    ylim(c(0,1)) +
    theme_bw() + 
    theme(axis.title.x=element_blank(),
          axis.text.x =element_blank(),
          axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(),
          axis.text.y =element_blank(),
          axis.ticks.y=element_blank()) +
    theme(axis.title = element_text(size=16, face="bold")) +
    theme(axis.text = element_text(size=14)) +
    #guides(fill='none') +
    labs(y='Prediction Accuracy (r)',
         fill="Secondary Traits")

  


  