
###
# This code loads the testData.csv file and analyses with a 3rd degree polynomial, 
# saves a mag and phase plot, and then prints out the p-values for the model
###

# install libraries and download the data
library(gsignal)
library(ggplot2)
library(gghighlight)
library(ggsignif)
library(grid)
library(gridExtra)
library(ggpubr)
library(plotrix)
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(readxl)
library(broom.mixed)
library(MuMIn)
library(R.matlab)
library(viridis)


###
# Functions
### 

# calculate standard error
stderror <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x))

# Analyze and plot magnitude data
plotMagFig<-function(dTablePred,dTable,degree,colormap,dTave) {
  ylimMag=c(20*log10(0.5),20*log10(300))
  xlimMag=c(1,15)
  nLevel=length(unique(dTable$level))
  cohorts=as.character(unique(dTable$cohort))
  if (nLevel==1) {
    fit<-lmer(magLog ~ cohort * stats::poly(freq,degree) + (1|id), data=dTable)
    eqn=sprintf("log(mag) ~ poly(freq,%d) * cohort + (1|id)",degree)
  } else {
    fit<-lmer(magLog ~ cohort*stats::poly(freq,degree)*level + (1|id), data=dTable)
    eqn=sprintf("log(mag) ~ poly(freq,%d) * level * cohort + (1|id)",degree)
  }
  pSum<-summary(fit)$coefficients
  R2<-r.squaredGLMM(fit)
  
  R2text=sprintf("R^2 = %3.2f",R2[2])
  grob1 <- grobTree(textGrob(eqn, x=0.99,  y=0.95, just='right', gp=gpar(col="black", fontsize=7, fontface="italic")))
  grob2 <- grobTree(textGrob(R2text, x=0.99,  y=0.85, just='right', gp=gpar(col="black", fontsize=7, fontface="italic")))
  grob3 <- grobTree(textGrob(cohorts[1], x=0.99,  y=0.95, just='right', gp=gpar(col=colormap[1], fontsize=7, fontface="italic")))
  grob4 <- grobTree(textGrob(cohorts[2], x=0.99,  y=0.85, just='right', gp=gpar(col=colormap[2], fontsize=7, fontface="italic")))
  
  dT=cbind(dTablePred, new.predict = predict(fit,newdata=dTablePred,allow.new.levels = TRUE))
  p1pred<-ggplot()+ 
    geom_point(data=dTave,aes(x=freq,y=magLog, group=interaction(cohort,level), color=cohort), size=0.7) +
    #    geom_errorbar(data=dTave,aes(x=freq,ymin=magLog-magSEM, ymax=magLog+magSEM, color=cohort), width=0.01) +
    geom_line(data = dT, aes(x = freq, y = new.predict, group=interaction(id,level), color=cohort), size=0.7) +
    scale_colour_manual(values = colormap) +
    annotation_custom(grob1) +
    annotation_custom(grob2) +
    theme_bw() +
    labs(x = "Freq (kHz)", y="Mag (dB re:1 nm)") +
    scale_x_continuous(trans='log10',breaks=c(1,2,5,10),limits=xlimMag) +
    ylim(ylimMag) +
    theme(legend.position = "none") +
    theme(
      axis.title=element_text(size=8),
      axis.text=element_text(size=8)
    )
  
  p1res<-ggplot(augment(fit), aes(x = .fitted, y = .resid, color=cohort)) +
    geom_point(size=0.7) +
    scale_colour_manual(values = colormap)+
    theme_bw() +
    labs(x = "Fitted Mag (dB re:1 nm)", y="Residual") +
    xlim(ylimMag) +
    annotation_custom(grob3) +
    annotation_custom(grob4) +
    geom_hline(yintercept = 0) +
    theme(legend.position = "none") +
    theme(
      axis.title=element_text(size=8),
      axis.text=element_text(size=8)
    )
  
  p<-ggarrange(p1pred, p1res, ncol = 2, nrow = 1)
  dReturn<-list(p=p,pSum=pSum,R2=R2,eqn=eqn)
  return(dReturn)
}

plotPhFig<-function(dTablePred,dTable,degree,colormap,dTave) {
  ylimPh=c(-5,0.5)
  xlimMag=c(1,15)
  nLevel=length(unique(dTable$level))
  cohorts=as.character(unique(dTable$cohort))
  fit<-lmer(phaseCycle ~ cohort*stats::poly(freq,degree)*level + (1|id), data=dTable)
  eqn=sprintf("phase ~ poly(freq,%d) * level * cohort + (1|id)",degree)
  pSum<-summary(fit)$coefficients
  R2<-r.squaredGLMM(fit)
  
  R2text=sprintf("R^2=^%3.2f",R2[2])
  grob1 <- grobTree(textGrob(eqn, x=0.99,  y=0.95, just='right', gp=gpar(col="black", fontsize=7, fontface="italic")))
  grob2 <- grobTree(textGrob(R2text, x=0.99,  y=0.85, just='right', gp=gpar(col="black", fontsize=7, fontface="italic")))
  grob3 <- grobTree(textGrob(cohorts[1], x=0.99,  y=0.95, just='right', gp=gpar(col=colormap[1], fontsize=7, fontface="italic")))
  grob4 <- grobTree(textGrob(cohorts[2], x=0.99,  y=0.85, just='right', gp=gpar(col=colormap[2], fontsize=7, fontface="italic")))
  
  dT=cbind(dTablePred, new.predict = predict(fit,newdata=dTablePred,allow.new.levels = TRUE))
  p1pred<-ggplot()+ 
    geom_point(data=dTave,aes(x=freq,y=phaseCycle, group=interaction(cohort,level), color=cohort), size=0.7) +
    #    geom_errorbar(data=dTave,aes(x=freq,ymin=phaseCycle-phaseSEM, ymax=phaseCycle+phaseSEM, color=cohort), width=0.01) +
    geom_line(data = dT, aes(x = freq, y = new.predict, group=interaction(id,level), color=cohort), size=0.7) +
    scale_colour_manual(values = colormap) +
    annotation_custom(grob1) +
    annotation_custom(grob2) +
    theme_bw() +
    labs(x = "Freq (kHz)", y="Phase (cycles)") +
    ylim(ylimPh) +
    scale_x_continuous(trans='log10',breaks=c(1,2,5,10),limits=xlimMag) +
    theme(legend.position = "none") +
    theme(
      axis.title=element_text(size=8),
      axis.text=element_text(size=8)
    )
  
  p1res<-ggplot(augment(fit), aes(x = .fitted, y = .resid, color=cohort)) +
    geom_point(size=0.7) +
    scale_colour_manual(values = colormap)+
    theme_bw() +
    labs(x = "Fitted Phase (cycles)", y="Residual") +
    xlim(ylimPh) +
    annotation_custom(grob3) +
    annotation_custom(grob4) +
    geom_hline(yintercept = 0) +
    theme(legend.position = "none") +
    theme(
      axis.title=element_text(size=8),
      axis.text=element_text(size=8)
    )
  
  p<-ggarrange(p1pred, p1res, ncol = 2, nrow = 1)
  dReturn<-list(p=p,pSum=pSum,R2=R2,eqn=eqn)
  return(dReturn)
}

# Average tuning curves in each cohort
aveDataTable<-function(dTable) {
  freqs=sort(unique(dTable$freq))
  levels=sort(unique(dTable$level))
  cohorts=sort(unique(dTable$cohort))
  nid1=length(unique(dTable$id[dTable$cohort==cohorts[1]]))
  nid2=length(unique(dTable$id[dTable$cohort==cohorts[2]]))
  nFreq=length(freqs)
  nLevel=length(levels)
  ncohort=length(cohorts)
  
  dTave=tibble(
    cohort=factor(levels=cohorts),
    freq=numeric(),
    level=numeric(),
    magLog=double(),
    phaseCycle=numeric(),
    magSEM=double(),
    phaseSEM=numeric(),
  )
  
  for (g in 1:ncohort){
    for (f in 1:nFreq) {
      for (l in 1:nLevel) {
        dT<-filter(dTable, cohort==cohorts[g] & freq==freqs[f] & level==levels[l])
        dTave<-add_row(dTave,
                       cohort=cohorts[g],
                       freq=freqs[f],
                       level=levels[l],
                       magLog=20*log10(mean(dT$mag,na.rm=TRUE)),
                       phaseCycle=mean(dT$phaseCycle,na.rm=TRUE),
                       magSEM=stderror(20*log10(dT$mag)),
                       phaseSEM=stderror(dT$phaseCycle)
        )
      }
    }
  }
  return(dTave)
}

# create table for predicted values
predDataTable<-function(dTave) {
  freqs=seq(min(dTave$freq),max(dTave$freq),length.out=20)
  levels=sort(unique(dTave$level))
  cohorts=sort(unique(dTave$cohort))
  nFreq=length(freqs)
  nLevel=length(levels)
  nCohort=length(cohorts)  

  dTablePred=tibble(
    id=character(),
    cohort=factor(levels=cohorts),
    freq=numeric(),
    level=numeric(),
    mag=double(),
    phase=numeric(),
  )

  for (g in 1:nCohort) {
    id=sprintf('Id%d',g)
    for (f in 1:nFreq){
      for (l in 1:nLevel){
        dTablePred<-add_row(dTablePred,
                             id=id,
                             cohort=cohorts[g],
                             freq=freqs[f],
                             level=levels[l],
                             mag=0,
                             phase=0
        )
      }
    }
  }
  return(dTablePred)
}


###
# main program starts here
###


dTable<-read_csv('testData.csv') # Sample data in CSV file
dTable$magLog<-20*log10(dTable$mag)
dTable$phaseCycle<-dTable$phase/(2*pi)

dTable<-dTable[!is.na(dTable$mag),]
dTave<-aveDataTable(dTable)
dTablePred<-predDataTable(dTave)
colormap=c(viridis(4)[1],viridis(4)[2])


degree=3
filename="Figure.pdf"
dSet1<-plotMagFig(dTablePred,dTable,degree,colormap,dTave) 
dSet2<-plotPhFig(dTablePred,dTable,degree,colormap,dTave) 
p<-ggarrange(dSet1$p,dSet2$p,ncol = 1, nrow = 2)
ggsave(filename, plot = p, device = NULL, path = NULL,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 1200, limitsize = TRUE)
print(dSet1$pSum) # magnitude data: all coefficients and p-values
print(dSet2$pSum) # phase data: all coefficients and p-values
print('')
print('')

print(sprintf('Magnitude data cohort comparison p-value: %5.5f',dSet1$pSum[2,5]))
print(sprintf('Phase data cohort comparison p-value: %5.5f',dSet2$pSum[2,5]))
