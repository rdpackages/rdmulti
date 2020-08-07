###################################################################
# rdmulti: analysis of RD designs with multiple cutoffs or scores
# RDMCPLOT illustration file
# 22-Apr-2020
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

rm(list = ls())
library(rdmulti)
library(ggplot2)

###################################################################
# Load data and generate plot variables (omitting plot)
###################################################################

data <- read.csv('simdata_multic.csv')
Y <- data$y
X <- data$x
C <- data$c

aux <- rdmcplot(Y,X,C,ci=95,nodraw=TRUE)

Xmean <- aux$Xmean
Ymean <- aux$Ymean
X0 <- aux$X0
X1 <- aux$X1
Yhat0 <- aux$Yhat0
Yhat1 <- aux$Yhat1
CI_l <- aux$CI_l
CI_r <- aux$CI_r

###################################################################
# Replicate default plot from rdmcplot output directly
###################################################################

aux$rdmc_plot

###################################################################
# Replicate default plot from rdmcplot by hand
###################################################################

# First cutoff: C=33

rdmc_plot <- ggplot() + theme_bw() + labs(x="Running variable", y="Outcome")
rdmc_plot <- rdmc_plot + geom_point(aes(x=Xmean[,1],y=Ymean[,1]),col="darkblue",shape=1,na.rm=TRUE) + 
                         geom_line(aes(x=X0[,1],y=Yhat0[,1]),col="darkblue",linetype=1,na.rm=TRUE) +
                         geom_line(aes(x=X1[,1],y=Yhat1[,1]),col="darkblue",linetype=1,na.rm=TRUE) +
                         geom_vline(xintercept=33,col="darkblue",linetype=1)
rdmc_plot

# Add second cutoff: C=66

rdmc_plot <- rdmc_plot + geom_point(aes(x=Xmean[,2],y=Ymean[,2]),col="darkred",shape=1,na.rm=TRUE) + 
                         geom_line(aes(x=X0[,2],y=Yhat0[,2]),col="darkred",linetype=1,na.rm=TRUE) +
                         geom_line(aes(x=X1[,2],y=Yhat1[,2]),col="darkred",linetype=1,na.rm=TRUE) +
                         geom_vline(xintercept=66,col="darkred",linetype=1)
rdmc_plot


# Add confidence intervals for both cutoffs

rdmc_plot <- rdmc_plot + geom_errorbar(aes(x=Xmean[,1], ymin=CI_l[,1], ymax=CI_r[,1]),col="darkblue",linetype=1) +
                         geom_errorbar(aes(x=Xmean[,2], ymin=CI_l[,2], ymax=CI_r[,2]),col="darkred",linetype=1)
rdmc_plot
