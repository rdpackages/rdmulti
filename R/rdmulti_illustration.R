###################################################################
# rdmulti: analysis of RD designs with multiple cutoffs or scores
# Illustration file
# 23-Aug-2020
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

rm(list = ls())

library(rdmulti)


###################################################################
# Multiple cutoffs
###################################################################

data <- read.csv('simdata_multic.csv')
Y <- data$y
X <- data$x
C <- data$c

aux <- rdmc(Y,X,C)
aux <- rdmc(Y,X,C,pooled_opt=paste('h=20','p=2',sep=','),verbose=TRUE)
aux <- rdmc(Y,X,C,h=c(11,10))
aux <- rdmc(Y,X,C,bwselect=c('msetwo','certwo'))


###################################################################
# Plots
###################################################################

aux <- rdmcplot(Y,X,C)
aux <- rdmcplot(Y,X,C,nobins=TRUE)
aux <- rdmcplot(Y,X,C,h=c(11,12),pvec=c(1,1))


###################################################################
# Cumulative cutoffs
###################################################################

data <- read.csv('simdata_cumul.csv')
Y <- data$y
X <- data$x
cvec <- c(data$c[1],data$c[2])

aux <- rdms(Y,X,cvec)
aux <- rdms(Y,X,cvec,h=c(11,8),kernel=c('uniform','triangular'))
aux <- rdms(Y,X,cvec,range = matrix(c(0,33.5,65.5,100),ncol=2))
cutoff <- cvec[1]*(X<=49.5) + cvec[2]*(X>49.5)
aux <- rdmc(Y,X,cutoff)

aux <- rdmcplot(Y,X,cutoff)


###################################################################
# Bivariate score
###################################################################

data <- read.csv('simdata_multis.csv')
Y <- data$y
X1 <- data$x1
X2 <- data$x2
zvar <- data$t
cvec <- c(data$c1[1],data$c1[2],data$c1[3])
cvec2 <- c(data$c2[1],data$c2[2],data$c2[3])

aux <- rdms(Y,X1,cvec,X2,zvar,cvec2)
aux <- rdms(Y,X1,cvec,X2,zvar,cvec2,h=c(15,13,17))
xnorm <- apply(cbind(abs(50-X1),abs(50-X2)),1,min)*(2*zvar-1)
aux <- rdms(Y,X1,cvec,X2,zvar,cvec2,xnorm=xnorm)

