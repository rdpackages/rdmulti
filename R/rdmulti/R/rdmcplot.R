###################################################################
# rdmcplot: RD plots with multiple cutoffs
# !version 1.1 20-Jun-2023
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' RD plots with multiple cutoffs.
#'
#' \code{rdmcplot()} RD plots with multiple cutoffs.
#'
#'
#' @author
#' Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}
#'
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}
#'
#' Gonzalo Vazquez-Bare, UC Santa Barbara. \email{gvazquez@econ.ucsb.edu}
#'
#' @references
#'
#' Cattaneo, M.D., R. Titiunik and G. Vazquez-Bare. (2020). \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf}{Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores}. \emph{Stata Journal}, forthcoming.
#'
#'
#' @param Y outcome variable.
#' @param X running variable.
#' @param C cutoff variable.
#' @param nbinsmat matrix of cutoff-specific number of bins. See \code{rdplot()}
#'   for details.
#' @param binselectvec vector of cutoff-specific bins selection method. See
#'   \code{rdplot()} for details.
#' @param scalevec vector of cutoff-specific scale factors. See \code{rdplot()}
#'   for details.
#' @param supportmat matrix of cutoff-specific support conditions. See
#'   \code{rdplot()} for details..
#' @param pvec vector of cutoff-specific polynomial orders. See \code{rdplot()}
#'   for details.
#' @param hmat matrix of cutoff-specific bandwidths. See \code{rdplot()} for
#'   details.
#' @param kernelvec vector of cutoff-specific kernels. See \code{rdplot()} for
#'   details.
#' @param weightsvec vector of cutoff-specific weights. See \code{rdplot()} for
#'   details.
#' @param covs_mat matrix of covariates. See \code{rdplot()} for
#'   details.
#' @param covs_list list of of covariates to be used in each cutoff.
#' @param covs_evalvec vector indicating the evaluation point for additional
#'   covariates. See \code{rdrobust()} for details.
#' @param covs_dropvec vector indicating whether collinear covariates should be
#'   dropped at each cutoff. See \code{rdrobust()} for details.
#' @param ci adds confidence intervals of the specified level to the plot. See
#'   \code{rdrobust()} for details.
#' @param col_bins vector of colors for bins.
#' @param pch_bins vector of characters (pch) type for bins.
#' @param col_poly vector of colors for polynomial curves.
#' @param lty_poly vector of lty for polynomial curves.
#' @param col_xline vector of colors for vertical lines.
#' @param lty_xline vector of lty for vertical lines.
#' @param nobins omits bins plot.
#' @param nopoly omits polynomial curve plot.
#' @param noxline omits vertical lines indicating the cutoffs.
#' @param nodraw omits plot.
#'
#'
#' @return
#' \item{clist}{list of cutoffs}
#' \item{cnum}{number of cutoffs}
#' \item{X0}{matrix of X values for control units}
#' \item{X1}{matrix of X values for treated units}
#' \item{Yhat0}{estimated polynomial for control units}
#' \item{Yhat1}{estimated polynomial for treated units}
#' \item{Xmean}{bin average of X values}
#' \item{Ymean}{bin average for Y values}
#' \item{CI_l}{lower end of confidence intervals}
#' \item{CI_r}{upper end of confidence intervals}
#' \item{cfail}{Cutoffs where rdrobust() encountered problems}

#'
#'
#' @examples
#' # Toy dataset
#' X <- runif(1000,0,100)
#' C <- c(rep(33,500),rep(66,500))
#' Y <- (1 + X + (X>=C))*(C==33)+(.5 + .5*X + .8*(X>=C))*(C==66) + rnorm(1000)
#' # rdmcplot with standard syntax
#' tmp <- rdmcplot(Y,X,C)
#'
#'
#' @export

rdmcplot <- function(Y,X,C,nbinsmat=NULL,binselectvec=NULL,scalevec=NULL,
                     supportmat=NULL,pvec=NULL,hmat=NULL,kernelvec=NULL,
                     weightsvec=NULL,covs_mat=NULL,covs_list=NULL,covs_evalvec=NULL,
                     covs_dropvec=NULL,ci=NULL,col_bins=NULL,pch_bins=NULL,
                     col_poly=NULL,lty_poly=NULL,col_xline=NULL,lty_xline=NULL,
                     nobins=FALSE,nopoly=FALSE,noxline=FALSE,nodraw=FALSE){

  #################################################################
  # Setup and error checking
  #################################################################

  if (!is.numeric(C)){stop('C has to be numeric')}
  if (max(C,na.rm=TRUE)>=max(X,na.rm=TRUE) | min(C,na.rm=TRUE)<=min(X,na.rm=TRUE)){stop('cutoff variable outside range of running variable')}

  clist <- sort(unique(C))
  cnum <- length(clist)

  D <- as.numeric(X>=C)

  if (is.null(pvec)){pvec = rep(4,cnum)}
  if (!is.null(hmat)){
    if (is.null(dim(hmat))){
      hmat <- matrix(hmat,nrow=cnum,ncol=2)
    }
    haux <- hmat
  }  else{
    haux <- matrix(Inf,ncol=2,nrow=cnum)
  }
  if (!is.null(nbinsmat)){
    if (is.null(dim(nbinsmat))){
      nbinsmat <- matrix(nbinsmat,nrow=cnum,ncol=2)
    }
  }
  if (is.null(binselectvec)) binselectvec <- rep('esmv',cnum)
  if (is.null(scalevec)) scalevec <- rep(1,cnum)
  if (is.null(kernelvec)) kernelvec <- rep('uni',cnum)
  if (is.null(covs_evalvec)) covs_evalvec <- rep(0,cnum)
  if (is.null(covs_dropvec)) covs_dropvec <- rep(TRUE,cnum)
  if (!is.null(covs_mat)){
    covs_mat <- as.matrix(covs_mat)
    if (!is.null(covs_list)){
      if (length(covs_list)!=cnum) stop('Elements in covs_list should equal number of cutoffs')
    }
  }

  X0 <- matrix(NA,nrow=length(Y),ncol=cnum)
  X1 <- matrix(NA,nrow=length(Y),ncol=cnum)
  YHAT0 <- matrix(NA,nrow=length(Y),ncol=cnum)
  YHAT1 <- matrix(NA,nrow=length(Y),ncol=cnum)
  XMEAN <- matrix(NA,nrow=length(Y),ncol=cnum)
  YMEAN <- matrix(NA,nrow=length(Y),ncol=cnum)
  CI_l <- matrix(NA,nrow=length(Y),ncol=cnum)
  CI_r <- matrix(NA,nrow=length(Y),ncol=cnum)
  Cfail <- numeric()


  #################################################################
  # Construct variables for plots
  #################################################################

  count <- 1
  count_fail <- 0
  for (c in clist){

    yc <- Y[C==c & X<=c+haux[count,2] & X>=c-haux[count,1]]
    xc <- X[C==c & X<=c+haux[count,2] & X>=c-haux[count,1]]
    dc <- D[C==c & X<=c+haux[count,2] & X>=c-haux[count,1]]
    yc0 <- yc[dc==0]
    yc1 <- yc[dc==1]
    xc0 <- xc[dc==0]
    xc1 <- xc[dc==1]

    if (!is.null(covs_mat)){
      covs_mat_c <- covs_mat[C==c & X<=c+haux[count,2] & X>=c-haux[count,1],]
      if (!is.null(covs_list)){
        covs_aux <- covs_mat_c[,covs_list[[count]]]
      } else{
        covs_aux <- covs_mat_c
      }
    } else {
      covs_aux <- NULL
    }

    aux <- try(rdrobust::rdplot(yc,xc,c=c,
                                nbins=nbinsmat[count,],
                                binselect=binselectvec[count],
                                scale=scalevec[count],
                                support=supportmat[count,],
                                p=pvec[count],
                                h=hmat[count,],
                                kernel=kernelvec[count],
                                weights=weightsvec[count],
                                covs=covs_aux,
                                covs_eval=covs_evalvec[count],
                                covs_drop=covs_dropvec[count],
                                ci=ci,
                                hide=TRUE),
               silent=TRUE)

    #if (class(aux)!="try-error"){
    if (!inherits(aux,'try-error')){
      xmean <- aux$vars_bins[,2]
      ymean <- aux$vars_bins[,3]

      xmean <- xmean[!is.na(xmean)]
      ymean <- ymean[!is.na(ymean)]

      x0 <- aux$vars_poly[aux$vars_poly[,1]<c,1]
      yhat0 <- aux$vars_poly[aux$vars_poly[,1]<c,2]
      x1 <- aux$vars_poly[aux$vars_poly[,1]>c,1]
      yhat1 <- aux$vars_poly[aux$vars_poly[,1]>c,2]

      length(xmean) <- length(Y)
      length(ymean) <- length(Y)
      length(x0) <- length(Y)
      length(yhat0) <- length(Y)
      length(x1) <- length(Y)
      length(yhat1) <- length(Y)

      XMEAN[,count] <- xmean
      YMEAN[,count] <- ymean
      X0[,count] <- x0
      X1[,count] <- x1
      YHAT0[,count] <- yhat0
      YHAT1[,count] <- yhat1

      if(!is.null(ci)){
        ci_l <- aux$vars_bins[,8]
        ci_r <- aux$vars_bins[,9]
        length(ci_l) <- length(Y)
        length(ci_r) <- length(Y)
        CI_l[,count] <- ci_l
        CI_r[,count] <- ci_r
      }
    } else{
      Cfail <- c(Cfail,c)
      count_fail <- count_fail + 1
    }

    count <- count + 1

  }

  Xmean <- data.frame(XMEAN)
  Xmean <- Xmean[1:max(colSums(!is.na(Xmean))),]
  names(Xmean) <- paste0(rep("Xmean"),1:cnum)
  Ymean <- data.frame(YMEAN)
  Ymean <- Ymean[1:max(colSums(!is.na(Ymean))),]
  names(Ymean) <- paste0(rep("Ymean"),1:cnum)
  X0    <- data.frame(X0)
  X0    <- X0[1:max(colSums(!is.na(X0))),]
  names(X0) <- paste0(rep("X0_"),1:cnum)
  X1    <- data.frame(X1)
  X1    <- X1[1:max(colSums(!is.na(X1))),]
  names(X1) <- paste0(rep("X1_"),1:cnum)
  Yhat0    <- data.frame(YHAT0)
  Yhat0    <- Yhat0[1:max(colSums(!is.na(Yhat0))),]
  names(Yhat0) <- paste0(rep("Yhat0_"),1:cnum)
  Yhat1    <- data.frame(YHAT1)
  Yhat1    <- Yhat1[1:max(colSums(!is.na(Yhat1))),]
  names(Yhat1) <- paste0(rep("Yhat1_"),1:cnum)
  if(!is.null(ci)){
    CI_l <- data.frame(CI_l)
    CI_l <- CI_l[1:max(colSums(!is.na(CI_l))),]
    names(CI_l) <- paste0(rep("CI_l_"),1:cnum)
    CI_r <- data.frame(CI_r)
    CI_r <- CI_r[1:max(colSums(!is.na(CI_r))),]
    names(CI_r) <- paste0(rep("CI_r_"),1:cnum)
  }


  #################################################################
  # Plots
  #################################################################

  colorlist <- c('darkblue','darkred','darkgreen','darkorange','gray50','khaki4','brown3','blue','darkgoldenrod4','cyan4')

  if (is.null(col_bins)){
    col_bins <- colorlist
  }
  if (is.null(pch_bins)){
    pch_bins <- rep(1,cnum)
  }
  if (is.null(col_poly)){
    col_poly <- colorlist
  }
  if (is.null(lty_poly)){
    lty_poly <- rep('solid',cnum)
  }
  if (is.null(col_xline)){
    col_xline <- colorlist
  }
  if (is.null(lty_xline)){
    lty_xline <- rep('dashed',cnum)
  }

  rdmc_plot <- ggplot() + theme_bw() + labs(x='Running variable', y='Outcome')

  if (nobins==FALSE){
    for (c in 1:cnum){
      rdmc_plot <- rdmc_plot + geom_point(aes_string(x=Xmean[,c],y=Ymean[,c]),col=col_bins[c],shape=pch_bins[c],na.rm=TRUE)
    }
  }
  if (!is.null(ci)){
    for (c in 1:cnum){
      rdmc_plot <- rdmc_plot + geom_errorbar(aes_string(x=Xmean[,c],ymin=CI_l[,c],ymax=CI_r[,c]),col=col_bins[c],linetype=1)
    }
  }
  if (nopoly==FALSE){
    for (c in 1:cnum){
      rdmc_plot <- rdmc_plot + geom_line(aes_string(x=X0[,c],y=Yhat0[,c]),col=col_poly[c],linetype=lty_poly[c],na.rm=TRUE) +
        geom_line(aes_string(x=X1[,c],y=Yhat1[,c]),col=col_poly[c],linetype=lty_poly[c],na.rm=TRUE)
    }
  }
  if (noxline==FALSE){
    for (c in 1:cnum){
      rdmc_plot <- rdmc_plot + geom_vline(xintercept=clist[c],col=col_xline[c],linetype=lty_xline[c])
    }
  }

  if (nodraw==FALSE){
    print(rdmc_plot)
  }

  if (count_fail>0){
    warning("rdplot() could not run in one or more cutoffs.")
  }

  #################################################################
  # Return values
  #################################################################

  if (is.null(ci)){
    output <- list(clist = clist,
                  cnum = cnum,
                  X0 = X0,
                  X1 = X1,
                  Yhat0 = Yhat0,
                  Yhat1 = Yhat1,
                  Xmean = Xmean,
                  Ymean = Ymean,
                  rdmc_plot = rdmc_plot,
                  cfail = Cfail)
  }  else {
    output <- list(clist = clist,
                  cnum = cnum,
                  X0 = X0,
                  X1 = X1,
                  Yhat0 = Yhat0,
                  Yhat1 = Yhat1,
                  Xmean = Xmean,
                  Ymean = Ymean,
                  rdmc_plot = rdmc_plot,
                  CI_l = CI_l,
                  CI_r = CI_r,
                  cfail = Cfail)
  }


  return(output)
}
