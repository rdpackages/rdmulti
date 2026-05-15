###################################################################
# rdms: point estimation and robust bias-corrected inference for multi-score designs
# !version 2.0.0 15-May-2026
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Point estimation and robust bias-corrected inference for multi-score designs
#'
#' \code{rdms()} implements point estimation and robust bias-corrected inference
#' for Regression Discontinuity (RD) designs with multiple scores, including
#' cumulative-cutoff designs and designs with two running variables.
#'
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{matias.d.cattaneo@gmail.com}
#'
#' Rocio Titiunik, Princeton University. \email{rocio.titiunik@gmail.com}
#'
#' Gonzalo Vazquez-Bare, UC Santa Barbara. \email{gvazquezbare@gmail.com}
#'
#' @references
#'
#' Cattaneo, M.D., R. Titiunik and G. Vazquez-Bare. (2020). \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf}{Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores}. \emph{Stata Journal}, forthcoming.
#'
#'
#' @param Y outcome variable.
#' @param X running variable.
#' @param C vector of cutoffs.
#' @param X2 if specified, second running variable.
#' @param zvar if X2 is specified, treatment indicator.
#' @param C2 if specified, second vector of cutoffs.
#' @param rangemat matrix of cutoff-specific ranges for the running variable.
#' @param xnorm normalized running variable to estimate pooled effect.
#' @param fuzzy specifies a fuzzy design. See \code{rdrobust()} for details.
#' @param derivvec vector of cutoff-specific order of derivatives. See
#'   \code{rdrobust()} for details.
#' @param pooled_opt options to be passed to \code{rdrobust()} to calculate
#'   pooled estimand.
#' @param pvec vector of cutoff-specific polynomial orders. See
#'   \code{rdrobust()} for details.
#' @param qvec vector of cutoff-specific polynomial orders for bias estimation.
#'   See \code{rdrobust()} for details.
#' @param hmat matrix of cutoff-specific bandwidths. See \code{rdrobust()} for
#'   details.
#' @param bmat matrix of cutoff-specific bandwidths for bias estimation. See
#'   \code{rdrobust()} for details.
#' @param rhovec vector of cutoff-specific values of rho. See \code{rdrobust()}
#'   for details.
#' @param covs_mat matrix of covariates. See \code{rdplot()} for details.
#' @param covs_list list of of covariates to be used in each cutoff.
#' @param covs_dropvec vector indicating whether collinear covariates should be
#'   dropped at each cutoff. See \code{rdrobust()} for details.
#' @param kernelvec vector of cutoff-specific kernels. See \code{rdrobust()} for
#'   details.
#' @param weightsvec vector of length equal to the number of cutoffs indicating
#'   the names of the variables to be used as weights in each cutoff. See \code{rdrobust()} for details.
#' @param bwselectvec vector of cutoff-specific bandwidth selection methods. See
#'   \code{rdrobust()} for details.
#' @param scaleparvec vector of cutoff-specific scale parameters. See
#'   \code{rdrobust()} for details.
#' @param scaleregulvec vector of cutoff-specific scale regularization
#'   parameters. See \code{rdrobust()} for details.
#' @param masspointsvec vector indicating how to handle repeated values at each
#'   cutoff. See \code{rdrobust()} for details.
#' @param bwcheckvec vector indicating the value of bwcheck at each cutoff. See
#'   \code{rdrobust()} for details.
#' @param bwrestrictvec vector indicating whether computed bandwidths are
#'   restricted to the range or runvar at each cutoff. See \code{rdrobust()} for
#'   details.
#' @param stdvarsvec vector indicating whether variables are standardized at
#'   each cutoff. See \code{rdrobust()} for details.
#' @param vcevec vector of cutoff-specific variance-covariance estimation
#'   methods. See \code{rdrobust()} for details.
#' @param nnmatchvec vector of cutoff-specific nearest neighbors for variance
#'   estimation. See \code{rdrobust()} for details.
#' @param cluster cluster ID variable. See \code{rdrobust()} for details.
#' @param ginv.tolvec vector of cutoff-specific tolerances for generalized
#'   inverse calculations. See \code{rdrobust()} for details.
#' @param sharpbwvec vector indicating whether fuzzy RD bandwidth selection
#'   uses the sharp RD model at each cutoff. See \code{rdrobust()} for details.
#' @param level confidence level for confidence intervals. See \code{rdrobust()}
#'   for details.
#' @param plot plots cutoff-specific and pooled estimates.
#' @param conventional reports conventional, instead of robust-bias corrected,
#'   p-values and confidence intervals.
#' @param subset optional subset of observations to use.
#' @param data optional data frame for resolving string variable names.
#'
#'
#' @return \item{B}{vector of bias-corrected coefficients}
#' \item{V}{variance-covariance matrix of the estimators} \item{Coefs}{vector of
#' conventional coefficients} \item{Nh}{vector of sample sizes within bandwidth
#' at each cutoff} \item{CI}{bias corrected confidence intervals}
#' \item{H}{bandwidth used at each cutoff} \item{Pv}{vector of robust p-values}
#'
#'
#' @examples
#' # Toy dataset: cumulative cutoffs
#' X <- runif(1000,0,100)
#' C <- c(33,66)
#' Y <- (1+X)*(X<C[1])+(0.8+0.8*X)*(X>=C[1]&X<C[2])+(1.2+1.2*X)*(X>=C[2]) + rnorm(1000)
#' # rmds: basic syntax
#' tmp <- rdms(Y,X,C)
#'
#'
#' @export

rdms <- function(Y,X,C,X2=NULL,zvar=NULL,C2=NULL,rangemat=NULL,xnorm=NULL,
                 fuzzy=NULL,derivvec=NULL,pooled_opt=NULL,pvec=NULL,qvec=NULL,
                 hmat=NULL,bmat=NULL,rhovec=NULL,covs_mat=NULL,covs_list=NULL,
                 covs_dropvec=NULL,kernelvec=NULL,weightsvec=NULL,
                 bwselectvec=NULL,scaleparvec=NULL,scaleregulvec=NULL,
                 masspointsvec=NULL,bwcheckvec=NULL,bwrestrictvec=NULL,
                 stdvarsvec=NULL,vcevec=NULL,nnmatchvec=NULL,cluster=NULL,
                 ginv.tolvec=NULL,sharpbwvec=NULL,level=95,plot=FALSE,
                 conventional=FALSE,subset=NULL,data=NULL){

  #################################################################
  # Setup and error checking
  #################################################################

  if (!is.null(data)){
    if (is.character(Y) && length(Y)==1) Y <- data[[Y]]
    if (is.character(X) && length(X)==1) X <- data[[X]]
    if (is.character(C) && length(C)==1) C <- unique(na.omit(data[[C]]))
    if (!is.null(X2) && is.character(X2) && length(X2)==1) X2 <- data[[X2]]
    if (!is.null(zvar) && is.character(zvar) && length(zvar)==1) zvar <- data[[zvar]]
    if (!is.null(C2) && is.character(C2) && length(C2)==1) C2 <- unique(na.omit(data[[C2]]))
    if (!is.null(xnorm) && is.character(xnorm) && length(xnorm)==1) xnorm <- data[[xnorm]]
    if (!is.null(fuzzy) && is.character(fuzzy) && length(fuzzy)==1) fuzzy <- data[[fuzzy]]
    if (!is.null(cluster) && is.character(cluster) && length(cluster)==1) cluster <- data[[cluster]]
    if (!is.null(covs_mat) && is.character(covs_mat)) covs_mat <- data[,covs_mat,drop=FALSE]
  }

  n.orig <- length(Y)
  if (!is.null(subset)){
    Y <- Y[subset]
    X <- X[subset]
    if (!is.null(X2) && length(X2)==n.orig) X2 <- X2[subset]
    if (!is.null(zvar) && length(zvar)==n.orig) zvar <- zvar[subset]
    if (!is.null(xnorm) && length(xnorm)==n.orig) xnorm <- xnorm[subset]
    if (!is.null(fuzzy) && length(fuzzy)==n.orig) fuzzy <- fuzzy[subset]
    if (!is.null(cluster) && length(cluster)==n.orig) cluster <- cluster[subset]
    if (!is.null(covs_mat) && nrow(as.matrix(covs_mat))==n.orig) covs_mat <- as.matrix(covs_mat)[subset,,drop=FALSE]
  }

  if (!is.null(X2) & is.null(zvar)) stop('Need to specify zvar when X2 is specified')
  if (!is.null(X2) & is.null(C2)) stop('Need to specify C2 if X2 is specified')

  cnum <- length(C)

  if (!is.null(C2)){
    if (cnum!=length(C2)) stop('cutoff coordinates incorrectly specified')
  }

  if (!is.null(rangemat)){
    if (is.null(dim(rangemat))){
      rangemat <- matrix(rangemat,nrow=cnum,ncol=2)
    }
  }  else {
    rangemat <- cbind(rep(-Inf,cnum),rep(Inf,cnum))
  }

  if (!is.null(hmat)){
    if (is.null(dim(hmat))){
      hmat <- matrix(hmat,nrow=cnum,ncol=2)
    }
  }
  if (!is.null(bmat)){
    if (is.null(dim(bmat))){
      bmat <- matrix(bmat,nrow=cnum,ncol=2)
    }
  }

  if (is.null(covs_dropvec)) covs_dropvec <- rep(TRUE,cnum)
  if (is.null(kernelvec)) kernelvec <- rep('tri',cnum)
  if (is.null(bwselectvec)) bwselectvec <- rep('mserd',cnum)
  if (is.null(vcevec)) vcevec <- rep('nn',cnum)
  if (is.null(nnmatchvec)) nnmatchvec <- rep(3,cnum)
  if (is.null(scaleparvec)) scaleparvec <- rep(1,cnum)
  if (is.null(scaleregulvec)) scaleregulvec <- rep(1,cnum)
  if (is.null(masspointsvec)) masspointsvec <- rep('adjust',cnum)
  if (is.null(bwrestrictvec)) bwrestrictvec <- rep(TRUE,cnum)
  if (is.null(stdvarsvec)) stdvarsvec <- rep(FALSE,cnum)
  if (is.null(ginv.tolvec)) ginv.tolvec <- rep(1e-20,cnum)
  if (is.null(sharpbwvec)) sharpbwvec <- rep(FALSE,cnum)
  if (!is.null(covs_mat)){
    covs_mat <- as.matrix(covs_mat)
    if (!is.null(covs_list)){
      if (length(covs_list)!=cnum) stop('Elements in covs_list should equal number of cutoffs')
    }
  }

  B <- matrix(NA,nrow=1,ncol=cnum+1)
  V <- matrix(NA,nrow=1,ncol=cnum+1)
  Coefs <- matrix(NA,nrow=1,ncol=cnum+1)
  V_cl <- matrix(NA,nrow=1,ncol=cnum+1)
  Nh <- matrix(NA,nrow=2,ncol=cnum+1)
  CI <- matrix(NA,nrow=2,ncol=cnum+1)
  CI_cl <- matrix(NA,nrow=2,ncol=cnum+1)
  Pv <- matrix(NA,nrow=1,ncol=cnum+1)
  Pv_cl <- matrix(NA,nrow=1,ncol=cnum+1)
  H <- matrix(NA,nrow=2,ncol=cnum+1)
  Bbw <- matrix(NA,nrow=2,ncol=cnum+1)

  c.disp <- NULL

  #################################################################
  # Calculate cutoff-specific estimates
  #################################################################

  if (is.null(X2)){

    for (c in 1:cnum){

      xc.full <- X - C[c]
      Rc <- rangemat - C

      use <- xc.full>=Rc[c,1] & xc.full<=Rc[c,2]
      yc <- Y[use]
      xc <- xc.full[use]

      if (!is.null(cluster)){
        cc <- cluster[use]
      } else {
        cc <- NULL
      }

      if (!is.null(weightsvec)){
        weightaux <- weightsvec[c]
        weightsc <- eval(parse(text=weightaux), parent.frame())[use]
      } else{
        weightsc <- NULL
      }

      if (!is.null(covs_mat)){
        covs_mat_c <- covs_mat[use,,drop=FALSE]
        if (!is.null(covs_list)){
          covs_aux <- covs_mat_c[,covs_list[[c]]]
        } else{
          covs_aux <- covs_mat_c
        }
      } else {
        covs_aux <- NULL
      }
      if (!is.null(fuzzy)){
        fuzzyc <- fuzzy[use]
      } else {
        fuzzyc <- NULL
      }

      rdr.tmp <- rdrobust::rdrobust(yc,xc,
                                   fuzzy=fuzzyc,
                                   deriv=derivvec[c],
                                   p=pvec[c],
                                   q=qvec[c],
                                   h=hmat[c,],
                                   b=bmat[c,],
                                   rho=rhovec[c],
                                   covs=covs_aux,
                                   covs_drop=covs_dropvec[c],
                                   ginv.tol=ginv.tolvec[c],
                                   kernel=kernelvec[c],
                                   weights=weightsc,
                                   bwselect=bwselectvec[c],
                                   scalepar=scaleparvec[c],
                                   scaleregul=scaleregulvec[c],
                                   masspoints=masspointsvec[c],
                                   bwcheck=bwcheckvec[c],
                                   bwrestrict=bwrestrictvec[c],
                                   stdvars=stdvarsvec[c],
                                   vce=vcevec[c],
                                   nnmatch=nnmatchvec[c],
                                   cluster=cc,
                                   sharpbw=sharpbwvec[c],
                                   level=level)

      B[1,c] <- rdr.tmp$Estimate[2]
      V[1,c] <- rdr.tmp$se[3]^2
      Coefs[1,c] <- rdr.tmp$Estimate[1]
      V_cl[1,c] <- rdr.tmp$se[1]^2
      CI[,c] <- rdr.tmp$ci[3,]
      CI_cl[,c] <- rdr.tmp$ci[1,]
      H[,c] <- rdr.tmp$bws[1,]
      Bbw[,c] <- rdr.tmp$bws[2,]
      Nh[,c] <- rdr.tmp$N_h
      Pv[1,c] <- rdr.tmp$pv[3]
      Pv_cl[1,c] <- rdr.tmp$pv[1]

      c.disp <- c(c.disp,round(C[c],2))

    }

  } else {

    for (c in 1:cnum){

      xc.full <- sqrt((X-C[c])^2+(X2-C2[c])^2)*(2*zvar-1)

      use <- xc.full>=rangemat[c,1] & xc.full<=rangemat[c,2]
      yc <- Y[use]
      xc <- xc.full[use]

      if (!is.null(cluster)){
        cc <- cluster[use]
      } else {
        cc <- NULL
      }

      if (!is.null(weightsvec)){
        weightaux <- weightsvec[c]
        weightsc <- eval(parse(text=weightaux), parent.frame())[use]
      } else{
        weightsc = NULL
      }

      if (!is.null(covs_mat)){
        covs_mat_c <- covs_mat[use,,drop=FALSE]
        if (!is.null(covs_list)){
          covs_aux <- covs_mat_c[,covs_list[[c]]]
        } else{
          covs_aux <- covs_mat_c
        }
      } else {
        covs_aux <- NULL
      }
      if (!is.null(fuzzy)){
        fuzzyc <- fuzzy[use]
      } else {
        fuzzyc <- NULL
      }

      rdr.tmp <- rdrobust::rdrobust(yc,xc,
                                   fuzzy=fuzzyc,
                                   deriv=derivvec[c],
                                   p=pvec[c],
                                   q=qvec[c],
                                   h=hmat[c,],
                                   b=bmat[c,],
                                   rho=rhovec[c],
                                   covs=covs_aux,
                                   covs_drop=covs_dropvec[c],
                                   ginv.tol=ginv.tolvec[c],
                                   kernel=kernelvec[c],
                                   weights=weightsc,
                                   bwselect=bwselectvec[c],
                                   scalepar=scaleparvec[c],
                                   scaleregul=scaleregulvec[c],
                                   masspoints=masspointsvec[c],
                                   bwcheck=bwcheckvec[c],
                                   bwrestrict=bwrestrictvec[c],
                                   stdvars=stdvarsvec[c],
                                   vce=vcevec[c],
                                   nnmatch=nnmatchvec[c],
                                   cluster=cc,
                                   sharpbw=sharpbwvec[c],
                                   level=level)

      B[1,c] <- rdr.tmp$Estimate[2]
      V[1,c] <- rdr.tmp$se[3]^2
      Coefs[1,c] <- rdr.tmp$Estimate[1]
      V_cl[1,c] <- rdr.tmp$se[1]^2
      CI[,c] <- rdr.tmp$ci[3,]
      CI_cl[,c] <- rdr.tmp$ci[1,]
      H[,c] <- rdr.tmp$bws[1,]
      Bbw[,c] <- rdr.tmp$bws[2,]
      Nh[,c] <- rdr.tmp$N_h
      Pv[1,c] <- rdr.tmp$pv[3]
      Pv_cl[1,c] <- rdr.tmp$pv[1]

      c.disp <- c(c.disp,paste0('(',round(C[c],2),',',round(C2[c],2),')'))

    }

  }


  #################################################################
  # Calculate pooled estimates
  #################################################################

  if (!is.null(xnorm)){

    aux1 <- paste0('rdrobust::rdrobust(Y,xnorm,fuzzy=fuzzy,',pooled_opt,')')

    rdr <- eval(parse(text=aux1))

    B[1,cnum+1] <- rdr$Estimate[2]
    V[1,cnum+1] <- rdr$se[3]^2
    Coefs[1,cnum+1] <- rdr$Estimate[1]
    V_cl[1,cnum+1] <- rdr$se[1]^2
    CI[,cnum+1] <- rdr$ci[3,]
    CI_cl[,cnum+1] <- rdr$ci[1,]
    H[,cnum+1] <- rdr$bws[1,]
    Bbw[,cnum+1] <- rdr$bws[2,]
    Nh[,cnum+1] <- rdr$N_h
    Pv[1,cnum+1] <- rdr$pv[3]
    Pv_cl[1,cnum+1] <- rdr$pv[1]

  }


  #################################################################
  # Display results
  #################################################################

  cat('\n')
  cat(paste0(rep('=',80),collapse='')); cat('\n')
  cat(format('Cutoff',             width=17))
  cat(format('Coef.',              width=9))
  cat(format('P-value',            width=17))
  cat(format(paste0(level,'% CI'), width=16))
  cat(format('hl',                 width=9))
  cat(format('hr',                 width=10))
  cat(format('Nh',                 width=10)); cat('\n')
  cat(paste0(rep('=',80),collapse='')); cat('\n')

  if(conventional==FALSE){
    for (k in 1:cnum){
      if (is.null(C2)){
        cat(format(sprintf('%4.2f',c.disp[k]),       width=17))
      } else {
        cat(paste0('(',format(sprintf('%4.2f',C[k])),',',format(sprintf('%4.2f',C2[k])),format(')',width=5)))
      }
      cat(format(sprintf('%7.3f',Coefs[k]),        width=10))
      cat(format(sprintf('%1.3f',Pv[k]),           width=10))
      cat(format(sprintf('%4.3f',CI[1,k]),         width=10))
      cat(format(sprintf('%4.3f',CI[2,k]),         width=10))
      cat(format(sprintf('%4.3f',H[1,k]),          width=9))
      cat(format(sprintf('%4.3f',H[2,k]),          width=10))
      cat(format(sprintf('%4.0f',Nh[1,k]+Nh[2,k]), width=10))
      cat('\n')
    }

    if (!is.null(xnorm)){
      cat(paste0(rep('-',80),collapse='')); cat('\n')
      cat(format('Pooled',                                   width=17))
      cat(format(sprintf('%7.3f',Coefs[cnum+1]),             width=10))
      cat(format(sprintf('%1.3f',Pv[cnum+1]),                width=10))
      cat(format(sprintf('%4.3f',CI[1,cnum+1]),              width=10))
      cat(format(sprintf('%4.3f',CI[2,cnum+1]),              width=10))
      cat(format(sprintf('%4.3f',H[1,cnum+1]),               width=9))
      cat(format(sprintf('%4.3f',H[2,cnum+1]),               width=10))
      cat(format(sprintf('%4.0f',Nh[1,cnum+1]+Nh[2,cnum+1]), width=10))

      cat('\n')
    }
  } else{
    for (k in 1:cnum){
      if (is.null(C2)){
        cat(format(sprintf('%4.2f',c.disp[k]),       width=17))
      } else {
        cat(paste0('(',format(sprintf('%4.2f',C[k])),',',format(sprintf('%4.2f',C2[k])),format(')',width=5)))
      }
      cat(format(sprintf('%7.3f',Coefs[k]),        width=10))
      cat(format(sprintf('%1.3f',Pv_cl[k]),        width=10))
      cat(format(sprintf('%4.3f',CI_cl[1,k]),      width=10))
      cat(format(sprintf('%4.3f',CI_cl[2,k]),      width=10))
      cat(format(sprintf('%4.3f',H[1,k]),          width=9))
      cat(format(sprintf('%4.3f',H[2,k]),          width=10))
      cat(format(sprintf('%4.0f',Nh[1,k]+Nh[2,k]), width=10))
      cat('\n')
    }

    if (!is.null(xnorm)){
      cat(paste0(rep('-',80),collapse='')); cat('\n')
      cat(format('Pooled',                                   width=17))
      cat(format(sprintf('%7.3f',Coefs[cnum+1]),             width=10))
      cat(format(sprintf('%1.3f',Pv_cl[cnum+1]),             width=10))
      cat(format(sprintf('%4.3f',CI_cl[1,cnum+1]),           width=10))
      cat(format(sprintf('%4.3f',CI_cl[2,cnum+1]),           width=10))
      cat(format(sprintf('%4.3f',H[1,cnum+1]),               width=9))
      cat(format(sprintf('%4.3f',H[2,cnum+1]),               width=10))
      cat(format(sprintf('%4.0f',Nh[1,cnum+1]+Nh[2,cnum+1]), width=10))

      cat('\n')
    }
  }

  cat(paste0(rep('=',80),collapse='')); cat('\n')


  #################################################################
  # Return values
  #################################################################

  colnames(B) <- c(1:cnum,"pooled")
  colnames(V) <- c(1:cnum,"pooled")
  colnames(Coefs) <- c(1:cnum,"pooled")
  colnames(CI) <- c(1:cnum,"pooled")
  colnames(Nh) <- c(1:cnum,"pooled")
  colnames(H) <- c(1:cnum,"pooled")
  colnames(Pv) <- c(1:cnum,"pooled")

  rownames(Nh) <- c("left","right")
  rownames(H) <- c("left","right")

  output <- list(B = B,
                V = V,
                Coefs = Coefs,
                V_cl = V_cl,
                Nh = Nh,
                CI = CI,
                CI_cl = CI_cl,
                H = H,
                Bbw = Bbw,
                Pv = Pv,
                Pv_cl = Pv_cl)

  return(output)

}
