###################################################################
# rdmulti: analysis of RD designs with multiple cutoffs or scores
# !version 1.2 22-May-2025
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' rdmulti: analysis of RD Designs with multiple cutoffs or scores
#'
#' The regression discontinuity (RD) design is a popular quasi-experimental design
#' for causal inference and policy evaluation. The \code{'rdmulti'} package provides tools
#' to analyze RD designs with multiple cutoffs or scores: \code{\link{rdmc}()} estimates
#' pooled and cutoff-specific effects in multi-cutoff designs, \code{\link{rdmcplot}()}
#' draws RD plots for multi-cutoff RD designs and \code{\link{rdms}()} estimates effects in
#' cumulative cutoffs or multi-score designs. For more details, and related \code{Stata} and
#' \code{R} packages useful for analysis of RD designs, visit \url{https://rdpackages.github.io/}.
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
#' Calonico, S., M.D. Cattaneo, M. Farrell and R. Titiunik. (2017). \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2017_Stata.pdf}{ \code{rdrobust}: Software for Regression Discontinuity Designs}. \emph{Stata Journal} 17(2): 372-404.
#'
#' Calonico, S., M.D. Cattaneo, and R. Titiunik. (2014). \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_Stata.pdf}{Robust Data-Driven Inference in the Regression-Discontinuity Design}. \emph{Stata Journal} 14(4): 909-946.
#'
#' Calonico, S., M.D. Cattaneo, and R. Titiunik. (2015). \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_R.pdf}{ \code{rdrobust}: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs}. \emph{R Journal} 7(1): 38-51.
#'
#' Cattaneo, M.D., L. Keele, R. Titiunik and G. Vazquez-Bare. (2016). \href{https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2016_JOP.pdf}{Interpreting Regression Discontinuity Designs with Multiple Cutoffs}. \emph{Journal of Politics} 78(4): 1229-1248.
#'
#' Cattaneo, M.D., L. Keele, R. Titiunik and G. Vazquez-Bare. (2020). \href{https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2021_JASA.pdf}{Extrapolating Treatment Effects in Multi-Cutoff Regression Discontinuity Designs}. \emph{Journal of the American Statistical Association} 116(536): 1941, 1952.
#'
#' Cattaneo, M.D., R. Titiunik and G. Vazquez-Bare. (2020). \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf}{Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores}. \emph{Stata Journal} 20(4): 866-891.
#'
#' Keele, L. and R. Titiunik. (2015). \href{https://rdpackages.github.io/references/Keele-Titiunik_2015_PA.pdf}{Geographic Boundaries as Regression Discontinuities}. \emph{Political Analysis} 23(1): 127-155
#'
#' @importFrom graphics abline
#' @importFrom graphics arrows
#' @importFrom graphics legend
#' @importFrom graphics barplot
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom grDevices rgb
#' @importFrom stats poly
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @import ggplot2
#' @import rdrobust
#'
#'
#' @aliases rdmulti_package
"_PACKAGE"
