*******************************************************************************
** RDMULTI: Analysis of Regression Discontinuity Designs 
** 		    with multiple cutoffs or scores
** RDMCPLOT Illustration file
** Date: 22-Apr-2020
** Authors: Matias Cattaneo, Rocío Titiunik, Gonzalo Vazquez-Bare
********************************************************************************
** net install rdmulti, from(https://sites.google.com/site/rdpackages/rdmulti/stata) replace
********************************************************************************

clear all

********************************************************************************
** Load data and generate plot variables (omitting plot)
********************************************************************************

use simdata_multic, clear
sum 
tab c

rdmcplot y x, c(c) ci(95) genvars nodraw

********************************************************************************
** Replicate default plot
********************************************************************************

twoway (scatter rdmcplot_mean_y_1 rdmcplot_mean_x_1, mcolor(navy) msize(small)) ///
	(line rdmcplot_hat_y_1 rdmcplot_mean_x_1 if t==1, sort lcolor(navy)) ///
	(line rdmcplot_hat_y_1 rdmcplot_mean_x_1 if t==0, sort lcolor(navy)) ///
	(scatter rdmcplot_mean_y_2 rdmcplot_mean_x_2, mcolor(maroon) msize(small)) ///
	(line rdmcplot_hat_y_2 rdmcplot_mean_x_2 if t==1, sort lcolor(maroon)) ///
	(line rdmcplot_hat_y_2 rdmcplot_mean_x_2 if t==0, sort lcolor(maroon)), ///
	xline(33, lcolor(navy) lpattern(dash)) ///
	xline(66, lcolor(maroon) lpattern(dash)) ///
	legend(off) 
	
********************************************************************************
** Replicate plot with confidence intervals
********************************************************************************

twoway (scatter rdmcplot_mean_y_1 rdmcplot_mean_x_1, mcolor(navy) msize(small)) ///
	(line rdmcplot_hat_y_1 rdmcplot_mean_x_1 if t==1, sort lcolor(navy)) ///
	(line rdmcplot_hat_y_1 rdmcplot_mean_x_1 if t==0, sort lcolor(navy)) ///
	(scatter rdmcplot_mean_y_2 rdmcplot_mean_x_2, mcolor(maroon) msize(small)) ///
	(line rdmcplot_hat_y_2 rdmcplot_mean_x_2 if t==1, sort lcolor(maroon)) ///
	(line rdmcplot_hat_y_2 rdmcplot_mean_x_2 if t==0, sort lcolor(maroon)) ///
	(rcap rdmcplot_ci_l_1 rdmcplot_ci_r_1 rdmcplot_mean_x_1, sort lcolor(navy)) ///
	(rcap rdmcplot_ci_l_2 rdmcplot_ci_r_2 rdmcplot_mean_x_2, sort lcolor(maroon)), ///
	xline(33, lcolor(navy) lpattern(dash)) ///
	xline(66, lcolor(maroon) lpattern(dash)) ///
	legend(off)
