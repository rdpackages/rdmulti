*******************************************************************************
** RDMULTI: Analysis of Regression Discontinuity Designs 
** 		    with multiple cutoffs or scores
** RDMCPLOT Illustration file
** Authors: Matias Cattaneo, Rocío Titiunik, Gonzalo Vazquez-Bare
********************************************************************************
** net install rdmulti, from(https://raw.githubusercontent.com/rdpackages/rdmulti/master/stata) replace
********************************************************************************

clear all

********************************************************************************
** Load data and generate plot variables (omitting plot)
********************************************************************************

use simdata_multic, clear
sum 
tab c

********************************************************************************
** Replicate default plot
********************************************************************************
	
capture drop rdmcplot_*
rdmcplot y x, c(c) genvars nodraw
twoway (function `r(eq_l_1)', range(`r(range_l_1)') lcolor(navy)) ///
	   (function `r(eq_r_1)', range(`r(range_r_1)') lcolor(navy)) ///
	   (function `r(eq_l_2)', range(`r(range_l_2)') lcolor(maroon)) ///
	   (function `r(eq_r_2)', range(`r(range_r_2)') lcolor(maroon))	///
	   (scatter rdmcplot_mean_y_1 rdmcplot_mean_x_1, mcolor(navy) msize(small)) ///
	   (scatter rdmcplot_mean_y_2 rdmcplot_mean_x_2, mcolor(maroon) msize(small)), ///
	   xline(33, lcolor(navy) lpattern(dash)) ///
	   xline(66, lcolor(maroon) lpattern(dash)) ///
	   legend(off) 
	
********************************************************************************
** Replicate plot with confidence intervals
********************************************************************************

capture drop rdmcplot_*
rdmcplot y x, c(c) ci(95) genvars nodraw
twoway (function `r(eq_l_1)', range(`r(range_l_1)') lcolor(navy)) ///
	   (function `r(eq_r_1)', range(`r(range_r_1)') lcolor(navy)) ///
	   (function `r(eq_l_2)', range(`r(range_l_2)') lcolor(maroon)) ///
	   (function `r(eq_r_2)', range(`r(range_r_2)') lcolor(maroon))	///
	   (scatter rdmcplot_mean_y_1 rdmcplot_mean_x_1, mcolor(navy) msize(small)) ///
	   (scatter rdmcplot_mean_y_2 rdmcplot_mean_x_2, mcolor(maroon) msize(small)) ///
	   (rcap rdmcplot_ci_l_1 rdmcplot_ci_r_1 rdmcplot_mean_x_1, sort lcolor(navy)) ///
	   (rcap rdmcplot_ci_l_2 rdmcplot_ci_r_2 rdmcplot_mean_x_2, sort lcolor(maroon)), ///
	   xline(33, lcolor(navy) lpattern(dash)) ///
	   xline(66, lcolor(maroon) lpattern(dash)) ///
	   legend(off) 
