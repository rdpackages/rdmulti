*******************************************************************************
** RDMULTI: Analysis of Regression Discontinuity Designs 
** 		    with multiple cutoffs or scores
** Illustration file
** Authors: Matias Cattaneo, Rocío Titiunik, Gonzalo Vazquez-Bare
********************************************************************************
** net install rdmulti, from(https://raw.githubusercontent.com/rdpackages/rdmulti/master/stata) replace
********************************************************************************

clear all
set more off
set linesize 90

********************************************************************************
** MULTIPLE NON-CUMULATIVE CUTOFFS: Setup and summary stats
********************************************************************************

use simdata_multic, clear
sum 
tab c

********************************************************************************
** RDMC
********************************************************************************

* Basic syntax

rdmc y x, c(c)

* rdrobust pooled options 

rdmc y x, c(c) pooled_opt(h(20) p(2)) verbose

* Cutoff-specific bandwidths

gen double h = 11 in 1
replace h = 10 in 2
rdmc y x, c(c) h(h)

* Different bandwidth selectors at each cutoff

gen bwselect = "msetwo" in 1
replace bwselect = "certwo" in 2
rdmc y x, c(c) bwselect(bwselect)

* Add plot

rdmc y x, c(c) plot
graph display coefs
graph display weights

* Post estimation testing

rdmc y x, c(c)
matlist e(b)
lincom c1-c2

********************************************************************************
** RDMCPLOT
********************************************************************************

* Basic syntax

rdmcplot y x, c(c)

* Omit bins plot

rdmcplot y x, c(c) nobins

* Plot TE

gen p = 1 in 1/2
rdmcplot y x, c(c) h(h) p(p)

* Plots by hand

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
** MULTIPLE CUMULATIVE CUTOFFS: Setup and summary stats
********************************************************************************

use simdata_cumul, clear
sum 
tab c

********************************************************************************
** RDMS
********************************************************************************

* Basic syntax

rdms y x, c(c)

* Cutoff-specific bandwidths

gen double h = 11 in 1
replace h = 8 in 2
gen kernel = "uniform" in 1
replace kernel = "triangular" in 2
rdms y x, c(c) h(h) kernel(kernel)

* Restricting the range

gen double range_l = 0 in 1
gen double range_r = 65.5 in 1
replace range_l = 33.5 in 2
replace range_r = 100 in 2
rdms y x, c(c) range(range_l range_r)

* Pooled estimate using rdmc

gen double cutoff = c[1]*(x<=49.5) + c[2]*(x>49.5)
rdmc y x, c(cutoff)

* Plot using rdmcplot

gen binsopt = "mcolor(navy)" in 1/2
gen xlineopt = "lcolor(navy) lpattern(dash)" in 1/2
rdmcplot y x, c(cutoff) binsoptvar(binsopt) xlineopt(xlineopt) nopoly

********************************************************************************
** BIVARIATE SCORE: Setup and summary stats
********************************************************************************

use simdata_multis, clear
sum 
list c1 c2 in 1/3

gen xaux = 50 in 1/50
gen yaux = _n in 1/50
twoway (scatter x2 x1 if t==0, msize(small) mfcolor(white) msymbol(X)) ///
	   (scatter x2 x1 if t==1, msize(small) mfcolor(white) msymbol(T)) ///
	   (function y = 50, range(0 50) lcolor(black) lwidth(medthick)) ///
	   (line yaux xaux, lcolor(black) lwidth(medthick)) ///
	   (scatteri 50 25, msize(large) mcolor(black)) ///
   	   (scatteri 50 50, msize(large) mcolor(black)) ///
   	   (scatteri 25 50, msize(large) mcolor(black)), ///
	   text(25 25 "Treated", size(vlarge)) ///
	   text(60 60 "Control", size(vlarge)) ///
	   legend(off)

********************************************************************************
** RDMS
********************************************************************************

* Basic syntax 

rdms y x1 x2 t, c(c1 c2)

* Cutoff specific bandwidths

gen double h = 15 in 1
replace h = 13 in 2
replace h = 17 in 3
rdms y x1 x2 t, c(c1 c2) h(h)

* Pooled effect

gen double aux1 = abs(50 - x1)
gen double aux2 = abs(50 - x2)
egen xnorm = rowmin(aux1 aux2)
replace xnorm = xnorm*(2*t-1)
rdms y x1 x2 t, c(c1 c2) xnorm(xnorm)

