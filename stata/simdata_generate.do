*******************************************************************************
** RDMULTI: analysis of Regression Discontinuity Designs 
** 		    with multiple cutoffs or scores
** Simulated datasets
** Authors: Matias Cattaneo, Rocío Titiunik, Gonzalo Vazquez-Bare
********************************************************************************

********************************************************************************
** Multiple noncumulative cutoffs
********************************************************************************

clear
set obs 2

gen double c = _n

expand 1000
sort c

gen double x = runiform()
recode c (1=.33)(2=.66)
gen t = x>=c

gen double Y0 = 10 + 18*(x-.5)^3 + .6*(x-.5)
gen double Y1 = 20 + 20*(x-.5)^3 + 6*(x-.5)

gen double y = Y0 + 5*t + rnormal() if c==0.33
replace y = Y1 + 3*t + rnormal() if c==0.66

replace x = x*100
replace y = y*100
replace c = c*100

*twoway (scatter y x if c==33)(scatter y x if c==66)

sort x
drop Y*

*save simdata_multic, replace
*export delimited simdata_multic, replace


********************************************************************************
** Multiple cumulative cutoffs
********************************************************************************

clear
set obs 1000

gen double x = runiform()

gen double Y0 = 10 + 18*(x-.5)^3 + .6*(x-.5)
gen double Y1 = 15.2 + 20*(x-.5)^3 + 6*(x-.5)
gen double Y2 = 20 + 40*(x-.5)^3 - 3*(x-.5)

gen double y = Y0*(x<.33) + Y1*(x>=.33 & x<.66) + Y2*(x>=.66) + rnormal(0,0.5)

replace y = y*100
replace x = x*100

*scatter y x

sort x
gen double c = 33 in 1
replace c = 66 in 2

drop Y*

*save simdata_cumul, replace
*export delimited simdata_cumul, replace


********************************************************************************
** Multiple scores 
********************************************************************************

clear 
set obs 1000

gen double x1 = runiform()
gen double x2 = runiform()

gen t = x1<=.5 & x2<=.5

gen double Y0 = 5 + x1 - .005*x1^2 + .8*x2 + .03*x2^3 + x1*x2 - .008*x1^2*x2^3

gen y = Y0 + 5*t + rnormal()

replace y = y*100
replace x1 = x1*100
replace x2 = x2*100

*twoway (scatter x1 x2 if t==0)(scatter x1 x2 if t==1)

sort x1 x2
gen double c1 = 25 in 1
replace c1 = 50 in 2
replace c1 = 50 in 3
gen double c2 = 50 in 1
replace c2 = 50 in 2
replace c2 = 25 in 3
drop Y0

*save simdata_multis, replace
*export delimited simdata_multis, replace
