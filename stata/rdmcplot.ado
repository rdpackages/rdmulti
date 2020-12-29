********************************************************************************
* RDMCPLOT: Regression discontinuity plots with multiple cutoffs
* !version 0.6 2020-12-28
* Authors: Matias Cattaneo, Rocío Titiunik, Gonzalo Vazquez-Bare
********************************************************************************

capture program drop rdmcplot
program define rdmcplot, rclass 
	syntax varlist (min=2 max=2) [if] [in], Cvar(string) [Pvar(string) NBINSvar(string) NBINSRightvar(string) ///
														  COVSvar(string) COVSEVALvar(string) COVSDROPvar(string) BINSELECTvar(string) ///
														  SCALEvar(string) SCALERightvar(string) ///
														  KERNELvar(string) WEIGHTSvar(string) ///
														  Hvar(string) HRightvar(string) ///
														  SUPPORTvar(string) SUPPORTRightvar(string) /// 
														  BINSOPTvar(string) LINEOPTvar(string) XLINEOPTvar(string) /// 
														  ci(real 0) NObins NOpoly NOscatter NOxline NOdraw genvars]
	
	
********************************************************************************
** Setup and error checking
********************************************************************************

	marksample touse, novarlist
	
	tokenize `varlist'
	local yvar `1' 
	local xvar `2' 
	
	capture confirm numeric variable `cvar'
	if _rc!=0 {
		di as error "cutoff variable has to be numeric"
		exit 108
	}
	
	qui sum `xvar' if `touse'
	local xmax = r(max)
	local xmin = r(min)
	qui sum `cvar' if `touse'
	local cmax = r(max)
	local cmin = r(min)
	
	if `cmax'>=`xmax' | `cmin'<=`xmin' {
		di as error "cutoff variable outside range of running variable"
		exit 125
	}

	tempvar treated
	gen double `treated' = `xvar'>=`cvar'

	qui levelsof `cvar' if `touse', local(clist)
	local n_cutoffs: word count `clist'
	
	if "`pvar'"!=""{
		capture confirm numeric variable `pvar'
		if _rc!=0 {
			di as error "p variable has to be numeric"
			exit 108
		}
		qui count if `pvar'!=.
		local n_pvar = r(N)
		if `n_pvar' != `n_cutoffs' {
			di as error "lengths of pvar and cvar have to coincide"
			exit 125
		}
	}

	if "`nbinsvar'"!=""{
		capture confirm numeric variable `nbinsvar'
		if _rc!=0 {
			di as error "nbins variable has to be numeric"
			exit 108
		}
		qui count if `nbinsvar'!=.
		local n_nbinsvar = r(N)
		if `n_nbinsvar' != `n_cutoffs' {
			di as error "lengths of nbins and cvar have to coincide"
			exit 125
		}
		if "`nbinsrightvar'"==""{
			tempvar nbinsrightvar
			qui gen `nbinsrightvar' = `nbinsvar'
		}
		else{
			capture confirm numeric variable `nbinsrightvar'
			if _rc!=0 {
				di as error "nbinsright variable has to be numeric"
				exit 108
			}
			qui count if `nbinsrightvar'!=.
			local n_nbinsrightvar = r(N)
			if `n_nbinsrightvar' != `n_cutoffs' {
				di as error "lengths of nbinsright and cvar have to coincide"
				exit 125
			}
		}
	}	
	
	if "`covsvar'"!=""{
		capture confirm string variable `covsvar'
		if _rc!=0 {
			di as error "covs variable has to be string"
			exit 108
		}
		qui count if `covsvar'!=""
		local n_covsvar = r(N)
		if `n_covsvar' != `n_cutoffs' {
			di as error "lengths of covs and cvar have to coincide"
			exit 125
		}
	}
	
	if "`covsevalvar'"!=""{
		capture confirm string variable `covsevalvar'
		if _rc!=0 {
			di as error "covseval variable has to be string"
			exit 108
		}
		qui count if `covsevalvar'!=""
		local n_covsevalvar = r(N)
		if `n_covsevalvar' != `n_cutoffs' {
			di as error "lengths of covseval and cvar have to coincide"
			exit 125
		}
	}
	
	if "`covsdropvar'"!=""{
		capture confirm string variable `covsdropvar'
		if _rc!=0 {
			di as error "covsdrop variable has to be string"
			exit 108
		}
		qui count if `covsdropvar'!=""
		local n_covsdropvar = r(N)
		if `n_covsdropvar' != `n_cutoffs' {
			di as error "lengths of covsdropvar and cvar have to coincide"
			exit 125
		}
	}
	
	if "`binselectvar'"!=""{
		capture confirm string variable `binselectvar'
		if _rc!=0 {
			di as error "binselect variable has to be string"
			exit 108
		}
		qui count if `binselectvar'!=""
		local n_binselectvar = r(N)
		if `n_binselectvar' != `n_cutoffs' {
			di as error "lengths of binselect and cvar have to coincide"
			exit 125
		}
	}
	
	if "`scalevar'"!=""{
		capture confirm numeric variable `scaleleftvar'
		if _rc!=0 {
			di as error "scale variable has to be numeric"
			exit 108
		}
		qui count if `scalevar'!=.
		local n_scalevar = r(N)
		if `n_scalevar' != `n_cutoffs' {
			di as error "lengths of scale and cvar have to coincide"
			exit 125
		}
		if "`scalerightvar'"==""{
			tempvar scalerightvar
			qui gen `scalerightvar' = `scalevar'
		}
		else{
			capture confirm numeric variable `scalerightvar'
			if _rc!=0 {
				di as error "scaleright variable has to be numeric"
				exit 108
			}
			qui count if `scalerightvar'!=.
			local n_scalerightvar = r(N)
			if `n_scalerightvar' != `n_cutoffs' {
				di as error "lengths of scaleright and cvar have to coincide"
				exit 125
			}
		}
	}	
	
	if "`kernelvar'"!=""{
		capture confirm string variable `kernelvar'
		if _rc!=0 {
			di as error "kernel variable has to be string"
			exit 108
		}
		qui count if `kernelvar'!=""
		local n_kernelvar = r(N)
		if `n_kernelvar' != `n_cutoffs' {
			di as error "lengths of kernelvar and cvar have to coincide"
			exit 125
		}
	}
	
	if "`weightsvar'"!=""{
		capture confirm string variable `weightsvar'
		if _rc!=0 {
			di as error "weights variable has to be string"
			exit 108
		}
		qui count if `weightsvar'!=""
		local weightsvar = r(N)
		if `n_weightsvar' != `n_cutoffs' {
			di as error "lengths of weightsvar and cvar have to coincide"
			exit 125
		}
	}	
	
	if "`hvar'"!=""{
		capture confirm numeric variable `hvar'
		if _rc!=0 {
			di as error "h variable has to be numeric"
			exit 108
		}		
		qui count if `hvar'!=.
		local n_hvar = r(N)
		if `n_hvar' != `n_cutoffs' {
			di as error "lengths of hvar and cvar have to coincide"
			exit 125
		}
		if "`hrightvar'"==""{
			tempvar hrightvar
			qui gen `hrightvar' = `hvar'
		}
		else {
			capture confirm numeric variable `hrightvar'
			if _rc!=0 {
				di as error "hright variable has to be numeric"
				exit 108
			}		
			qui count if `hrightvar'!=.
			local n_hrightvar = r(N)
			if `n_hrightvar' != `n_cutoffs' {
				di as error "lengths of hrightvar and cvar have to coincide"
				exit 125
			}
		}
	}
	
	if "`supportvar'"!=""{
		capture confirm numeric variable `supportvar'
		if _rc!=0 {
			di as error "support variable has to be numeric"
			exit 108
		}
		qui count if `supportvar'!=.
		local n_supportvar = r(N)
		if `n_supportvar' != `n_cutoffs' {
			di as error "lengths of support and cvar have to coincide"
			exit 125
		}
		if "`supportrightvar'"==""{
			tempvar supportrightvar
			qui gen `supportrightvar' = `supportvar'
		}
		else{
			capture confirm numeric variable `supportrightvar'
			if _rc!=0 {
				di as error "supportright variable has to be numeric"
				exit 108
			}
			qui count if `supportrightvar'!=.
			local n_supportrightvar = r(N)
			if `n_supportrightvar' != `n_cutoffs' {
				di as error "lengths of supportright and cvar have to coincide"
				exit 125
			}
		}
	}
	
	if "`binsoptvar'"!=""{
		capture confirm string variable `binsoptvar'
		if _rc!=0 {
			di as error "binsopt variable has to be string"
			exit 108
		}
		qui count if `binsoptvar'!=""
		local n_binsoptvar = r(N)
		if `n_binsoptvar' != `n_cutoffs' {
			di as error "lengths of binsopt and cvar have to coincide"
			exit 125
		}
	}
	
	if "`lineoptvar'"!=""{
		capture confirm string variable `lineoptvar'
		if _rc!=0 {
			di as error "lineopt variable has to be string"
			exit 108
		}
		qui count if `lineoptvar'!=""
		local n_lineoptvar = r(N)
		if `n_lineoptvar' != `n_cutoffs' {
			di as error "lengths of lineopt and cvar have to coincide"
			exit 125
		}
	}
	
	if "`xlineoptvar'"!=""{
		capture confirm string variable `xlineoptvar'
		if _rc!=0 {
			di as error "xlineopt variable has to be string"
			exit 108
		}
		qui count if `xlineoptvar'!=""
		local n_xlineoptvar = r(N)
		if `n_xlineoptvar' != `n_cutoffs' {
			di as error "lengths of xlineopt and cvar have to coincide"
			exit 125
		}
	}
	
	if ("`nobins'"!="" | "`noscatter'"!="") & "`nopoly'"!="" {
		di as error "cannot specify nobins and nopoly simultaneously"
		exit 198
	}
	
	local colorlist "navy maroon green dkorange gray khaki cranberry blue gold cyan"
	local colorlist1 "`colorlist'"
	
	
********************************************************************************
** Generate variables to plot
********************************************************************************	

	local i = 1
	local count = 0
	foreach c of numlist `clist'{
		
		tempvar yhat_`i' ybar_`i' xbar_`i' cileft_`i' ciright_`i'
		
		if "`colorlist1'"==""{
			local colorlist1 "`colorlist'"
			local count = `count' + 1
		}
		
		gettoken color colorlist1 : colorlist1
		
		local f = 1 - 0.2*`count'
		local color "`color'*`f'"
		
		if "`pvar'"!=""{
			local p = `pvar'[`i']
			local p_opt "p(`p')"
		}
		
		if "`nbinsvar'"!=""{
			local nbl = `nbinsvar'[`i']
			local nbr = `nbinsrightvar'[`i']
			local nbins_opt "nbins(`nbl' `nbr')"
		}
		
		if "`covsvar'"!=""{
			local covs = `covsvar'[`i']
			local covs_opt "covs(`covs')"
		}
		
		if "`covsevalvar'"!=""{
			local covseval = `covsevalvar'[`i']
			local covseval_opt "covseval(`covs_eval')"
		}
		
		if "`covsdropvar'"!=""{
			local covsdrop = `covsdropvar'[`i']
			local covsdrop_opt "covsdrop(`covs_drop')"
		}
		
		if "`binselectvar'"!=""{
			local binselect = `binselectvar'[`i']
			local binselect_opt "binselect(`binselect')"
		}
		
		if "`scalevar'"!=""{
			local scalel = `scalevar'[`i']
			local scaler = `scalerightvar'[`i']
			local scale_opt "scale(`scalel' `scaler')"
		}
		
		if "`kernelvar'"!=""{
			local kernel = `kernelvar'[`i']
			local k_opt "kernel(`kernel')"
		}
		
		if "`weightsvar'"!=""{
			local weights = `weightsvar'[`i']
			local weights_opt "weights(`weights')"
		}	
		
		if "`hvar'"!=""{
			local hleft = `hvar'[`i']		
			local hright = `hrightvar'[`i']			
			local h_opt "h(`hleft' `hright')"
			local range_cond "& `c'-`hleft'<=`xvar' & `xvar'<=`c'+`hright'"
		}
		
		if "`supportvar'"!=""{
			local supleft = `supporvar'[`i']		
			local supright = `supportrightvar'[`i']			
			local support_opt "support(`supleft' `supright')"
		}
		
		if "`binsoptvar'"!=""{
			local binsopt = `binsoptvar'[`i']
			local bins_opt "`binsopt'"
		}
		
		if "`lineoptvar'"!=""{
			local lineopt = `lineoptvar'[`i']
			local line_opt "`lineopt'"
		}
		
		if "`xlineoptvar'"!=""{
			local xlineopt = `xlineoptvar'[`i']
			local xline_opt "`xlineopt'"
		}
		
		if `ci'==0{

			qui {
				capture drop rdplot_*
				capture rdplot `yvar' `xvar' if abs(`cvar'-`c')<=c(epsfloat) & `touse' `range_cond', c(`c') `p_opt' `nbins_opt' `covs_opt' `covseval_opt' `covsdrop_opt' `binselect_opt' ///
					`scale_opt' `kernel_opt' `weights_opt' `h_opt' `support_opt' genvars hide
				
				if _rc!=0{
					di as error "rdplot could not run in cutoff `c'. Please check this cutoff manually." 
					exit 2001
				}
				
				gen double `yhat_`i'' = rdplot_hat_y
				gen double `ybar_`i'' = rdplot_mean_y	
				gen double `xbar_`i'' = rdplot_mean_x

			}

			if "`genvars'" != ""{
				qui gen double rdmcplot_hat_y_`i' = `yhat_`i''
				qui gen double rdmcplot_mean_y_`i' = `ybar_`i''
				qui gen double rdmcplot_mean_x_`i' = `xbar_`i''

				label variable rdmcplot_hat_y_`i' "Predicted polynomial for c=`c'"
				label variable rdmcplot_mean_y_`i' "Bin mean of y for c=`c'"
				label variable rdmcplot_mean_x_`i' "Bin mean of x for c=`c'"

			}

			local scat_plots "`scat_plots' (scatter `ybar_`i'' `xbar_`i'', msize(small) mcolor(`color') `bins_opt')"
			local line_plots "`line_plots' (line `yhat_`i'' `xvar' if `cvar'==`c' & `treated'==1, sort lwidth(medthin) lcolor(`color') `line_opt')(line `yhat_`i'' `xvar' if `cvar'==`c' & `treated'==0, sort lwidth(medthin) lcolor(`color') `line_opt')"

			if "`noxline'" == ""{
				if "`xlineoptvar'" == ""{
					local xline_plots "`xline_plots' xline(`c', lcolor(`color') lwidth(medthin) lpattern(shortdash) `xline_opt')"
				}
				else {
					local xline_plots "`xline_plots' xline(`c', `xline_opt')"
				}
			}

		}
		
		else {
			qui {
				capture drop rdplot_*
				rdplot `yvar' `xvar' if `cvar'==`c' & `touse' `range_cond', c(`c') `p_opt' `nbins_opt' `covs_opt' `covseval_opt' `covsdrop_opt' `binselect_opt' ///
					`scale_opt' `kernel_opt' `weights_opt' `h_opt' `support_opt' genvars hide ci(`ci')
				gen double `yhat_`i'' = rdplot_hat_y
				gen double `ybar_`i'' = rdplot_mean_y	
				gen double `xbar_`i'' = rdplot_mean_x
				gen double `cileft_`i'' = rdplot_ci_l
				gen double `ciright_`i'' = rdplot_ci_r 

			}

			if "`genvars'" != ""{
				qui gen double rdmcplot_hat_y_`i' = `yhat_`i''
				qui gen double rdmcplot_mean_y_`i' = `ybar_`i''
				qui gen double rdmcplot_mean_x_`i' = `xbar_`i''
				qui gen double rdmcplot_ci_l_`i' = `cileft_`i''
				qui gen double rdmcplot_ci_r_`i' = `ciright_`i''
				
				label variable rdmcplot_hat_y_`i' "Predicted polynomial for c=`c'"
				label variable rdmcplot_mean_y_`i' "Bin mean of y for c=`c'"
				label variable rdmcplot_mean_x_`i' "Bin mean of x for c=`c'"
				label variable rdmcplot_ci_l_`i' "Left CI for c=`c'"
				label variable rdmcplot_ci_r_`i' "Right CI for c=`c'"
			}

			local scat_plots "`scat_plots' (scatter `ybar_`i'' `xbar_`i'', msize(small) mcolor(`color') `bins_opt')"
			local line_plots "`line_plots' (line `yhat_`i'' `xvar' if `cvar'==`c' & `treated'==1, sort lwidth(medthin) lcolor(`color') `line_opt')(line `yhat_`i'' `xvar' if `cvar'==`c' & `treated'==0, sort lwidth(medthin) lcolor(`color') `line_opt')"
			local ci_plots "`ci_plots' (rcap `cileft_`i'' `ciright_`i'' `xbar_`i'', lcolor(`color'))"
			
			if "`noxline'" == ""{
				if "`xlineoptvar'" == ""{
					local xline_plots "`xline_plots' xline(`c', lcolor(`color') lwidth(medthin) lpattern(shortdash) `xline_opt')"
				}
				else {
					local xline_plots "`xline_plots' xline(`c', `xline_opt')"
				}
			}
		
		}
		

		
		local ++i
	}
	
	
********************************************************************************
** Plot
********************************************************************************

	if "`nodraw'"==""{
		if "`noscatter'"=="" & "`nobins'"=="" & "`nopoly'"==""{
			twoway `scat_plots' `line_plots' `ci_plots', `xline_plots' legend(off)
		}
		else if ("`noscatter'"!="" | "`nobins'"!="") & "`nopoly'"==""{
			twoway `line_plots', `xline_plots' legend(off)
		}
		else if "`noscatter'"=="" & "`nobins'"=="" & "`nopoly'"!=""{
			twoway `scat_plots' `ci_plots', `xline_plots' legend(off)
		}
	}
	
	
********************************************************************************
** Return values
********************************************************************************
	
	drop rdplot_*
	
	ret local clist `clist'
	ret local cvar `cvar'
	ret scalar n_cutoffs = `n_cutoffs'

end
