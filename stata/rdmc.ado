********************************************************************************
* RDMC: analysis of Regression Discontinuity Designs with multiple cutoffs
* !version 0.6 2020-12-28
* Authors: Matias Cattaneo, Rocío Titiunik, Gonzalo Vazquez-Bare
********************************************************************************

capture program drop rdmc
program define rdmc, eclass sortpreserve

	syntax varlist (min=2 max=2) [if] [in], Cvar(string) [pooled_opt(string) ///
														  DERIVvar(string) Pvar(string) Qvar(string) ///
														  Hvar(string) HRightvar(string) Bvar(string) BRightvar(string) ///
														  RHOvar(string) COVSvar(string) COVSDROPvar(string) KERNELvar(string) ////
														  WEIGHTSvar(string) BWSELECTvar(string) VCEvar(string) level(real 95) ///
														  SCALEPARvar(string) SCALEREGULvar(string) fuzzy(string) ///
														  MASSPOINTSvar(string) BWCHECKvar(string) BWRESTRICTvar(string) STDVARSvar(string) ///
														  plot graph_opt(string) verbose]

	
	
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
	
	tempvar rv_norm
	qui gen double `rv_norm' = `xvar' - `cvar'
	
	qui levelsof `cvar' if `touse', local(cutoff_list)
	local n_cutoffs: word count `cutoff_list'
	
	if "`derivvar'"!=""{
		capture confirm numeric variable `bvar'
		if _rc!=0 {
			di as error "deriv variable has to be numeric"
			exit 108
		}
		qui count if `derivvar'!=.
		local n_derivvar = r(N)
		if `n_derivvar' != `n_cutoffs' {
			di as error "lengths of derivvar and cvar have to coincide"
			exit 125
		}
	}
	
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
	
	if "`qvar'"!=""{
		capture confirm numeric variable `qvar'
		if _rc!=0 {
			di as error "q variable has to be numeric"
			exit 108
		}
		qui count if `qvar'!=.
		local n_qvar = r(N)
		if `n_qvar' != `n_cutoffs' {
			di as error "lengths of qvar and cvar have to coincide"
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
	
	if "`bvar'"!=""{
		capture confirm numeric variable `bvar'
		if _rc!=0 {
			di as error "b variable has to be numeric"
			exit 108
		}
		qui count if `bvar'!=.
		local n_bvar = r(N)
		if `n_bvar' != `n_cutoffs' {
			di as error "lengths of bvar and cvar have to coincide"
			exit 125
		}
		if "`brightvar'"==""{
			tempvar brightvar
			qui gen `brightvar' = `bvar'
		}
		else {
			capture confirm numeric variable `brightvar'
			if _rc!=0 {
				di as error "bright variable has to be numeric"
				exit 108
			}		
			qui count if `brightvar'!=.
			local n_brightvar = r(N)
			if `n_brightvar' != `n_cutoffs' {
				di as error "lengths of brightvar and cvar have to coincide"
				exit 125
			}
		}		
	}
	
	if "`rhovar'"!=""{
		capture confirm numeric variable `rhovar'
		if _rc!=0 {
			di as error "rho variable has to be numeric"
			exit 108
		}		
		qui count if `rhovar'!=.
		local n_rhovar = r(N)
		if `n_rhovar' != `n_cutoffs' {
			di as error "lengths of rhovar and cvar have to coincide"
			exit 125
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
			di as error "lengths of covsvar and cvar have to coincide"
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
	
	if "`bwselectvar'"!=""{
		capture confirm string variable `bwselectvar'
		if _rc!=0 {
			di as error "bwselect variable has to be string"
			exit 108
		}
		qui count if `bwselectvar'!=""
		local n_bwselectvar = r(N)
		if `n_bwselectvar' != `n_cutoffs' {
			di as error "lengths of bwselect and cvar have to coincide"
			exit 125
		}
	}
	
	if "`vcevar'"!=""{
		capture confirm string variable `vcevar'
		if _rc!=0 {
			di as error "vce variable has to be string"
			exit 108
		}
		qui count if `vcevar'!=""
		local n_vcevar = r(N)
		if `n_vcevar' != `n_cutoffs' {
			di as error "lengths of vce and cvar have to coincide"
			exit 125
		}
	}
	
	if "`scaleparvar'"!=""{
		capture confirm numeric variable `scaleparvar'
		if _rc!=0 {
			di as error "scalepar variable has to be numeric"
			exit 108
		}
		qui count if `scaleparvar'!=.
		local n_scaleparvar = r(N)
		if `n_scaleparvar' != `n_cutoffs' {
			di as error "lengths of scaleparvar and cvar have to coincide"
			exit 125
		}
	}
	
	if "`scaleregulvar'"!=""{
		capture confirm numeric variable `scaleregulvar'
		if _rc!=0 {
			di as error "scaleregul variable has to be numeric"
			exit 108
		}
		qui count if `scaleregulvar'!=.
		local n_scaleregulvar = r(N)
		if `n_scaleregulvar' != `n_cutoffs' {
			di as error "lengths of scaleregulvar and cvar have to coincide"
			exit 125
		}
	}
	
	if "`masspointsvar'"!=""{
		capture confirm string variable `masspointsvar'
		if _rc!=0 {
			di as error "masspoints variable has to be string"
			exit 108
		}
		qui count if `masspointsvar'!=""
		local n_masspointsvar = r(N)
		if `n_masspointsvar' != `n_cutoffs' {
			di as error "lengths of masspointsvar and cvar have to coincide"
			exit 125
		}
	}
	
	if "`bwcheckvar'"!=""{
		capture confirm string variable `bwcheckvar'
		if _rc!=0 {
			di as error "bwcheck variable has to be string"
			exit 108
		}
		qui count if `bwcheckvar'!=""
		local n_bwcheckvar = r(N)
		if `n_bwcheckvar' != `n_cutoffs' {
			di as error "lengths of bwcheckvar and cvar have to coincide"
			exit 125
		}
	}
	
	if "`bwrestrictvar'"!=""{
		capture confirm string variable `bwrestrictvar'
		if _rc!=0 {
			di as error "bwrestrict variable has to be string"
			exit 108
		}
		qui count if `bwrestrictvar'!=""
		local n_bwrestrictvar = r(N)
		if `n_bwrestrictvar' != `n_cutoffs' {
			di as error "lengths of bwrestrictvar and cvar have to coincide"
			exit 125
		}
	}
	
	if "`stdvarsvar'"!=""{
		capture confirm string variable `stdvarsvar'
		if _rc!=0 {
			di as error "stdvars variable has to be string"
			exit 108
		}
		qui count if `stdvarsvar'!=""
		local n_bwrestrictvar = r(N)
		if `n_stdvarsvar' != `n_cutoffs' {
			di as error "lengths of stdvarsvar and cvar have to coincide"
			exit 125
		}
	}
	
	if "`fuzzy'"!=""{
		local fuzzy_opt "fuzzy(`fuzzy')"
	}
	
	tempname b V
	mat `b' = J(1,`n_cutoffs'+2,.)
	mat `V' = J(`n_cutoffs'+2,`n_cutoffs'+2,0)
	mat sampsis = J(2,`n_cutoffs'+2,.)
	mat weights = J(1,`n_cutoffs',.)
	mat coefs = J(1,`n_cutoffs'+2,.)
	mat CI_rb = J(2,`n_cutoffs'+2,.)
	mat H = J(2,`n_cutoffs'+2,.)
	mat pv_rb = J(1,`n_cutoffs'+2,.)


********************************************************************************
*** Calculate pooled estimate
********************************************************************************
	
	if "`verbose'"==""{
		local quietlylocal "quietly"
	}
	
	`quietlylocal' rdrobust `yvar' `rv_norm' if `touse', `pooled_opt' `fuzzy_opt' level(`level')
	
	local se_rb_pooled = e(se_tau_rb)
	
	mat `b'[1,`n_cutoffs'+2] = e(tau_bc)
	mat `V'[`n_cutoffs'+2,`n_cutoffs'+2] = e(se_tau_rb)^2
	mat coefs[1,`n_cutoffs'+2] = e(tau_cl)
	mat CI_rb[1,`n_cutoffs'+2] = e(ci_l_rb)
	mat CI_rb[2,`n_cutoffs'+2] = e(ci_r_rb)
	mat H[1,`n_cutoffs'+2] = e(h_l)
	mat H[2,`n_cutoffs'+2] = e(h_r)
	mat sampsis[1,`n_cutoffs'+2] = e(N_h_l)
	mat sampsis[2,`n_cutoffs'+2] = e(N_h_r)
	mat pv_rb[1,`n_cutoffs'+2] = e(pv_rb)

	
********************************************************************************	
** Calculate cutoff-specific estimates and weights
********************************************************************************
	
	local count = 1
	foreach cutoff of local cutoff_list{

		* Compute cutoff-specific estimates

		if "`derivvar'"!=""{
			local deriv = `derivvar'[`count']
			local deriv_opt "deriv(`deriv')"
		}
		
		if "`pvar'"!=""{
			local p = `pvar'[`count']
			local p_opt "p(`p')"
		}
		
		if "`qvar'"!=""{
			local q = `qvar'[`count']
			local q_opt "q(`q')"
		}		
		
		if "`hvar'"!=""{
			local hleft = `hvar'[`count']
			local hright = `hrightvar'[`count']
			local h_opt "h(`hleft' `hright')"
		}
		
		if "`bvar'"!=""{
			local bleft = `bvar'[`count']
			local bright = `brightvar'[`count']
			local b_opt "b(`bleft' `bright')"
		}
		
		if "`rhovar'"!=""{
			local rho = `rhovar'[`count']
			local rho_opt "rho(`rho')"
		}

		if "`covsvar'"!=""{
			local covs = `covsvar'[`count']
			local covs_opt "covs(`covs')"
		}
		
		if "`covsdropvar'"!=""{
			local covsdrop = `covsdropvar'[`count']
			local covsdrop_opt "covs_drop(`covsdrop')"
		}
		
		if "`kernelvar'"!=""{
			local kernel = `kernelvar'[`count']
			local k_opt "kernel(`kernel')"
		}
		
		if "`weightsvar'"!=""{
			local weights = `weightsvar'[`count']
			local weights_opt "weights(`weights')"
		}
		
		if "`bwselectvar'"!=""{
			local bwselect = `bwselectvar'[`count']
			local bwselect_opt "bwselect(`bwselect')"
		}
		
		if "`vcevar'"!=""{
			local vce = `vcevar'[`count']
			local vce_opt "vce(`vce')"
		}
		
		if "`scaleparvar'"!=""{
			local scalepar = `scaleparvar'[`count']
			local scalepar_opt "scalepar(`scalepar')"
		}		
		
		if "`scaleregulvar'"!=""{
			local scaleregul = `scaleregulvar'[`count']
			local scaleregul_opt "scaleregul(`scaleregul')"
		}
		
		if "`masspointsvar'"!=""{
			local masspoints = `masspointsvar'[`count']
			local masspoints_opt "masspoints(`masspoints')"
		}
		
		if "`bwcheckvar'"!=""{
			local bwcheck = `bwcheckvar'[`count']
			local bwcheck_opt "bwcheck(`bwcheck')"
		}
		
		if "`bwrestrictvar'"!=""{
			local bwrestrict = `bwrestrict'[`count']
			local bwrestrict_opt "bwrestrict(`bwrestrict')"
		}
		
		if "`stdvarsvar'"!=""{
			local stdvars = `stdvars'[`count']
			local stdvars_opt "stdvars(`stdvars')"
		}

		
		capture rdrobust `yvar' `rv_norm' if abs(`cvar'-`cutoff')<=c(epsfloat) & `touse', `deriv_opt' `p_opt' `q_opt' `h_opt' `b_opt' `rho_opt' `covs_opt' `covsdrop_opt' ///
																		         `k_opt' `weights_opt' `bwselect_opt' `vce_opt' `scalepar_opt' `scaleregul_opt' ///
																		         `fuzzy_opt'  level(`level') `masspoints_opt' `bwcheck_opt' `bwrestrict_opt' `stdvars_opt'
		
		if _rc!=0{
			di as error "rdrobust could not run in cutoff `cutoff'. Please check this cutoff manually." 
			exit 2001
		}
		
		local colname "`colname' c`count'"
		
		mat `b'[1,`count'] = e(tau_bc)
		mat `V'[`count',`count'] = e(se_tau_rb)^2
		mat coefs[1,`count'] = e(tau_cl)
		mat CI_rb[1,`count'] = e(ci_l_rb)
		mat CI_rb[2,`count'] = e(ci_r_rb)
		mat H[1,`count'] = e(h_l)
		mat H[2,`count'] = e(h_r)
		mat sampsis[1,`count'] = e(N_h_l)
		mat sampsis[2,`count'] = e(N_h_r)
		mat pv_rb[1,`count'] = e(pv_rb)
		
		local ++count
	}
	
	* Compute weights
		
	mata: st_matrix("weights",colsum(st_matrix("sampsis")[,1..cols(st_matrix("sampsis"))-2])/sum(colsum(st_matrix("sampsis")[,1..cols(st_matrix("sampsis"))-2])))
	
	* Compute weighted average of cutoff-specific effects
	
	mata: st_numscalar("tweight",st_matrix("coefs")[1,1..cols(st_matrix("coefs"))-2]*st_matrix("weights")')
	mata: st_numscalar("tweight_bc",st_matrix("`b'")[1,1..cols(st_matrix("`b'"))-2]*st_matrix("weights")')
	mata: st_numscalar("Vweight_bc",(st_matrix("weights"):^2)*diagonal(st_matrix("`V'")[1..rows(st_matrix("`V'"))-2,1..cols(st_matrix("`V'"))-2]))	
	mata: st_numscalar("Nweight_l",rowsum(st_matrix("sampsis")[1,1..cols(st_matrix("sampsis"))-2]))
	mata: st_numscalar("Nweight_r",rowsum(st_matrix("sampsis")[2,1..cols(st_matrix("sampsis"))-2]))

	local tstat_weight = tweight_bc/sqrt(Vweight_bc)
	local pval_weight = 2*(1-normal(abs(`tstat_weight')))
	local ci_l_weight = tweight_bc - invnormal(1-(1-`level'/100)/2)*sqrt(Vweight_bc)
	local ci_r_weight = tweight_bc + invnormal(1-(1-`level'/100)/2)*sqrt(Vweight_bc)
	
	mat `b'[1,`n_cutoffs'+1] = tweight_bc
	mat `V'[`n_cutoffs'+1,`n_cutoffs'+1] = Vweight_bc
	mat coefs[1,`n_cutoffs'+1] = tweight
	mat CI_rb[1,`n_cutoffs'+1] = `ci_l_weight'
	mat CI_rb[2,`n_cutoffs'+1] = `ci_r_weight'
	mat sampsis[1,`n_cutoffs'+1] = Nweight_l
	mat sampsis[2,`n_cutoffs'+1] = Nweight_r
	mat pv_rb[1,`n_cutoffs'+1] = `pval_weight'
	
	
********************************************************************************
** Display results
********************************************************************************

	di _newline
	di as text "Cutoff-specific RD estimation with robust bias-corrected inference"
	local count = 1
	di as text "{hline 12}{c TT}{hline 67}"
	di as text "{ralign 12:Cutoff}" as text _col(10) "{c |}"	_col(18) "Coef." 					_col(27) "P>|z|"  				_col(34)  "[`level'% Conf. Int.]"	_col(55) "hl" 	_col(62) "hr"		_col(70) "Nh"				_col(75) "Weight"
	di as text "{hline 12}{c +}{hline 67}"

	foreach c of local cutoff_list{

		di as res %12.3f `c'  		as text _col(10) "{c |}"	as res	_col(13) %9.3f coefs[1,`count'] 		_col(20)  %8.2f pv_rb[1,`count']	_col(33) %8.2f CI_rb[1,`count'] %8.2f CI_rb[2,`count']				_col(51) %7.2f H[1,`count'] 			_col(58) %7.2f H[2,`count']			_col(66) %6.0f sampsis[1,`count']+sampsis[2,`count'] 				_col(72) %9.3f weights[1,`count']			
		local ++count

	}
	
	di as text "{hline 12}{c +}{hline 67}"
	di as res "{ralign 12:Weighted}"  		as text _col(10) "{c |}"	as res	_col(13) %9.3f coefs[1,`n_cutoffs'+1]	_col(20)  %8.2f pv_rb[1,`n_cutoffs'+1]	_col(33) %8.2f CI_rb[1,`n_cutoffs'+1] %8.2f CI_rb[2,`n_cutoffs'+1] 	_col(51) %7.2f .					_col(58) %7.2f .					_col(66) %6.0f sampsis[1,`n_cutoffs'+1]+sampsis[2,`n_cutoffs'+1]	_col(80) "."				
	di as res "{ralign 12:Pooled}"  		as text _col(10) "{c |}"	as res	_col(13) %9.3f coefs[1,`n_cutoffs'+2] 	_col(20)  %8.2f pv_rb[1,`n_cutoffs'+2]	_col(33) %8.2f CI_rb[1,`n_cutoffs'+2] %8.2f CI_rb[2,`n_cutoffs'+2] 	_col(51) %7.2f H[1,`n_cutoffs'+2] 	_col(58) %7.2f H[2,`n_cutoffs'+2]	_col(66) %6.0f sampsis[1,`n_cutoffs'+2]+sampsis[2,`n_cutoffs'+2] 	_col(80) "."			

	di as text "{hline 12}{c BT}{hline 67}"
	
	if "`covsvar'" != ""{
		di as text "Covariate-adjusted estimates."
	}
	
	
********************************************************************************	
** Plots
********************************************************************************	
	
	if "`plot'"!=""{
	
		capture drop _aux_*
		tempvar aux_count aux_ci_l aux_ci_r aux_pooled aux_ci_l_w aux_ci_r_w aux_w aux_cutoffs aux_tag
		
		qui egen `aux_tag' = tag(`cvar')
		qui gen `aux_cutoffs' = `cvar' if `aux_tag'==1
		sort `aux_cutoffs'
		
		qui gen `aux_count' = _n in 1/`n_cutoffs'		
		
		* Plot coefficients
				
		qui gen `aux_ci_l' = CI_rb[1,`n_cutoffs'+2] in 1/`n_cutoffs'
		qui gen `aux_ci_r' = CI_rb[2,`n_cutoffs'+2] in 1/`n_cutoffs'
		qui gen `aux_pooled' = coefs[1,`n_cutoffs'+2] in 1/`n_cutoffs'
		qui gen `aux_ci_l_w' = CI_rb[1,`n_cutoffs'+1] in 1/`n_cutoffs'
		qui gen `aux_ci_r_w' = CI_rb[2,`n_cutoffs'+1] in 1/`n_cutoffs'
		qui gen `aux_w' = coefs[1,`n_cutoffs'+1] in 1/`n_cutoffs'
		mat Ct = coefs[1,2...]
		mat Ct = Ct'
		mat CIt = CI_rb[1...,2...]
		mat CIt = CIt'
		svmat Ct, names(_aux_coefs)
		svmat CIt, names(_aux_ci)
		twoway (rarea `aux_ci_r' `aux_ci_l' `aux_cutoffs', sort color(gs11%40)) ///
			   (rarea `aux_ci_r_w' `aux_ci_l_w' `aux_cutoffs', sort color(maroon%30)) ///
			   (rcap _aux_ci1 _aux_ci2 `aux_cutoffs', lcolor(navy)) ///
			   (scatter _aux_coefs1 `aux_cutoffs', mcolor(navy)) ///
			   (line `aux_pooled' `aux_cutoffs', lcolor(gs6)) ///
   			   (line `aux_w' `aux_cutoffs', lcolor(maroon)), ///
			   name(coefs, replace) ///
			   xtitle("Cutoff") ytitle("Treatment effect") ///
			   legend(order(4 "Cutoff-specific estimates" 3 "`level'% CI" 5 "Pooled with `level'% CI" 6 "Weighted with `level'% CI" )) `graph_opt'
		
		* Plot weights
		
		mat Wt = weights'
		svmat Wt, names(_aux_weights)
		twoway bar _aux_weights1 `aux_count', xtitle("Cutoff") ytitle("Weight") ///
											  ylabel(0(.2)1) ///
											  name(weights, replace) barwidth(.5)
		
		drop _aux_*
	}
	
	
********************************************************************************	
** Return values
********************************************************************************

	matname weights `colname', columns(.) explicit

	local colname "`colname' weighted pooled"
	
	matname `b' `colname', columns(.) explicit
	matname coefs `colname', columns(.) explicit
	matname sampsis `colname', columns(.) explicit
	matname sampsis "left right", rows(.) explicit
	matname H `colname', columns(.) explicit
	matname H "left right", rows(.) explicit
	matname CI_rb `colname', columns(.) explicit
	matname pv_rb `colname', columns(.) explicit

	matname `V' `colname', explicit
	
	ereturn post `b' `V'
	
	ereturn scalar tau_weight = coefs[1,`n_cutoffs'+1]
	ereturn scalar se_weight_rb = sqrt(Vweight_bc)
	ereturn scalar pv_weight_rb = pv_rb[1,`n_cutoffs'+1]
	ereturn scalar ci_weight_l = CI_rb[1,`n_cutoffs'+1]
	ereturn scalar ci_weight_r = CI_rb[2,`n_cutoffs'+1]
	ereturn scalar N_h_l_weight = sampsis[1,`n_cutoffs'+1]
	ereturn scalar N_h_r_weight = sampsis[2,`n_cutoffs'+1]
	ereturn scalar tau_pool = coefs[1,`n_cutoffs'+2]
	ereturn scalar se_pool_rb = `se_rb_pooled'
	ereturn scalar pv_pool_rb = pv_rb[1,`n_cutoffs'+2]
	ereturn scalar ci_weight_l = CI_rb[1,`n_cutoffs'+2]
	ereturn scalar ci_weight_r = CI_rb[2,`n_cutoffs'+2]
	ereturn scalar h_l = H[1,`n_cutoffs'+2]
	ereturn scalar h_r = H[2,`n_cutoffs'+2]
	ereturn scalar N_h_r_pool = sampsis[1,`n_cutoffs'+2]
	ereturn scalar N_h_l_pool = sampsis[2,`n_cutoffs'+2]
	
	ereturn matrix sampsis = sampsis
	ereturn matrix weights = weights
	ereturn matrix H = H
	ereturn matrix CI_rb = CI_rb
	ereturn matrix pv_rb = pv_rb
	ereturn matrix coefs = coefs
	
	ereturn local cmd "rdmc"
	
end
