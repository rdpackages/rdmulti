********************************************************************************
* RDMS: analysis of Regression Discontinuity Designs with multiple scores
* !version 0.4 22-Apr-2020
* Authors: Matias Cattaneo, Rocío Titiunik, Gonzalo Vazquez-Bare
********************************************************************************

capture program drop rdms
program define rdms, eclass sortpreserve

	syntax varlist (min=2 max=4) [if] [in], Cvar(string) [range(string) xnorm(string) pooled_opt(string) ///
														  DERIVvar(string) Pvar(string) Qvar(string) ///
														  Hvar(string) HRightvar(string) Bvar(string) BRightvar(string) ///
														  RHOvar(string) COVSvar(string) COVSDROPvar(string) KERNELvar(string) ///
														  WEIGHTSvar(string) BWSELECTvar(string) VCEvar(string) level(real 95) ///
														  SCALEPARvar(string) SCALEREGULvar(string) fuzzy(string) ///
														  MASSPOINTSvar(string) BWCHECKvar(string) BWRESTRICTvar(string) STDVARSvar(string) ///
														  plot graph_opt(string)]

	
********************************************************************************
** Setup and error checking
********************************************************************************
	
	marksample touse, novarlist
	
	tokenize `varlist'
	local yvar `1' 
	local xvar `2' 
	local xvar2 `3'
	local zvar `4'
	
	if "`xvar2'"!="" & "`zvar'"==""{
		di as error "Need to specify zvar when xvar2 is specified"
		exit 102
	}
	
	tokenize `cvar'
	local cvar `1'
	local cvar2 `2'
	
	if "`cvar2'"=="" & "`xvar2'"!=""{
		di as error "Too few variables specified in cvar"
		exit 102
	}
		
	qui count if `cvar'!=.
	local n_cutoffs = r(N)
	
	if "`cvar2'"!=""{
		qui count if `cvar2'!=.
		local n_cutoffs2 = r(N)
		if `n_cutoffs'!=`n_cutoffs2' {
			di as error "cutoffs coordinates incorrectly specified"
			exit 198
		}
	}
	
	if "`range'"!=""{
		tempvar range_c range_t
		tokenize `range'
		qui gen double `range_c' = `1' - `cvar'
		if "`2'"==""{
			local range_t `1'
		} 
		else {
			qui gen double `range_t' = `2' - `cvar'
		}
	}

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
			di as error "lengths of bwselectvar and cvar have to coincide"
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
			di as error "lengths of vcevar and cvar have to coincide"
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
	
	if "`xnorm'"==""{
		mat `b' = J(1,`n_cutoffs',.)
		mat `V' = J(`n_cutoffs',`n_cutoffs',0)
		mat sampsis = J(2,`n_cutoffs',.)
		mat coefs = J(1,`n_cutoffs',.)
		mat CI_rb = J(2,`n_cutoffs',.)
		mat H = J(2,`n_cutoffs',.)
		mat pv_rb = J(1,`n_cutoffs',.)
	} 
	else {
		mat `b' = J(1,`n_cutoffs'+1,.)
		mat `V' = J(`n_cutoffs'+1,`n_cutoffs'+1,0)
		mat sampsis = J(2,`n_cutoffs'+1,.)
		mat coefs = J(1,`n_cutoffs'+1,.)
		mat CI_rb = J(2,`n_cutoffs'+1,.)	
		mat H = J(2,`n_cutoffs'+1,.)
		mat pv_rb = J(1,`n_cutoffs'+1,.)
	}


********************************************************************************	
** Calculate cutoff-specific estimates
********************************************************************************

	if "`xvar2'"==""{
		
		forvalues c = 1/`n_cutoffs'{
		
			local cutoff_`c' = round(`cvar'[`c'],.001)
		
			tempvar xc_`c'
			qui gen double `xc_`c'' = `xvar' - `cvar'[`c']

			if "`range'"==""{
				tempvar range_c range_t
				qui gen double `range_c' = .
				qui gen double `range_t' = .
				qui sum `xc_`c''
				qui replace `range_c' = r(min) in `c'
				qui replace `range_t' = r(max) in `c'
			}
			
			if "`derivvar'"!=""{
				local deriv = `derivvar'[`c']
				local deriv_opt "deriv(`deriv')"
			}
			
			if "`pvar'"!=""{
				local p = `pvar'[`c']
				local p_opt "p(`p')"
			}
			
			if "`qvar'"!=""{
				local q = `qvar'[`c']
				local q_opt "q(`q')"
			}		
			
			if "`hvar'"!=""{
				local hleft = `hvar'[`c']
				local hright = `hrightvar'[`c']
				local h_opt "h(`hleft' `hright')"
			}
			
			if "`bvar'"!=""{
				local bleft = `bvar'[`c']
				local bright = `brightvar'[`c']
				local b_opt "b(`bleft' `bright')"
			}
			
			if "`rhovar'"!=""{
				local rho = `rhovar'[`c']
				local rho_opt "rho(`rho')"
			}

			if "`covsvar'"!=""{
				local covs = `covsvar'[`c']
				local covs_opt "covs(`covs')"
			}
			
			if "`covsdropvar'"!=""{
				local covsdrop = `covsdropvar'[`count']
				local covsdrop_opt "covs_drop(`covsdrop')"
			}
			
			if "`kernelvar'"!=""{
				local kernel = `kernelvar'[`c']
				local k_opt "kernel(`kernel')"
			}
			
			if "`weightsvar'"!=""{
				local weights = `weightsvar'[`c']
				local weights_opt "weights(`weights')"
			}
			
			if "`bwselectvar'"!=""{
				local bwselect = `bwselectvar'[`c']
				local bwselect_opt "bwselect(`bwselect')"
			}
			
			if "`vcevar'"!=""{
				local vce = `vcevar'[`c']
				local vce_opt "vce(`vce')"
			}
			
			if "`scaleparvar'"!=""{
				local scalepar = `scaleparvar'[`c']
				local scalepar_opt "scalepar(`scalepar')"
			}		
			
			if "`scaleregulvar'"!=""{
				local scaleregul = `scaleregulvar'[`c']
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
			
			qui rdrobust `yvar' `xc_`c'' if `range_c'[`c']<=`xc_`c'' & `xc_`c''<=`range_t'[`c'] & `touse', `deriv_opt' `p_opt' `q_opt' `h_opt' `b_opt' `rho_opt' `covs_opt' `covsdrop_opt' ///
																											`k_opt' `weights_opt' `bwselect_opt' `vce_opt' `scalepar_opt' `scaleregul_opt' ///
																											`fuzzy_opt'  level(`level') `masspoints_opt' `bwcheck_opt' `bwrestrict_opt' `stdvars_opt'
			
			local h_`c' = e(h_l)
			local n_h_`c' = e(N_h_l) + e(N_h_r)
			local tau_`c' = e(tau_cl)
			local se_rb_`c' = e(se_tau_rb)
			local pv_rb_`c' = e(pv_rb)
			local ci_l_`c' = e(ci_l_rb)
			local ci_r_`c' = e(ci_r_rb)
			
			local colname "`colname' c`c'"
			
			mat `b'[1,`c'] = e(tau_bc)
			mat `V'[`c',`c'] = e(se_tau_rb)^2
			mat coefs[1,`c'] = e(tau_cl)
			mat CI_rb[1,`c'] = e(ci_l_rb)
			mat CI_rb[2,`c'] = e(ci_r_rb)
			mat sampsis[1,`c'] = e(N_h_l)
			mat sampsis[2,`c'] = e(N_h_r)
			mat H[1,`c'] = e(h_l)
			mat H[2,`c'] = e(h_r)
			mat pv_rb[1,`c'] = e(pv_rb)

		}
	}
	
	else {
		
		forvalues c = 1/`n_cutoffs'{
		
			local cutoff_`c'_1 = round(`cvar'[`c'],.001)
			local cutoff_`c'_2 = round(`cvar2'[`c'],.001)
			local cutoff_`c' = abbrev("(`cutoff_`c'_1',`cutoff_`c'_2')",19)
		
			* Calculate (Euclidean) distance to cutoff
			
			tempvar xc_`c'
			qui gen double `xc_`c'' = sqrt((`xvar'-`cvar'[`c'])^2+(`xvar2'-`cvar2'[`c'])^2)*(2*`zvar'-1)
			
			if "`range'"==""{
				tempvar range_c range_t
				qui gen double `range_c' = .
				qui gen double `range_t' = .
				qui sum `xc_`c''
				qui replace `range_c' = abs(r(min)) in `c'
				qui replace `range_t' = abs(r(max)) in `c'
			}
			
			if "`derivvar'"!=""{
				local deriv = `derivvar'[`c']
				local deriv_opt "deriv(`deriv')"
			}
			
			if "`pvar'"!=""{
				local p = `pvar'[`c']
				local p_opt "p(`p')"
			}
			
			if "`qvar'"!=""{
				local q = `qvar'[`c']
				local q_opt "q(`q')"
			}		
			
			if "`hvar'"!=""{
				local hleft = `hvar'[`c']
				local hright = `hrightvar'[`c']
				local h_opt "h(`hleft' `hright')"
			}
			
			if "`bvar'"!=""{
				local bleft = `bvar'[`c']
				local bright = `brightvar'[`c']
				local b_opt "b(`bleft' `bright')"
			}
			
			if "`rhovar'"!=""{
				local rho = `rhovar'[`c']
				local rho_opt "rho(`rho')"
			}

			if "`covsvar'"!=""{
				local covs = `covsvar'[`c']
				local covs_opt "covs(`covs')"
			}
			
			if "`kernelvar'"!=""{
				local kernel = `kernelvar'[`c']
				local k_opt "kernel(`kernel')"
			}
			
			if "`weightsvar'"!=""{
				local weights = `weightsvar'[`c']
				local weights_opt "weights(`weights')"
			}
			
			if "`bwselectvar'"!=""{
				local bwselect = `bwselectvar'[`c']
				local bwselect_opt "bwselect(`bwselect')"
			}
			
			if "`vcevar'"!=""{
				local vce = `vcevar'[`c']
				local vce_opt "vce(`vce')"
			}
			
			if "`scaleparvar'"!=""{
				local scalepar = `scaleparvar'[`c']
				local scalepar_opt "scalepar(`scalepar')"
			}		
			
			if "`scaleregulvar'"!=""{
				local scaleregul = `scaleregulvar'[`c']
				local scaleregul_opt "scaleregul(`scaleregul')"
			}

			qui rdrobust `yvar' `xc_`c'' if -`range_c'[`c']<=`xc_`c'' & `xc_`c''<=`range_t'[`c'] & `touse', `deriv_opt' `p_opt' `q_opt' `h_opt' `b_opt' `rho_opt' `covs_opt' ///
																											`k_opt' `weights_opt' `bwselect_opt' `vce_opt' `scalepar_opt' `scaleregul_opt' ///
																											`fuzzy_opt'  level(`level')
			
			local h_`c' = e(h_l)
			local n_h_`c' = e(N_h_l) + e(N_h_r)
			local tau_`c' = e(tau_cl)
			local se_rb_`c' = e(se_tau_rb)
			local pv_rb_`c' = e(pv_rb)
			local ci_l_`c' = e(ci_l_rb)
			local ci_r_`c' = e(ci_r_rb)
			
			local colname "`colname' c`c'"
			
			mat `b'[1,`c'] = e(tau_bc)
			mat `V'[`c',`c'] = e(se_tau_rb)^2
			mat coefs[1,`c'] = e(tau_cl)
			mat CI_rb[1,`c'] = e(ci_l_rb)
			mat CI_rb[2,`c'] = e(ci_r_rb)
			mat sampsis[1,`c'] = e(N_h_l)
			mat sampsis[2,`c'] = e(N_h_r)
			mat H[1,`c'] = e(h_l)
			mat H[2,`c'] = e(h_r)
			mat pv_rb[1,`c'] = e(pv_rb)
			
			}
		
	}
	
	
********************************************************************************	
** Calculate pooled estimate
********************************************************************************
	
	if "`xnorm'"!=""{
		
		qui rdrobust `yvar' `xnorm' if `touse', `pooled_opt'
		
		local tau_pooled = e(tau_cl)
		local se_rb_pooled = e(se_tau_rb)
		local pv_rb_pooled = e(pv_rb)
		local ci_l_pooled = e(ci_l_rb)
		local ci_r_pooled = e(ci_r_rb)
		local h_l_pooled = e(h_l)
		local h_r_pooled = e(h_r)
		local N_l_pooled = e(N_l)
		local N_r_pooled = e(N_r)
		local N_eff_l_pooled = e(N_h_l)
		local N_eff_r_pooled = e(N_h_r)
		
		local colname "`colname' pooled"
		
		mat `b'[1,`n_cutoffs'+1] = e(tau_bc)
		mat `V'[`n_cutoffs'+1,`n_cutoffs'+1] = e(se_tau_rb)^2
		mat coefs[1,`n_cutoffs'+1] = e(tau_cl)
		mat CI_rb[1,`n_cutoffs'+1] = e(ci_l_rb)
		mat CI_rb[2,`n_cutoffs'+1] = e(ci_r_rb)
		mat sampsis[1,`n_cutoffs'+1] = e(N_h_l)
		mat sampsis[2,`n_cutoffs'+1] = e(N_h_r)	
		mat H[1,`n_cutoffs'+1] = e(h_l)
		mat H[2,`n_cutoffs'+1] = e(h_r)	
		mat pv_rb[1,`n_cutoffs'+1] = e(pv_rb)
	}
	
********************************************************************************
** Display results
********************************************************************************

	di _newline
	di as text "Cutoff-specific RD estimation with robust bias-corrected inference"
	di as text "{hline 15}{c TT}{hline 64}"
	di as text "{ralign 15:Cutoff}" as text _col(14) "{c |}"	_col(23) "Coef." 					_col(33) "P>|z|"  				_col(43)  "[95% Conf. Int.]"	_col(64) "hl"	_col(71) "hr"	_col(79) "Nh"
	di as text "{hline 15}{c +}{hline 64}"

	forvalues c = 1/`n_cutoffs'{
		di as res "{ralign 15: `cutoff_`c''}"	as text _col(14) "{c |}"	as res	_col(19) %9.3f coefs[1,`c'] 			_col(29)  %9.3f pv_rb[1,`c']			_col(40) %9.2f CI_rb[1,`c'] %9.2f CI_rb[2,`c']						_col(60) %7.2f H[1,`c'] %7.2f H[2,`c'] 						_col(75) %6.0f sampsis[1,`c']+sampsis[2,`c']
	}
		
	if "`xnorm'"!=""{
		di as text "{hline 15}{c +}{hline 64}"
		di as res "{ralign 15:Pooled}"  		as text _col(14) "{c |}"	as res	_col(19) %9.3f coefs[1,`n_cutoffs'+1] 	_col(29)  %9.3f pv_rb[1,`n_cutoffs'+1]	_col(40) %9.2f CI_rb[1,`n_cutoffs'+1] %9.2f CI_rb[2,`n_cutoffs'+1] 	_col(60) %7.2f H[1,`n_cutoffs'+1] %7.2f H[2,`n_cutoffs'+1] 	_col(75) %6.0f sampsis[1,`n_cutoffs'+1]+sampsis[2,`n_cutoffs'+1]
		di as text "{hline 15}{c BT}{hline 64}"
	}
	else {
		di as text "{hline 15}{c BT}{hline 64}"
	}

********************************************************************************
** Plots
********************************************************************************
	
	if "`plot'"!=""{
	
		if "`xnorm'"!=""{	
			capture drop _aux_*
			tempvar aux_count aux_ci_l aux_ci_r aux_pooled aux_cutoffs aux_tag
			
			qui gen `aux_count' = _n in 1/`n_cutoffs'
			
			local xmax_range = `n_cutoffs'+.2
			
			* Plot coefficients
					
			qui gen `aux_ci_l' = `ci_l_pooled' in 1/`n_cutoffs'
			qui gen `aux_ci_r' = `ci_r_pooled' in 1/`n_cutoffs'
			qui gen `aux_pooled' = `tau_pooled' in 1/`n_cutoffs'
			mat Ct = coefs[1,2...]
			mat Ct = Ct'
			mat CIt = CI_rb[1...,2...]
			mat CIt = CIt'
			svmat Ct, names(_aux_coefs)
			svmat CIt, names(_aux_ci)
			twoway (rarea `aux_ci_r' `aux_ci_l' `aux_count', sort color(gs11)) ///
				   (rcap _aux_ci1 _aux_ci2 `aux_count', lcolor(navy)) ///
				   (scatter _aux_coefs1 `aux_count', mcolor(navy)) ///
				   (line `aux_pooled' `aux_count', lcolor(gs6)), ///
				   yline(0, lpattern(shortdash) lcolor(black)) ///
				   xtitle("Cutoff") ytitle("Treatment effect") ///
				   legend(order(3 "Estimate" 2 "95% CI" 4 "Pooled estimate" 1 "95% CI for pooled estimate" )) ///
				   xlabel(1(1)`n_cutoffs') xscale(range(0.8 `xmax_range')) `graph_opt'
			
			drop _aux_*
		}
		else {
			capture drop _aux_*
			tempvar aux_count
			qui gen `aux_count' = _n in 1/`n_cutoffs'
			local xmax_range = `n_cutoffs'+.2
			
			* Plot coefficients
					
			mat Ct = coefs
			mat Ct = Ct'
			mat CIt = CI_rb
			mat CIt = CIt'
			svmat Ct, names(_aux_coefs)
			svmat CIt, names(_aux_ci)
			twoway (rcap _aux_ci1 _aux_ci2 `aux_count', lcolor(navy)) ///
				   (scatter _aux_coefs1 `aux_count', mcolor(navy)), ///
				   yline(0, lpattern(shortdash) lcolor(black)) ///
				   xtitle("Cutoff") ytitle("Treatment effect") ///
				   legend(order(2 "Estimate" 1 "95% CI")) ///
				   xlabel(1(1)`n_cutoffs') xscale(range(0.8 `xmax_range')) `graph_opt'
			
			drop _aux_*
		}
	}
	
	
********************************************************************************
** Return values
********************************************************************************

	matname `b' `colname', columns(.) explicit
	matname `V' `colname', explicit
	
	matname H "left right", rows(.) explicit
	matname sampsis "left right", rows(.) explicit
	
	ereturn post `b' `V'
	
	ereturn matrix sampsis = sampsis
	ereturn matrix H = H
	ereturn matrix CI_rb = CI_rb
	ereturn matrix pv_rb = pv_rb
	ereturn matrix coefs = coefs
	
	ereturn local cmd "rdms"

end
