{smcl}
{* *! version 0.6 2021-01-04}{...}
{viewerjumpto "Syntax" "rdms##syntax"}{...}
{viewerjumpto "Description" "rdms##description"}{...}
{viewerjumpto "Options" "rdms##options"}{...}
{viewerjumpto "Examples" "rdms##examples"}{...}
{viewerjumpto "Saved results" "rdms##saved_results"}{...}
{viewerjumpto "References" "rdmc##references"}{...}
{viewerjumpto "Authors" "rdmc##authors"}{...}

{title:Title}

{p 4 8}{cmd:rdms} {hline 2} Analysis of Regression Discontinuity Designs with Multiple Scores.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdms} {it:depvar} {it:runvar1} [{it:runvar2 treatvar}] {ifin}{cmd:,}
{cmd:{opt c:var}(}{it:cvar1 [cvar2]}{cmd:)}
[
{cmd:range(}{it:range1 [range2]}{cmd:)} 
{cmd:xnorm(}{it:string}{cmd:)} 
{cmd:fuzzy(}{it:string}{cmd:)} 
{cmd:{opt deriv:var}(}{it:string}{cmd:)} 
{cmd:pooled_opt(}{it:string}{cmd:)} 
{cmd:{opt p:var}(}{it:string}{cmd:)} 
{cmd:{opt q:var}(}{it:string}{cmd:)} 
{cmd:{opt h:var}(}{it:string}{cmd:)} 
{cmd:{opt hr:ightvar}(}{it:string}{cmd:)} 
{cmd:{opt b:var}(}{it:string}{cmd:)} 
{cmd:{opt br:ightvar}(}{it:string}{cmd:)} 
{cmd:{opt rho:var}(}{it:string}{cmd:)} 
{cmd:{opt covs:var}(}{it:string}{cmd:)} 
{cmd:{opt covsdrop:var}(}{it:string}{cmd:)} 
{cmd:{opt kernel:var}(}{it:string}{cmd:)} 
{cmd:{opt weights:var}(}{it:string}{cmd:)} 
{cmd:{opt bwselect:var}(}{it:string}{cmd:)} 
{cmd:{opt scalepar:var}(}{it:string}{cmd:)} 
{cmd:{opt scaleregul:var}(}{it:string}{cmd:)} 
{cmd:{opt masspoints:var}(}{it:string}{cmd:)} 
{cmd:{opt bwcheck:var}(}{it:string}{cmd:)} 
{cmd:{opt bwrestrict:var}(}{it:string}{cmd:)} 
{cmd:{opt stdvars:var}(}{it:string}{cmd:)} 
{cmd:{opt vce:var}(}{it:string}{cmd:)} 
{cmd:level({it:#})}
{cmd:plot} 
{cmd:graph_opt(}{it:string}{cmd:)} 
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rdms} provides tools to analyze regression discontinuity (RD) designs with multiple scores.
For methodological background see
{browse "https://rdpackages.github.io/references/Keele-Titiunik_2015_PA.pdf":Keele and Titiunik (2015)},
{browse "https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2016_JOP.pdf":Cattaneo, Keele, Titiunik and Vazquez-Bare (2016)}, and
{browse "https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2021_JASA.pdf":Cattaneo, Keele, Titiunik and Vazquez-Bare (2021)}.
It also computes alternative estimation and inference procedures available in the literature.{p_end}

{p 8 8}If only {it:runvar1} is specified, {cmd:rdms} analyzes an RD design with
cumulative cutoffs in which a unit gets different dosages of a treatment depending on the value of {it:runvar1}. If {it:runvar1}, {it:runvar2} and {it:treatvar}
are specified, {cmd:rdms} analyzes an RD design with two running variables in which units with {it:treatvar} equal to one are treated.{p_end}

{p 8 8}Companion commands are: {help rdmc:rdmc} for multi-cutoff RD estimation and inference, and {help rdmcplot:rdmcplot} for multi-cutoff RD plots.{p_end}

{p 8 8}A detailed introduction to this command is given in
{browse "https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_rdmulti.pdf": Cattaneo, Titiunik and Vazquez-Bare (2020)}.{p_end}

{p 8 8}Companion {browse "www.r-project.org":R} functions are also available {browse "https://rdpackages.github.io/rdmulti":here}.{p_end}

{p 8 8}This command employs the Stata (and R) package {help rdrobust:rdrobust} for underlying calculations. See
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_Stata.pdf":Calonico, Cattaneo and Titiunik (2014)},
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_R.pdf":Calonico, Cattaneo and Titiunik (2015)},
and
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2017_Stata.pdf":Calonico, Cattaneo, Farrell and Titiunik (2017)}
for more details.{p_end}

{p 4 8}Related Stata and R packages useful for inference in RD designs are described in the following website:{p_end}

{p 8 8}{browse "https://rdpackages.github.io/":https://rdpackages.github.io/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimand}

{p 4 8}{cmd:{opt c:var}(}{it:string}{cmd:)} specifies the numeric variable {it:cvar1} containing the RD cutoff for {it:indepvar} in a cumulative cutoffs setting,
or the two scores {it:cvar1} and {it:cvar2} in a two-score setting.{p_end}

{p 4 8}{cmd:range(}{it:range1 [range2]}{cmd:)} specifies the range of the running variable to be used for estimation around each cutoff.  Specifying only one variable implies using the same range at each side of the cutoff.{p_end}

{p 4 8}{cmd:xnorm(}{it:string}{cmd:)} specifies the normalized running variable to estimate pooled effect.{p_end}

{p 4 8}{cmd:fuzzy(}{it:string}{cmd:)} indicates a fuzzy design.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt deriv:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the order of the derivative for {cmd:rdrobust} to calculate cutoff-specific estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{dlgtab:rdrobust Options}

{p 4 8}{cmd:pooled_opt(}{it:string}{cmd:)} specifies the options to be passed to {cmd:rdrobust} to calculate pooled estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{dlgtab:Local Polynomial Regression}

{p 4 8}{cmd:{opt p:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the order of the polynomials for {cmd:rdrobust} to calculate cutoff-specific estimates. 
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt q:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the order of the polynomials for bias estimation for {cmd:rdrobust} to calculate cutoff-specific estimates. 
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt h:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the bandwidths for {cmd:rdrobust} to calculate cutoff-specific estimates. 
When {cmd:hrightvar} is specified, {cmd:hvar} indicates the bandwidth to the left of the cutoff.
When {cmd:hrightvar} is not specified, the same bandwidths are used at each side. See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt hr:ightvar}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the bandwidths to the right of the cutoff for {cmd:rdrobust} to calculate cutoff-specific estimates. 
When {cmd:hrightvar} is not specified, the bandwidths in {cmd:hvar} are used at each side.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt b:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the bandwidths for bias estimation for {cmd:rdrobust} to calculate cutoff-specific estimates. 
When {cmd:brightvar} is specified, {cmd:bvar} indicates the bandwidth to the left of the cutoff.
When {cmd:brightvar} is not specified, the same bandwidths are used at each side. See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt br:ightvar}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the bandwidths for bias estimation to the right of the cutoff for {cmd:rdrobust} to calculate cutoff-specific estimates. 
When {cmd:brightvar} is not specified, the bandwidths in {cmd:bvar} are used at each side.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt rho:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the value of rho for {cmd:rdrobust} to calculate cutoff-specific estimates. 
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt covs:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the covariates for {cmd:rdrobust} to calculate cutoff-specific estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt covsdrop:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies whether collinear covariates should be dropped.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt kernel:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the kernels for {cmd:rdrobust} to calculate cutoff-specific estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt weights:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the weights for {cmd:rdrobust} to calculate cutoff-specific estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{dlgtab:Bandwidth Selection}

{p 4 8}{cmd:{opt bwselect:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the bandwidth selection method for {cmd:rdrobust} to calculate cutoff-specific estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt scalepar:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the value of scalepar for {cmd:rdrobust} to calculate cutoff-specific estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt scaleregul:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the value of scaleregul for {cmd:rdrobust} to calculate cutoff-specific estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt masspoints:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies how to handle repeated values in the running variable.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt bwcheck:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the value of {cmd:bwcheck}.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt bwrestrict:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies whether computed bandwidths are restricted to the range of {it:runvar}.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:{opt stdvars:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies whether {it:depvar} and {it:runvar} are standardized.
See {help rdrobust:rdrobust} for details.{p_end}

{dlgtab:Variance-Covariance Estimation}

{p 4 8}{cmd:{opt vce:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the variance-covariance matrix estimation method for {cmd:rdrobust} to calculate cutoff-specific estimates.
See {help rdrobust:rdrobust} for details.{p_end}

{p 4 8}{cmd:level(}{it:#}{cmd:)} specifies the confidence level for confidence intervals.
See {help rdrobust:rdrobust} for details.{p_end}

{dlgtab:Plot}

{p 4 8}{cmd:plot} plots the pooled and cutoff-specific estimates.{p_end}

{p 4 8}{cmd:graph_opt(}{it:string}{cmd:)} options to be passed to the graph when {cmd:plot} is specified.{p_end}


    {hline}
	
		
{marker examples}{...}
{title:Examples}

{p 4 8}Standard use of rdms for cumulative cutoffs{p_end}
{p 8 8}{cmd:. rdms yvar xvar, c(cvar)}{p_end}

{p 4 8}rdms with plot{p_end}
{p 8 8}{cmd:. rdms yvar xvar, c(cvar) plot}{p_end}

{p 4 8}Standard use of rdms for multiple scores{p_end}
{p 8 8}{cmd:. rdms yvar xvar1 xvar2 treatvar, c(cvar)}{p_end}


{marker saved_results}{...}
{title:Saved results}

{p 4 8}{cmd:rdms} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}bias-corrected coefficient vector{p_end}
{synopt:{cmd:e(V)}}robust variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(coefs)}}conventional coefficient vector{p_end}
{synopt:{cmd:e(pv_rb)}}robust p-value vector{p_end}
{synopt:{cmd:e(CI_rb)}}bias-corrected confidence intervals{p_end}
{synopt:{cmd:e(H)}}vector of bandwidths at each side of each cutoff{p_end}
{synopt:{cmd:e(sampsis)}}vector of sample sizes at each side of each cutoff{p_end}


{marker references}{...}
{title:References}

{p 4 8}Calonico, S., M. D. Cattaneo, M. H. Farrell, and R. Titiunik. 2017.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2017_Stata.pdf":rdrobust: Software for Regression Discontinuity Designs}.{p_end}
{p 8 8}{it:Stata Journal} 17(2): 372-404.{p_end}

{p 4 8}Calonico, S., M. D. Cattaneo, and R. Titiunik. 2014.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_Stata.pdf":Robust Data-Driven Inference in the Regression-Discontinuity Design}.{p_end}
{p 8 8}{it:Stata Journal} 14(4): 909-946.{p_end}

{p 4 8}Calonico, S., M. D. Cattaneo, and R. Titiunik. 2015.
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_R.pdf":rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs}.{p_end}
{p 8 8}{it:R Journal} 7(1): 38-51.{p_end}

{p 4 8}Cattaneo, M. D., L. Keele, R. Titiunik, and G. Vazquez-Bare. 2016.
{browse "https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2016_JOP.pdf":Interpreting Regression Discontinuity Designs with Multiple Cutoffs}.{p_end}
{p 8 8}{it:Journal of Politics} 78(4): 1229-1248.{p_end}

{p 4 8}Cattaneo, M. D., L. Keele, R. Titiunik, and G. Vazquez-Bare. 2021.
{browse "https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2021_JASA.pdf":Extrapolating Treatment Effects in Multi-Cutoff Regression Discontinuity Designs}.{p_end}
{p 8 8}{it:Journal of American Statistical Association}, forthcoming.{p_end}

{p 4 8}Cattaneo, M. D., R. Titiunik, and G. Vazquez-Bare. 2020.
{browse "https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf":Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores}.{p_end}
{p 8 8}{it:Stata Journal}, forthcoming.{p_end}

{p 4 8}Keele, L., and R. Titiunik. 2015.
{browse "https://rdpackages.github.io/references/Keele-Titiunik_2015_PA.pdf":Geographic Boundaries as Regression Discontinuities}.{p_end}
{p 8 8}{it:Political Analysis} 23(1): 127-155.{p_end}


{marker authors}{...}
{title:Authors}

{p 4 8}Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.{p_end}

{p 4 8}Rocio Titiunik, Princeton University, Princeton, NJ.
{browse "mailto:titiunik@princeton.edu":titiunik@princeton.edu}.{p_end}

{p 4 8}Gonzalo Vazquez-Bare, UC Santa Barbara, Santa Barbara, CA.
{browse "mailto:gvazquez@econ.ucsb.edu":gvazquez@econ.ucsb.edu}.{p_end}



