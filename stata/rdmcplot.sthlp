{smcl}
{* *! version 0.5 2020-08-22}{...}
{viewerjumpto "Syntax" "rdms##syntax"}{...}
{viewerjumpto "Description" "rdms##description"}{...}
{viewerjumpto "Options" "rdms##options"}{...}
{viewerjumpto "Examples" "rdms##examples"}{...}
{viewerjumpto "Saved results" "rdms##saved_results"}{...}
{viewerjumpto "References" "rdmc##references"}{...}
{viewerjumpto "Authors" "rdmc##authors"}{...}

{title:Title}

{p 4 8}{cmd:rdmcplot} {hline 2} RD Plots for Regression Discontinuity Designs with Multiple Cutoffs.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdmcplot} {it:depvar} {it:runvar} {ifin}{cmd:,}
{cmd:{opt c:var}(}{it:string}{cmd:)} 
[
{cmd:{opt nbins:var}(}{it:string}{cmd:)}
{cmd:{opt nbinsr:ightvar}(}{it:string}{cmd:)}
{cmd:{opt binselect:var}(}{it:string}{cmd:)} 
{cmd:{opt scale:var}(}{it:string}{cmd:)} 
{cmd:{opt scaler:ightvar}(}{it:string}{cmd:)} 
{cmd:{opt support:var}(}{it:string}{cmd:)} 
{cmd:{opt supportr:ightvar}(}{it:string}{cmd:)} 
{cmd:{opt p:var}(}{it:string}{cmd:)} 
{cmd:{opt h:var}(}{it:string}{cmd:)} 
{cmd:{opt hr:ightvar}(}{it:string}{cmd:)} 
{cmd:{opt kernel:var}(}{it:string}{cmd:)} 
{cmd:{opt weights:var}(}{it:string}{cmd:)} 
{cmd:{opt covs:var}(}{it:string}{cmd:)} 
{cmd:{opt covseval:var}(}{it:string}{cmd:)} 
{cmd:{opt covsdrop:var}(}{it:string}{cmd:)} 
{cmd:{opt binsopt:var}(}{it:string}{cmd:)} 
{cmd:{opt lineopt:var}(}{it:string}{cmd:)} 
{cmd:{opt xlineopt:var}(}{it:string}{cmd:)} 
{cmd:ci(}{it:cilevel}{cmd:)} 
{cmd:nobins} 
{cmd:nopoly}
{cmd:noxline}  
{cmd:nodraw}
{cmd:genvars} 
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rdmcplot} plots estimated regression functions at each cutoff in regression discontinuity (RD) designs with multiple cutoffs.
For methodological background see
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_JASA.pdf":Calonico, Cattaneo and Titiunik (2015a)}, 
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_JASA.pdf":Calonico, Cattaneo and Titiunik (2015a)}, 
{browse "https://rdpackages.github.io/references/Keele-Titiunik_2015_PA.pdf":Keele and Titiunik (2015)},
{browse "https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2016_JOP.pdf":Cattaneo, Keele, Titiunik and Vazquez-Bare (2016)}, and
{browse "https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2021_JASA.pdf":Cattaneo, Keele, Titiunik and Vazquez-Bare (2021)}.{p_end}

{p 8 8}Companion commands are: {help rdmc:rdmc} for multi-cutoff RD estimation and inference, and {help rdms:rdms} for multi-score RD estimation and inference.{p_end}

{p 8 8}A detailed introduction to this command is given in
{browse "https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf": Cattaneo, Titiunik and Vazquez-Bare (2020)}.{p_end}

{p 8 8}This command employs the Stata (and R) package {help rdrobust:rdrobust} for underlying calculations. See
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_Stata.pdf":Calonico, Cattaneo and Titiunik (2014)},
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_R.pdf":Calonico, Cattaneo and Titiunik (2015b)},
and
{browse "https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2017_Stata.pdf":Calonico, Cattaneo, Farrell and Titiunik (2017)}
for more details.{p_end}

{p 4 8}Related Stata and R packages useful for inference in RD designs are described in the following website:{p_end}

{p 8 8}{browse "https://rdpackages.github.io/":https://rdpackages.github.io/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimand}

{p 4 8}{cmd:{opt c:var}(}{it:string}{cmd:)} specifies the numeric variable containing the RD cutoff for {it:indepvar} for each unit in the sample.{p_end}

{dlgtab:Bin Selection}

{p 4 8}{cmd:{opt nbins:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the number of bins for {cmd:rdplot}. 
When {cmd:nbinsrightvar} is specified, {cmd:nbinsvar} indicates the number of bins to the left of the cutoff.
When {cmd:nbinsrightvar} is not specified, the same number of bins is used at each side.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt nbinsr:ightvar}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the number of bins to the right of the cutoff for {cmd:rdplot}. 
When {cmd:nbinsrightvar} is not specified, the number of bins in {cmd:nbinsvar} is used at each side.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt binselect:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the bins selection method for {cmd:rdplot}.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt scale:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the scale for {cmd:rdplot}.
When {cmd:scalerightvar} is specified, {cmd:nbinsvar} indicates the scale to the left of the cutoff.
When {cmd:scalerightvar} is not specified, the same scale is used at each side. See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt scaler:ightvar}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the scale to the right of the cutoff for {cmd:rdplot}.
When {cmd:scalerightvar} is not specified, the scale in {cmd:scalevar} is used at each side.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt support:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the support for {cmd:rdplot}. 
When {cmd:supportrightvar} is specified, {cmd:supportvar} indicates the support to the left of the cutoff.
When {cmd:supportrightvar} is not specified, the same support is used at each side.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt supportr:ightvar}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the support to the right of the cutoff for {cmd:rdplot}. 
When {cmd:supportrightvar} is not specified, the support in {cmd:supportvar} are used at each side.
See {help rdplot:rdplot} for details.{p_end}

{dlgtab:Polynomial Fit}

{p 4 8}{cmd:{opt p:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the order of the polynomials for {cmd:rplot}. 
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt h:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the bandwidths for {cmd:rdplot}. 
When {cmd:hrightvar} is specified, {cmd:hvar} indicates the bandwidth to the left of the cutoff.
When {cmd:hrightvar} is not specified, the same bandwidths are used at each side.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt hr:ightvar}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the bandwidths to the right of the cutoff for {cmd:rdplot}. 
When {cmd:hrightvar} is not specified, the bandwidths in {cmd:hvar} are used at each side.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt kernel:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the kernels for {cmd:rdplot}.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt weights:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the weights for {cmd:rdplot}.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt covs:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the covariates for {cmd:rdplot}.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt covseval:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies the evaluation points for additional covariates.
See {help rdplot:rdplot} for details.{p_end}

{p 4 8}{cmd:{opt covsdrop:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies whether collinear covariates should be dropped.
See {help rdplot:rdplot} for details.{p_end}

{dlgtab:Plot}

{p 4 8}{cmd:{opt binsopt:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies options for the bins plots.{p_end}

{p 4 8}{cmd:{opt lineopt:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies options for the polynomial plots.{p_end}

{p 4 8}{cmd:{opt xlineopt:var}(}{it:string}{cmd:)} a variable of length equal to the number of different cutoffs that specifies options for the vertical lines indicating the cutoffs.{p_end}

{p 4 8}{cmd:ci(}{it:cilevel}{cmd:)} adds confidence intervals of level {it:cilevel} to the plot.{p_end}

{p 4 8}{cmd:nobins} omits the bins plot.{p_end}

{p 4 8}{cmd:nopoly} omits the polynomial curve plot.{p_end}

{p 4 8}{cmd:noxline} omits the vertical lines indicating the cutoffs.{p_end}

{p 4 8}{cmd:nodraw} omits the plot.{p_end}

{dlgtab:Generate Variables}

{p 4 8}{cmd:genvars} generates variables to replicate plots by hand. Variable labels indicate the corresponding cutoff.{p_end}
{p 8 8}{cmd:rdmcplot_hat_y_{it:c}} predicted value of the outcome variable given by the global polynomial estimator in cutoff number {it:c}.{p_end}
{p 8 8}{cmd:rdmcplot_mean_x_{it:c}} sample mean of the running variable within the corresponding bin for each observation in cutoff number {it:c}.{p_end}
{p 8 8}{cmd:rdmcplot_mean_y_{it:c}} sample mean of the outcome variable within the corresponding bin for each observation in cutoff number {it:c}.{p_end}
{p 8 8}{cmd:rdmcplot_ci_l_{it:c}} lower end value of the confidence interval for the sample mean of the outcome variable within the corresponding bin for each observation in cutoff number {it:c}.{p_end}
{p 8 8}{cmd:rdmcplot_ci_r_{it:c}} upper end value of the confidence interval for the sample mean of the outcome variable within the corresponding bin for each observation in cutoff number {it:c}.{p_end}

    {hline}
	
		
{marker examples}{...}
{title:Examples}

{p 4 8}Standard use of rdmcplot{p_end}
{p 8 8}{cmd:. rdmcplot yvar xvar, c(cvar)}{p_end}

{p 4 8}rdmcplot without bins plot{p_end}
{p 8 8}{cmd:. rdmcplot yvar xvar, c(cvar) nobins}{p_end}


{marker saved_results}{...}
{title:Saved results}

{p 4 8}{cmd:rdmcplot} saves the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(p)}}order of the polynomial{p_end}
{synopt:{cmd:r(cnum)}}number of cutoffs{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(cvar)}}cutoff variable{p_end}
{synopt:{cmd:r(clist)}}cutoff list{p_end}


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
{p 8 8}{it:Journal of Politics} 78(4): 1229-124.{p_end}

{p 4 8}Cattaneo, M. D., L. Keele, R. Titiunik, and G. Vazquez-Bare. 2021.
{browse "https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2021_JASA.pdf":Extrapolating Treatment Effects in Multi-Cutoff Regression Discontinuity Designs}.{p_end}
{p 8 8}{it:Journal of American Statistical Association}, forthcoming.{p_end}

{p 4 8}Cattaneo, M. D., R. Titiunik, and G. Vazquez-Bare. 2020.
{browse "https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf":Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores}.{p_end}
{p 8 8}{it:Working paper}.{p_end}

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


