### RDMULTI

The rdmulti package provides Stata and R implementation of RD plots, estimation, inference and extrapolation methods for RD designs with multiple cutoffs and multiple scores.

This work was supported in part by the National Science Foundation through grants [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561).

## Stata Implementation

To install/update in Stata type:
```
net install rdmulti, from(https://raw.githubusercontent.com/rdpackages/rdmulti/master/stata) replace
```

- Help: [rdmc](stata/rdmc.pdf), [rdmcplot](stata/rdmcplot.pdf), [rdms](stata/rdms.pdf).

- Replication: [do-file](stata/rdmulti_illustration.do), [rdmcplot illustration](stata/rdmcplot_illustration.do), [dataset1](stata/simdata_multic.dta), [dataset2](stata/simdata_cumul.dta), [dataset3](stata/simdata_multis.dta).

## R Implementation

To install/update in R type:
```
install.packages('rdmulti')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rdmulti/rdmulti.pdf), [CRAN repository](https://cran.r-project.org/package=rdmulti).

- Replication: [R-script](R/rdmulti_illustration.R), [rdmcplot illustration](R/rdmcplot_illustration.R), [dataset1](stata/simdata_multic.csv), [dataset2](stata/simdata_cumul.csv), [dataset3](stata/simdata_multis.csv), [R illustration](R/rdmulti_illustration.pdf)

## References

For overviews and introductions, see [rdpackages website]().

### Software and Implementation

- Calonico, Cattaneo and Titiunik (2014): [Robust Data-Driven Inference in the Regression-Discontinuity Design](references/Calonico-Cattaneo-Titiunik_2014_Stata.pdf), _Stata Journal_ 14(4): 909-946.

- Calonico, Cattaneo and Titiunik (2015): [rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs](references/Calonico-Cattaneo-Titiunik_2015_R.pdf), _R Journal_ 7(1): 38-51.

- Calonico, Cattaneo, Farrell and Titiunik (2017): [rdrobust: Software for Regression Discontinuity Designs](references/Calonico-Cattaneo-Farrell-Titiunik_2017_Stata.pdf), _Stata Journal_ 17(2): 372-404.

### Technical and Methodological

- Calonico, Cattaneo and Titiunik (2014): [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf), _Econometrica_ 82(6): 2295-2326. [Supplemental Appendix](references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf).

- Calonico, Cattaneo and Titiunik (2015): [Optimal Data-Driven Regression Discontinuity Plots](references/Calonico-Cattaneo-Titiunik_2015_JASA.pdf), _Journal of the American Statistical Association_ 110(512): 1753-1769. [Supplemental Appendix](references/Calonico-Cattaneo-Titiunik_2015_JASA--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2018): [On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference](references/Calonico-Cattaneo-Farrell_2018_JASA.pdf), _Journal of the American Statistical Association_ 113(522): 767-779. [Supplemental Appendix](references/Calonico-Cattaneo-Farrell_2018_JASA--Supplement.pdf).

- Calonico, Cattaneo, Farrell and Titiunik (2019): [Regression Discontinuity Designs Using Covariates](references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf), _Review of Economics and Statistics_ 101(3): 442-451. [Supplemental Appendix](references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2020): [Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs](references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf), _Econometrics Journal_ 23(2): 192-210. [Supplemental Appendix](references/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2020): [Coverage Error Optimal Confidence Intervals for Local Polynomial Regression](references/Calonico-Cattaneo-Farrell_2020_CEopt.pdf), working paper. [Supplemental Appendix](references/Calonico-Cattaneo-Farrell_2020_CEopt--Supplement.pdf).

<br>
<br>
