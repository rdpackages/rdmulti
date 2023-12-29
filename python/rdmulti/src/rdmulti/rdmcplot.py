#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdrobust import rdplot

def rdmcplot(Y, X, C, nbinsmat=None, binselectvec=None, scalevec=None,
             supportmat=None, pvec=None, hmat=None, kernelvec=None,
             weightsvec=None, covs_mat=None, covs_list=None, covs_evalvec=None,
             covs_dropvec=None, ci=None, col_bins=None, pch_bins=None,
             col_poly=None, lty_poly=None, col_xline=None, lty_xline=None,
             nobins=False, nopoly=False, noxline=False, nodraw=False):
    
    """
    RD plots with multiple cutoffs

    rdmcplot() generates RD plots with multiple cutoffs.

    Author:
    Matias Cattaneo, Princeton University. Email: cattaneo@princeton.edu
    Rocio Titiunik, Princeton University. Email: titiunik@princeton.edu
    Ricardo Masini, Princeton University. Email: rmasini@princeton.edu
    Gonzalo Vazquez-Bare, UC Santa Barbara. Email: gvazquez@econ.ucsb.edu

    References:
    Cattaneo, M.D., R. Titiunik and G. Vazquez-Bare. (2020). [Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores](https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf). Stata Journal, forthcoming.

    Parameters:
        Y: Outcome variable.
        X: Running variable.
        C: Cutoff variable.
        nbinsmat: Matrix of cutoff-specific number of bins. See rdplot() for details.
        binselectvec: Vector of cutoff-specific bins selection method. See rdplot() for details.
        scalevec: Vector of cutoff-specific scale factors. See rdplot() for details.
        supportmat: Matrix of cutoff-specific support conditions. See rdplot() for details.
        pvec: Vector of cutoff-specific polynomial orders. See rdplot() for details.
        hmat: Matrix of cutoff-specific bandwidths. See rdplot() for details.
        kernelvec: Vector of cutoff-specific kernels. See rdplot() for details.
        weightsvec: Vector of cutoff-specific weights. See rdplot() for details.
        covs_mat: Matrix of covariates. See rdplot() for details.
        covs_list: List of covariates to be used in each cutoff.
        covs_evalvec: Vector indicating the evaluation point for additional covariates. See rdrobust() for details.
        covs_dropvec: Vector indicating whether collinear covariates should be dropped at each cutoff. See rdrobust() for details.
        ci: Adds confidence intervals of the specified level to the plot. See rdrobust() for details.
        col_bins: Vector of colors for bins.
        pch_bins: Vector of characters (pch) type for bins.
        col_poly: Vector of colors for polynomial curves.
        lty_poly: Vector of lty for polynomial curves.
        col_xline: Vector of colors for vertical lines.
        lty_xline: Vector of lty for vertical lines.
        nobins: Omits bins plot.
        nopoly: Omits polynomial curve plot.
        noxline: Omits vertical lines indicating the cutoffs.
        nodraw: Omits plot.

    Returns:
        clist: List of cutoffs.
        cnum: Number of cutoffs.
        X0: Matrix of X values for control units.
        X1: Matrix of X values for treated units.
        Yhat0: Estimated polynomial for control units.
        Yhat1: Estimated polynomial for treated units.
        Xmean: Bin average of X values.
        Ymean: Bin average for Y values.
        CI_l: Lower end of confidence intervals.
        CI_r: Upper end of confidence intervals.
        cfail: Cutoffs where rdrobust() encountered problems.

    Examples:
    # Toy dataset
    X = np.random.uniform(0, 100, 1000)
    C = np.concatenate([np.repeat(33, 500), np.repeat(66, 500)])
    Y = (1 + X + (X >= C)) * (C == 33) + (0.5 + 0.5 * X + 0.8 * (X >= C)) * (C == 66) + np.random.normal(size=1000)
    # rdmcplot with standard syntax
    tmp = rdmcplot(Y, X, C)
    """
  
    #################################################################
    # Setup and error checking
    #################################################################

    try: C = np.array(C, dtype = 'float64')
    except: raise Exception('C has to be numeric')

    if np.max(C) >= np.max(X) or np.min(C) <= np.min(X):
        raise ValueError("Cutoff variable outside range of running variable")

    clist = np.sort(np.unique(C))
    cnum = len(clist)

    D = np.asarray(X >= C, dtype=float)

    if pvec is None:
        pvec = np.repeat(4, cnum)

    if hmat is not None:
        if np.isscalar(hmat): hmat = np.full((cnum, 2), hmat)
        elif len(hmat)==cnum : hmat = np.column_stack((hmat,hmat))
        elif len(hmat)==2: hmat = np.tile(hmat,(cnum,1))
        haux = hmat
    else:
        hmat = np.repeat(None,cnum)
        haux = np.full((cnum, 2), np.inf)

    if nbinsmat is None:
        nbinsmat = np.repeat(None,cnum)
    else:
        if nbinsmat.ndim == 1:
            nbinsmat = np.repeat(nbinsmat, 2).reshape(cnum, 2)

    if supportmat is None:
        supportmat = np.repeat(None,cnum)
    else:
        if supportmat.ndim == 1:
            supportmat = np.repeat(supportmat, 2).reshape(cnum, 2)

    if binselectvec is None:
        binselectvec = np.repeat('esmv', cnum)

    if scalevec is None:
        scalevec = np.repeat(1, cnum)

    if kernelvec is None:
        kernelvec = np.repeat('uni', cnum)

    if covs_evalvec is None:
        covs_evalvec = np.repeat(0, cnum)

    if covs_dropvec is None:
        covs_dropvec = np.repeat(True, cnum)
    
    if weightsvec is None:  
        weightsvec = np.repeat(None,cnum)

    if covs_mat is not None:
        covs_mat = np.asarray(covs_mat)
        if covs_list is not None:
            if len(covs_list) != cnum:
                raise ValueError("Elements in covs_list should equal number of cutoffs")

    X0 = np.full((len(Y), cnum), np.nan)
    X1 = np.full((len(Y), cnum), np.nan)
    YHAT0 = np.full((len(Y), cnum), np.nan)
    YHAT1 = np.full((len(Y), cnum), np.nan)
    XMEAN = np.full((len(Y), cnum), np.nan)
    YMEAN = np.full((len(Y), cnum), np.nan)
    CI_l = np.full((len(Y), cnum), np.nan)
    CI_r = np.full((len(Y), cnum), np.nan)
    Cfail = np.array([])

    #################################################################
    # Construct variables for plots
    #################################################################

    count = 0
    count_fail = 0

    for c in clist:
        yc = Y[(C == c) & (X <= c + haux[count, 1]) & (X >= c - haux[count, 0])]
        xc = X[(C == c) & (X <= c + haux[count, 1]) & (X >= c - haux[count, 0])]
        dc = D[(C == c) & (X <= c + haux[count, 1]) & (X >= c - haux[count, 0])]
        yc0 = yc[dc == 0]
        yc1 = yc[dc == 1]
        xc0 = xc[dc == 0]
        xc1 = xc[dc == 1]

        if covs_mat is not None:
            covs_mat_c = covs_mat[(C == c) & (X <= c + haux[count, 1]) & (X >= c - haux[count, 0]), :]
            if covs_list is not None:
                covs_aux = covs_mat_c[:, covs_list[count]]
            else:
                covs_aux = covs_mat_c
        else:
            covs_aux = None

        rdplot(yc, xc, c=c,
                        nbins=nbinsmat[count],
                        binselect=binselectvec[count],
                        scale=scalevec[count],
                        support=supportmat[count],
                        p=pvec[count],
                        h=hmat[count],
                        kernel=kernelvec[count],
                        weights=weightsvec[count],
                        covs=covs_aux,
                        covs_eval=covs_evalvec[count],
                        covs_drop=covs_dropvec[count],
                        ci=ci,
                        hide=True)

        try:
            aux = rdplot(yc, xc, c=c,
                        nbins=nbinsmat[count],
                        binselect=binselectvec[count],
                        scale=scalevec[count],
                        support=supportmat[count],
                        p=pvec[count],
                        h=hmat[count],
                        kernel=kernelvec[count],
                        weights=weightsvec[count],
                        covs=covs_aux,
                        covs_eval=covs_evalvec[count],
                        covs_drop=covs_dropvec[count],
                        ci=ci,
                        hide=True)
        except:
            Cfail = np.append(Cfail, c)
            count_fail += 1
        
        else:
            xmean = aux.vars_bins.iloc[:,1].values
            ymean = aux.vars_bins.iloc[:,2].values
            
            xmean = xmean[~np.isnan(xmean)]
            ymean = ymean[~np.isnan(ymean)]
            
            x0 = aux.vars_poly.iloc[:,0].values
            smaller = x0 < c
            greater = x0 > c
            x0 = x0[smaller]
            
            yhat0 = aux.vars_poly.iloc[:,1].values
            yhat0 = yhat0[smaller]
            
            x1 = aux.vars_poly.iloc[:,0].values
            x1 = x1[greater]
            
            yhat1 = aux.vars_poly.iloc[:,1].values
            yhat1 = yhat1[greater]
            
            XMEAN[:len(xmean), count] = xmean
            YMEAN[:len(ymean), count] = ymean
            X0[:len(x0), count] = x0
            X1[:len(x1), count] = x1
            YHAT0[:len(yhat0), count] = yhat0
            YHAT1[:len(yhat1), count] = yhat1

            if ci is not None:
                ci_l = aux.vars_bins.iloc[:, 7]
                ci_r = aux.vars_bins.iloc[:, 8]
                ci_l = np.resize(ci_l, len(Y))
                ci_r = np.resize(ci_r, len(Y))
                CI_l.iloc[:, count] = ci_l
                CI_r.iloc[:, count] = ci_r
        finally:
            count += 1


    Xmean = pd.DataFrame(XMEAN)
    Ymean = pd.DataFrame(YMEAN)
    X0 = pd.DataFrame(X0)
    X1 = pd.DataFrame(X1)
    Yhat0 = pd.DataFrame(YHAT0)
    Yhat1 = pd.DataFrame(YHAT1)
    CI_l = pd.DataFrame(CI_l)
    CI_r = pd.DataFrame(CI_r)

    Xmean = Xmean.iloc[:, :np.sum(~np.isnan(Xmean.sum()))]
    Xmean.columns = [f"Xmean{i + 1}" for i in range(Xmean.shape[1])]
    Ymean = Ymean.iloc[:, :np.sum(~np.isnan(Ymean.sum()))]
    Ymean.columns = [f"Ymean{i + 1}" for i in range(Ymean.shape[1])]
    X0 = X0.iloc[:, :np.sum(~np.isnan(X0.sum()))]
    X0.columns = [f"X0_{i + 1}" for i in range(X0.shape[1])]
    X1 = X1.iloc[:, :np.sum(~np.isnan(X1.sum()))]
    X1.columns = [f"X1_{i + 1}" for i in range(X1.shape[1])]
    Yhat0 = Yhat0.iloc[:, :np.sum(~np.isnan(Yhat0.sum()))]
    Yhat0.columns = [f"Yhat0_{i + 1}" for i in range(Yhat0.shape[1])]
    Yhat1 = Yhat1.iloc[:, :np.sum(~np.isnan(Yhat1.sum()))]
    Yhat1.columns = [f"Yhat1_{i + 1}" for i in range(Yhat1.shape[1])]

    if ci is not None:
        CI_l = CI_l.iloc[:, :np.sum(~np.isnan(CI_l.sum()))]
        CI_l.columns = [f"CI_l_{i + 1}" for i in range(CI_l.shape[1])]
        CI_r = CI_r.iloc[:, :np.sum(~np.isnan(CI_r.sum()))]
        CI_r.columns = [f"CI_r_{i + 1}" for i in range(CI_r.shape[1])]

    #################################################################
    # Plots
    #################################################################
    
    colorlist = ['darkblue', 'darkred', 'darkgreen', 'darkorange', 'gray50', 'khaki4', 'brown3', 'blue', 'darkgoldenrod4', 'cyan4']

    if col_bins is None:
        col_bins = colorlist

    if pch_bins is None:
        pch_bins = [1] * cnum

    if col_poly is None:
        col_poly = colorlist

    if lty_poly is None:
        lty_poly = ['solid'] * cnum

    if col_xline is None:
        col_xline = colorlist

    if lty_xline is None:
        lty_xline = ['dashed'] * cnum

    rdmc_plot = plt.figure().gca()
    rdmc_plot.set_xlabel('Running variable')
    rdmc_plot.set_ylabel('Outcome')

    if not nobins:
        for c in range(cnum):
            rdmc_plot.plot(Xmean.iloc[:, c].values,
                            Ymean.iloc[:, c].values,
                            marker='o',
                            mfc='none',
                            color=col_bins[c],
                            linestyle='None')

    if ci is not None:
        for c in range(cnum):
            rdmc_plot.errorbar(Xmean.iloc[:, c].values,
                                Ymean.iloc[:, c].values,
                                yerr=[CI_l.iloc[:, c].values,CI_r.iloc[:, c].values],
                                color=col_bins[c],
                                linestyle='-', 
                                linewidth=1)

    if not nopoly:
        for c in range(cnum):
            rdmc_plot.plot(X0.iloc[:, c].values,
                            Yhat0.iloc[:, c].values,
                            color=col_poly[c],
                            linestyle=lty_poly[c])
            rdmc_plot.plot(X1.iloc[:, c].values,
                            Yhat1.iloc[:, c].values,
                            color=col_poly[c],
                            linestyle=lty_poly[c])

    if not noxline:
        for c in range(cnum):
            rdmc_plot.axvline(x=clist[c], color=col_xline[c], linestyle=lty_xline[c])

    if not nodraw:
        plt.show()

    if count_fail > 0:
        print("rdplot() could not run in one or more cutoffs.")

    #################################################################
    # Return values
    #################################################################

    if ci is None:
        output = {
            'clist': clist,
            'cnum': cnum,
            'X0': X0,
            'X1': X1,
            'Yhat0': Yhat0,
            'Yhat1': Yhat1,
            'Xmean': Xmean,
            'Ymean': Ymean,
            'rdmc_plot': rdmc_plot,
            'cfail': Cfail
        }
    else:
        output = {
            'clist': clist,
            'cnum': cnum,
            'X0': X0,
            'X1': X1,
            'Yhat0': Yhat0,
            'Yhat1': Yhat1,
            'Xmean': Xmean,
            'Ymean': Ymean,
            'rdmc_plot': rdmc_plot,
            'CI_l': CI_l,
            'CI_r': CI_r,
            'cfail': Cfail
        }
    return output