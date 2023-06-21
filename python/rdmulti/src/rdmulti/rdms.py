#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas  as pd
from rdrobust import rdrobust

def rdms(Y, X, C, X2=None, zvar=None, C2=None, rangemat=None, xnorm=None,
         fuzzy=None, derivvec=None, pooled_opt=None, pvec=None, qvec=None,
         hmat=None, bmat=None, rhovec=None, covs_mat=None, covs_list=None,
         covs_dropvec=None, kernelvec=None, weightsvec=None,
         bwselectvec=None, scaleparvec=None, scaleregulvec=None,
         masspointsvec=None, bwcheckvec=None, bwrestrictvec=None,
         stdvarsvec=None, vcevec=None, nnmatchvec=None, cluster=None,
         level=95, plot=False, conventional=False):
    
    """
    Analysis of RD designs with cumulative cutoffs or two running variables

    rdms() analyzes RD designs with cumulative cutoffs or two running variables.

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
        C: Vector of cutoffs.
        X2: If specified, second running variable.
        zvar: If X2 is specified, treatment indicator.
        C2: If specified, second vector of cutoffs.
        rangemat: Matrix of cutoff-specific ranges for the running variable.
        xnorm: Normalized running variable to estimate pooled effect.
        fuzzy: Specifies a fuzzy design. See rdrobust() for details.
        derivvec: Vector of cutoff-specific order of derivatives. See rdrobust() for details.
        pooled_opt: Options to be passed to rdrobust() to calculate pooled estimand.
        pvec: Vector of cutoff-specific polynomial orders. See rdrobust() for details.
        qvec: Vector of cutoff-specific polynomial orders for bias estimation. See rdrobust() for details.
        hmat: Matrix of cutoff-specific bandwidths. See rdrobust() for details.
        bmat: Matrix of cutoff-specific bandwidths for bias estimation. See rdrobust() for details.
        rhovec: Vector of cutoff-specific values of rho. See rdrobust() for details.
        covs_mat: Matrix of covariates. See rdplot() for details.
        covs_list: List of covariates to be used in each cutoff.
        covs_dropvec: Vector indicating whether collinear covariates should be dropped at each cutoff. See rdrobust() for details.
        kernelvec: Vector of cutoff-specific kernels. See rdrobust() for details.
        weightsvec: Vector of cutoff-specific weights. See rdrobust() for details.
        bwselectvec: Vector of cutoff-specific bandwidth selection methods. See rdrobust() for details.
        scaleparvec: Vector of cutoff-specific scale parameters. See rdrobust() for details.
        scaleregulvec: Vector of cutoff-specific scale regularization parameters. See rdrobust() for details.
        masspointsvec: Vector indicating how to handle repeated values at each cutoff. See rdrobust() for details.
        bwcheckvec: Vector indicating the value of bwcheck at each cutoff. See rdrobust() for details.
        bwrestrictvec: Vector indicating whether computed bandwidths are restricted to the range or runvar at each cutoff. See rdrobust() for details.
        stdvarsvec: Vector indicating whether variables are standardized at each cutoff. See rdrobust() for details.
        vcevec: Vector of cutoff-specific variance-covariance estimation methods. See rdrobust() for details.
        nnmatchvec: Vector of cutoff-specific nearest neighbors for variance estimation. See rdrobust() for details.
        cluster: Cluster ID variable. See rdrobust() for details.
        level: Confidence level for confidence intervals. See rdrobust() for details.
        plot: Plots cutoff-specific and pooled estimates.
        conventional: Reports conventional, instead of robust-bias corrected, p-values and confidence intervals.

    Returns:
        B: Vector of bias-corrected coefficients.
        V: Variance-covariance matrix of the estimators.
        Coefs: Vector of conventional coefficients.
        Nh: Vector of sample sizes within bandwidth at each cutoff.
        CI: Bias corrected confidence intervals.
        H: Bandwidth used at each cutoff.
        Pv: Vector of robust p-values.

    Examples:
    # Toy dataset: cumulative cutoffs
    X = np.random.uniform(0, 100, 1000)
    C = [33, 66]
    Y = (1 + X) * (X < C[0]) + (0.8 + 0.8 * X) * (X >= C[0] & X < C[1]) + (1.2 + 1.2 * X) * (X >= C[1]) + np.random.normal(size=1000)
    # rdms: basic syntax
    tmp = rdms(Y, X, C)
    """

    # Setup and error checking
    if X2 is not None and zvar is None:
        raise ValueError("Need to specify zvar when X2 is specified")
    if X2 is not None and C2 is None:
        raise ValueError("Need to specify C2 if X2 is specified")

    cnum = len(C)

    if C2 is not None:
        if cnum != len(C2):
            raise ValueError("Cutoff coordinates incorrectly specified")

    if rangemat is not None:
        if len(rangemat.shape) == 1:
            rangemat = rangemat.reshape(cnum, 2)
    else:
        rangemat = np.column_stack((np.full(cnum, -np.inf), np.full(cnum, np.inf)))

    if hmat is not None:
        if np.isscalar(hmat): hmat = np.full((cnum, 2), hmat)
        elif len(hmat)==cnum : hmat = np.column_stack((hmat,hmat))
        elif len(hmat)==2: hmat = np.tile(hmat,(cnum,1))   
    else: hmat = np.repeat(None,cnum)
    
    if bmat is not None:
        if np.isscalar(bmat): hmat = np.full((cnum, 2), bmat)
        elif len(bmat)==cnum : bmat = np.column_stack((bmat,bmat))
        elif len(bmat)==2: bmat = np.tile(bmat,(cnum,1))
    else: bmat = np.repeat(None,cnum)

    if covs_dropvec is None:
        covs_dropvec = np.full(cnum, True)
    if kernelvec is None:
        kernelvec = np.full(cnum, 'tri')
    if bwselectvec is None:
        bwselectvec = np.full(cnum, 'mserd')
    if vcevec is None:
        vcevec = np.full(cnum, 'nn')
    if nnmatchvec is None:
        nnmatchvec = np.full(cnum, 3)
    if scaleparvec is None:
        scaleparvec = np.full(cnum, 1)
    if scaleregulvec is None:
        scaleregulvec = np.full(cnum, 1)
    if masspointsvec is None:
        masspointsvec = np.full(cnum, 'adjust')
    if bwrestrictvec is None:
        bwrestrictvec = np.full(cnum, True)
    if stdvarsvec is None:
        stdvarsvec = np.full(cnum, False)
    if covs_mat is not None:
        covs_mat = np.array(covs_mat)
        if covs_list is not None:
            if len(covs_list) != cnum:
                raise ValueError("Elements in covs_list should equal number of cutoffs")
    
    
    if derivvec is None: derivvec = np.repeat(None,cnum)
    if pvec is None: pvec = np.repeat(None,cnum)
    if qvec is None: qvec = np.repeat(None,cnum)
    if rhovec is None: rhovec = np.repeat(None,cnum)
    if bwcheckvec is None: bwcheckvec = np.repeat(None,cnum)
   

    B = np.full((1, cnum+1), np.nan)
    V = np.full((1, cnum+1), np.nan)
    Coefs = np.full((1, cnum+1), np.nan)
    V_cl = np.full((1, cnum+1), np.nan)
    Nh = np.full((2, cnum+1), np.nan)
    CI = np.full((2, cnum+1), np.nan)
    CI_cl = np.full((2, cnum+1), np.nan)
    Pv = np.full((1, cnum+1), np.nan)
    Pv_cl = np.full((1, cnum+1), np.nan)
    H = np.full((2, cnum+1), np.nan)
    Bbw = np.full((2, cnum+1), np.nan)

    c_disp = []

    #################################################################
    # Calculate cutoff-specific estimates
    #################################################################
   
    if X2 is None:
        Rc = rangemat - C.reshape(cnum,-1)
        for c in range(cnum):
            xc = X - C[c]
            yc = Y[(xc >= Rc[c, 0]) & (xc <= Rc[c, 1])]
            xc = xc[(xc >= Rc[c, 0]) & (xc <= Rc[c, 1])]

            weightsc = None
            if weightsvec is not None:
                weightaux = weightsvec[c]
                weightsc = eval(weightaux + '[xc >= Rc[c, 0] & xc <= Rc[c, 1]]')

            covs_aux = None
            if covs_mat is not None:
                covs_mat_c = covs_mat[(xc >= Rc[c, 0]) & (xc <= Rc[c, 1]), :]
                if covs_list is not None:
                    covs_aux = covs_mat_c[:, covs_list[c]]
                else:
                    covs_aux = covs_mat_c

            if cluster is not None:
                cc  = cluster[xc>=Rc[c,0] & xc<=Rc[c,1]]
            else:
                cc = None

            rdr_tmp = rdrobust(yc, xc,
                                fuzzy=fuzzy,
                                deriv=derivvec[c],
                                p=pvec[c],
                                q=qvec[c],
                                h=hmat[c],
                                b=bmat[c],
                                rho=rhovec[c],
                                covs=covs_aux,
                                covs_drop=covs_dropvec[c],
                                kernel=kernelvec[c],
                                weights=weightsc,
                                bwselect=bwselectvec[c],
                                scalepar=scaleparvec[c],
                                scaleregul=scaleregulvec[c],
                                masspoints=masspointsvec[c],
                                bwcheck=bwcheckvec[c],
                                bwrestrict=bwrestrictvec[c],
                                stdvars=stdvarsvec[c],
                                vce=vcevec[c],
                                nnmatch=nnmatchvec[c],
                                cluster=cc,
                                level=level)
            
            B[0,c] = rdr_tmp.Estimate.iloc[0,1]
            V[0,c] = rdr_tmp.se.iloc[2,0]**2
            Coefs[0,c] = rdr_tmp.Estimate.iloc[0,0]
            V_cl[0,c] = rdr_tmp.se.iloc[0,0]**2
            CI[:,c] = rdr_tmp.ci.iloc[2]
            CI_cl[:,c] = rdr_tmp.ci.iloc[0]
            H[:,c] = rdr_tmp.bws.iloc[0,:]
            Bbw[:,c] = rdr_tmp.bws.iloc[1,:]
            Nh[:,c] = rdr_tmp.N_h
            Pv[0,c] = rdr_tmp.pv.iloc[2,0]
            Pv_cl[0,c] = rdr_tmp.pv.iloc[0,0]

            c_disp.append(round(C[c], 2))
    else:

        for c in range(cnum):
            xc = np.sqrt((X - C[c])**2 + (X2 - C2[c])**2) * (2 * zvar - 1)

            yc = Y[(xc >= rangemat[c, 0]) & (xc <= rangemat[c, 1])]
            xc = xc[(xc >= rangemat[c, 0]) & (xc <= rangemat[c, 1])]

            weightsc = None
            if weightsvec is not None:
                weightaux = weightsvec[c]
                weightsc = eval(weightaux + '[xc >= rangemat[c, 0] & xc <= rangemat[c, 1]]')

            covs_aux = None
            if covs_mat is not None:
                covs_mat_c = covs_mat[(xc >= rangemat[c, 0]) & (xc <= rangemat[c, 1]), :]
                if covs_list is not None:
                    covs_aux = covs_mat_c[:, covs_list[c]]
                else:
                    covs_aux = covs_mat_c
            
            if cluster is not None:
                cc  = cluster[xc>=rangemat[c,0] & xc<=rangemat[c,1]]
            else:
                cc = None

            rdr_tmp = rdrobust(yc, xc,
                                fuzzy=fuzzy,
                                deriv=derivvec[c],
                                p=pvec[c],
                                q=qvec[c],
                                h=hmat[c],
                                b=bmat[c],
                                rho=rhovec[c],
                                covs=covs_aux,
                                covs_drop=covs_dropvec[c],
                                kernel=kernelvec[c],
                                weights=weightsc,
                                bwselect=bwselectvec[c],
                                scalepar=scaleparvec[c],
                                scaleregul=scaleregulvec[c],
                                masspoints=masspointsvec[c],
                                bwcheck=bwcheckvec[c],
                                bwrestrict=bwrestrictvec[c],
                                stdvars=stdvarsvec[c],
                                vce=vcevec[c],
                                nnmatch=nnmatchvec[c],
                                cluster=cc,
                                level=level)

            B[0,c] = rdr_tmp.Estimate.iloc[0,1]
            V[0,c] = rdr_tmp.se.iloc[2,0]**2
            Coefs[0,c] = rdr_tmp.Estimate.iloc[0,0]
            V_cl[0,c] = rdr_tmp.se.iloc[0,0]**2
            CI[:,c] = rdr_tmp.ci.iloc[2]
            CI_cl[:,c] = rdr_tmp.ci.iloc[0]
            H[:,c] = rdr_tmp.bws.iloc[0,:]
            Bbw[:,c] = rdr_tmp.bws.iloc[1,:]
            Nh[:,c] = rdr_tmp.N_h
            Pv[0,c] = rdr_tmp.pv.iloc[2,0]
            Pv_cl[0,c] = rdr_tmp.pv.iloc[0,0]

            c_disp.append('('+str(round(C[c], 2))+','+str(round(C2[c], 2))+')')

    #################################################################
    # Calculate pooled estimates
    #################################################################

    if xnorm is not None:

        if pooled_opt is None: pooled_opt =''
        else: pooled_opt = "," + pooled_opt
        rdr_tmp = eval(f'rdrobust(Y, xnorm, fuzzy=fuzzy{pooled_opt})')
        
        B[0,cnum] = rdr_tmp.Estimate.iloc[0,1]
        V[0,cnum] = rdr_tmp.se.iloc[2,0]**2
        Coefs[0,cnum] = rdr_tmp.Estimate.iloc[0,0]
        V_cl[0,cnum] = rdr_tmp.se.iloc[0,0]**2
        CI[:,cnum] = rdr_tmp.ci.iloc[2]
        CI_cl[:,cnum] = rdr_tmp.ci.iloc[0]
        H[:,cnum] = rdr_tmp.bws.iloc[0,:]
        Bbw[:,cnum] = rdr_tmp.bws.iloc[1,:]
        Nh[:,cnum] = rdr_tmp.N_h
        Pv[0,cnum] = rdr_tmp.pv.iloc[2,0]
        Pv_cl[0,cnum] = rdr_tmp.pv.iloc[0,0]

    #################################################################
    # Display results
    #################################################################

    print('')
    print('='*85)
    print('Cutoff'.ljust(16),
            'Coef.'.ljust(8),
            'P-value'.ljust(16),
            '95% CI'.ljust(16),
            'hl'.ljust(9),
            'hr'.ljust(9),
            'Nh'.ljust(5)), 
    print('='*85)
    
    if not conventional:
        for k in range(cnum):
            if C2 is None:
                print("{:4.2f}".format(c_disp[k]).ljust(17), end="")
            else:
                print('(', C[k], ' , ', C2[k],')'.ljust(4), sep="",end="")
            print(
                "{:7.3f}".format(Coefs[0,k]).ljust(9),
                "{:1.3f}".format(Pv[0,k]).ljust(9),
                "{:4.3f}".format(CI[0,k]).ljust(10),
                "{:4.3f}".format(CI[1,k]).ljust(10),
                "{:4.3f}".format(H[0,k]).ljust(9),
                "{:4.3f}".format(H[1,k]).ljust(9),
                "{:4.0f}".format(Nh[0,k]+Nh[1,k]).ljust(8))

        if xnorm is not None:
            print('-'*85)

            print('Pooled'.ljust(15),
                "{:7.3f}".format(Coefs[0,cnum]).ljust(9),
                "{:1.3f}".format(Pv[0,cnum]).ljust(9),
                "{:4.3f}".format(CI[0,cnum]).ljust(10),
                "{:4.3f}".format(CI[1,cnum]).ljust(10),
                "{:4.3f}".format(H[0,cnum]).ljust(9),
                "{:4.3f}".format(H[1,cnum]).ljust(9),
                "{:4.0f}".format(Nh[0,cnum]+Nh[1,cnum]).ljust(8))
    else:
        for k in range(cnum):
            if C2 is None:
                print("{:4.2f}".format(c_disp[k]).ljust(17), end="")
            else:
                print('(', C[k], ' , ', C2[k],')'.ljust(4), sep="",end="")
            print(
                "{:7.3f}".format(Coefs[0,k]).ljust(9),
                "{:1.3f}".format(Pv_cl[0,k]).ljust(9),
                "{:4.3f}".format(CI_cl[0,k]).ljust(10),
                "{:4.3f}".format(CI_cl[1,k]).ljust(10),
                "{:4.3f}".format(H[0,k]).ljust(9),
                "{:4.3f}".format(H[1,k]).ljust(9),
                "{:4.0f}".format(Nh[0,k]+Nh[1,k]).ljust(8))
        
        if xnorm is not None:
            print('-'*85)
            print('Pooled'.ljust(15),
            "{:7.3f}".format(Coefs[0,cnum]).ljust(9),
            "{:1.3f}".format(Pv_cl[0,cnum]).ljust(9),
            "{:4.3f}".format(CI_cl[0,cnum]).ljust(10),
            "{:4.3f}".format(CI_cl[1,cnum]).ljust(10),
            "{:4.3f}".format(H[0,cnum]).ljust(9),
            "{:4.3f}".format(H[1,cnum]).ljust(9),
            "{:4.0f}".format(Nh[0,cnum]+Nh[1,cnum]).ljust(8))

    print('='*85)

    #################################################################
    # Return values
    ################################################################

    colnames = list((np.arange(cnum)+1).astype('str')) + ["pooled"]
    rownames = pd.Index(["left","right"])

    output = {
        "B": pd.DataFrame(B,columns = colnames),
        "V": pd.DataFrame(V, columns = colnames),
        "Coefs": pd.DataFrame(Coefs, columns = colnames),
        "V_cl": pd.DataFrame(V_cl, columns = colnames),
        "Nh": pd.DataFrame(Nh, columns = colnames, index = rownames),
        "CI": pd.DataFrame(CI, columns = colnames),
        "CI_cl": pd.DataFrame(CI_cl, columns = colnames),
        "H": pd.DataFrame(H, columns = colnames, index = rownames),
        "Bbw": pd.DataFrame(Bbw, columns = colnames, index = rownames),
        "Pv": pd.DataFrame(Pv, columns = colnames),
        "Pv_cl": pd.DataFrame(Pv_cl, columns = colnames)
    }

    return output