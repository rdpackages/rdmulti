#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas  as pd
from scipy.stats import norm
import matplotlib.pyplot as plt
import warnings
from rdrobust import rdrobust

def rdmc(Y,X,C,fuzzy=None,derivvec=None,pooled_opt=None,verbose=False,
            pvec=None,qvec=None,hmat=None,bmat=None,rhovec=None,
            covs_mat=None,covs_list=None,covs_dropvec=None,kernelvec=None,weightsvec=None,
            bwselectvec=None,scaleparvec=None,scaleregulvec=None,
            masspointsvec=None,bwcheckvec=None,bwrestrictvec=None,
            stdvarsvec=None,vcevec=None,nnmatchvec=None,cluster=None,
            level=95,plot=False,conventional=False):
    
    """
    Analysis of RD designs with multiple cutoffs.

    Arguments:
    - Y: outcome variable.
    - X: running variable.
    - C: cutoff variable.
    - fuzzy: specifies a fuzzy design.
    - derivvec: vector of cutoff-specific order of derivatives.
    - pooled_opt: options to be passed to rdrobust() to calculate pooled estimand.
    - verbose: displays the output from rdrobust for estimating the pooled estimand.
    - pvec: vector of cutoff-specific polynomial orders.
    - qvec: vector of cutoff-specific polynomial orders for bias estimation.
    - hmat: matrix of cutoff-specific bandwidths.
    - bmat: matrix of cutoff-specific bandwidths for bias estimation.
    - rhovec: vector of cutoff-specific values of rho.
    - covs_mat: matrix of covariates.
    - covs_list: list of covariates to be used in each cutoff.
    - covs_dropvec: vector indicating whether collinear covariates should be dropped at each cutoff.
    - kernelvec: vector of cutoff-specific kernels.
    - weightsvec: vector of cutoff-specific weights.
    - bwselectvec: vector of cutoff-specific bandwidth selection methods.
    - scaleparvec: vector of cutoff-specific scale parameters.
    - scaleregulvec: vector of cutoff-specific scale regularization parameters.
    - masspointsvec: vector indicating how to handle repeated values at each cutoff.
    - bwcheckvec: vector indicating the value of bwcheck at each cutoff.
    - bwrestrictvec: vector indicating whether computed bandwidths are restricted to the range or runvar at each cutoff.
    - stdvarsvec: vector indicating whether variables are standardized at each cutoff.
    - vcevec: vector of cutoff-specific variance-covariance estimation methods.
    - nnmatchvec: vector of cutoff-specific nearest neighbors for variance estimation.
    - cluster: cluster ID variable.
    - level: confidence level for confidence intervals.
    - plot: plots cutoff-specific estimates and weights.
    - conventional: reports conventional, instead of robust-bias corrected, p-values and confidence intervals.

    Returns:
    - tau: pooled estimate.
    - se_rb: robust bias-corrected standard error for pooled estimate.
    - pv_rb: robust bias-corrected p-value for pooled estimate.
    - ci_rb_l: left limit of robust bias-corrected CI for pooled estimate.
    - ci_rb_r: right limit of robust bias-corrected CI for pooled estimate.
    - hl: bandwidth to the left of the cutoff for pooled estimate.
    - hr: bandwidth to the right of the cutoff for pooled estimate.
    - Nhl: sample size within bandwidth to the left of the cutoff for pooled estimate.
    - Nhr: sample size within bandwidth to the right of the cutoff for pooled estimate.
    - B: vector of bias-corrected estimates.
    - V: vector of robust variances of the estimates.
    - Coefs: vector of conventional estimates.
    - W: vector of weights for each cutoff-specific estimate.
    - Nh: vector of sample sizes within bandwidth.
    - CI: robust bias-corrected confidence intervals.
    - H: matrix of bandwidths.
    - Pv: vector of robust p-values.
    - rdrobust_results: results from rdrobust for pooled estimate.
    - cfail: Cutoffs where rdrobust() encountered problems.

    See Also
    --------
    rdms, rdmcplot
    
    Example
    ------- 
    X = numpy.random.uniform(0, 100, size=1000)
    C = numpy.repeat([33, 66], repeats=[500, 500])
    Y = (1 + X + (X >= C)) * (C == 33) + (.5 + .5 * X + .8 * (X >= C)) * (C == 66) + numpy.random.standard_normal(1000)
    rdmc(Y, X, C)
    """

    #################################################################
    # Setup and error checking
    #################################################################

    try: C = np.array(C, dtype = 'float64')
    except: raise Exception('C has to be numeric')
    if (np.nanmax(C)>=np.nanmax(X)) or (np.nanmin(C)<=np.nanmin(X)):
        raise Exception('Cutoff variable outside range of running variable')

    clist = np.sort(np.unique(C))
    cnum = len(clist)

    Xc = X - C

    if derivvec is None: derivvec = np.repeat(None,cnum)
    if pvec is None: pvec = np.repeat(None,cnum)
    if qvec is None: qvec = np.repeat(None,cnum) 
    if rhovec is None: rhovec = np.repeat(None,cnum)
    if bwcheckvec is None: bwcheckvec = np.repeat(None,cnum)


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

    if weightsvec is not None:
        weightsvec = weightsvec.reshape(len(Y),-1)

    if covs_dropvec is None: covs_dropvec = np.repeat(True,cnum)
    if kernelvec is None: kernelvec = np.repeat('tri',cnum)
    if bwselectvec is None: bwselectvec = np.repeat('mserd',cnum)
    if vcevec is None: vcevec = np.repeat('nn',cnum)
    if nnmatchvec is None: nnmatchvec = np.repeat(3,cnum)
    if scaleparvec is None: scaleparvec = np.repeat(1,cnum)
    if scaleregulvec is None: scaleregulvec = np.repeat(1,cnum)
    if masspointsvec is None: masspointsvec = np.repeat('adjust',cnum)
    if bwrestrictvec is None: bwrestrictvec = np.repeat(True,cnum)
    if stdvarsvec is None: stdvarsvec = np.repeat(False,cnum)
    if covs_mat is not None:
        covs_mat = np.array(covs_mat).reshape(len(covs_mat),-1)
        if (covs_list is not None) and (len(covs_list)!=cnum):
            raise Exception('Elements in covs_list should equal number of cutoffs')
  
    B = np.full((1,cnum+2),np.nan) 
    V = np.full((1,cnum+2),np.nan) 
    V_cl = np.full((1,cnum+2),np.nan) 
    Coefs = np.full((1,cnum+2),np.nan)
    Nh = np.full((2,cnum+2),np.nan)
    CI = np.full((2,cnum+2),np.nan)
    CI_cl = np.full((2,cnum+2),np.nan)
    Pv = np.full((1,cnum+2),np.nan)
    Pv_cl = np.full((1,cnum+2),np.nan)
    H = np.full((2,cnum+2),np.nan)
    Bbw = np.full((2,cnum+2),np.nan)
    W = np.full((1,cnum),np.nan)
    Cfail = []

    #################################################################
    # Calculate pooled estimate
    #################################################################

    if pooled_opt is None: pooled_opt =''
    else: pooled_opt = "," + pooled_opt
    rdr = eval(f'rdrobust(Y,Xc,fuzzy=fuzzy{pooled_opt})')

    hl, hr = rdr.bws.iloc[0].values
    Nhl, Nhr = rdr.N_h
    B[0,cnum+1] = rdr.Estimate.iloc[0,1]
    V[0,cnum+1] = rdr.se.iloc[2,0]**2
    Coefs[0,cnum+1] = rdr.Estimate.iloc[0,0]
    V_cl[0,cnum+1] = rdr.se.iloc[0,0]**2
    CI[:,cnum+1] = rdr.ci.iloc[2]
    CI_cl[:,cnum+1] = rdr.ci.iloc[0]
    H[:,cnum+1] = rdr.bws.iloc[0,:]
    Bbw[:,cnum+1] = rdr.bws.iloc[1,:]
    Nh[:,cnum+1] = rdr.N_h
    Pv[0,cnum+1] = rdr.pv.iloc[2,0]
    Pv_cl[0,cnum+1] = rdr.pv.iloc[0,0]

    #################################################################
    # Calculate cutoff-specific estimates and weights
    #################################################################

    count  = count_fail = 0
    epsilon = np.finfo(np.float64).eps
    for c in clist:

        yc = Y[np.abs(C-c)<=epsilon]
        xc = Xc[np.abs(C-c)<=epsilon]

        if weightsvec is not None:
            weightsc = weightsvec[np.abs(C-c)<=epsilon, count] 
        else:
            weightsc = None

        if covs_mat is not None:
            covs_mat_c = covs_mat[np.abs(C-c)<=epsilon,:]
            if covs_list is not None:
                covs_aux = covs_mat_c[:,covs_list[count]]
            else:
                covs_aux = covs_mat_c
        else:
            covs_aux = None

        if cluster is not None:
            cc  = cluster[np.abs(C-c)<=epsilon]
        else:
            cc = None
        try:
            rdr_tmp = rdrobust(yc,xc,
                                fuzzy=fuzzy,
                                deriv=derivvec[count],
                                p=pvec[count],
                                q=qvec[count],
                                h=hmat[count],
                                b=bmat[count],
                                rho=rhovec[count],
                                covs=covs_aux,
                                covs_drop=covs_dropvec[count],
                                kernel=kernelvec[count],
                                weights=weightsc,
                                bwselect=bwselectvec[count],
                                scalepar=scaleparvec[count],
                                scaleregul=scaleregulvec[count],
                                masspoints=masspointsvec[count],
                                bwcheck=bwcheckvec[count],
                                bwrestrict=bwrestrictvec[count],
                                stdvars=stdvarsvec[count],
                                vce=vcevec[count],
                                nnmatch=nnmatchvec[count],
                                cluster=cc,
                                level=level)
        except:
            Cfail = np.append(Cfail,c) 
            count_fail +=1
        else:
            B[0,count] = rdr_tmp.Estimate.iloc[0,1]
            V[0,count] = rdr_tmp.se.iloc[2,0]**2
            Coefs[0,count] = rdr_tmp.Estimate.iloc[0,0]
            V_cl[0,count] = rdr_tmp.se.iloc[0,0]**2
            CI[:,count] = rdr_tmp.ci.iloc[2]
            CI_cl[:,count] = rdr_tmp.ci.iloc[0]
            H[:,count] = rdr_tmp.bws.iloc[0,:]
            Bbw[:,count] = rdr_tmp.bws.iloc[1,:]
            Nh[:,count] = rdr_tmp.N_h
            Pv[0,count] = rdr_tmp.pv.iloc[2,0]
            Pv_cl[0,count] = rdr_tmp.pv.iloc[0,0]
        finally:
            count += 1


    # Weights

    W[0,:] = np.nansum(Nh[:,:cnum], axis = 0)/np.nansum(Nh[:,:cnum])
    W[np.isnan(W)] = 0

    # Weighted estimate

    Baux = B.copy()
    Vaux = V.copy()
    Coefsaux = Coefs.copy()
    Vaux_cl = V_cl.copy()
    Baux[np.isnan(Baux)] = 0
    Vaux[np.isnan(Vaux)] = 0
    Coefsaux[np.isnan(Coefsaux)] = 0
    Vaux_cl[np.isnan(Vaux_cl)] = 0

    B[0,cnum] = np.matmul(Baux[0,:cnum], W.T)
    V[0,cnum] = np.matmul(Vaux[0,:cnum], (W**2).T)
    Coefs[0,cnum] = np.matmul(Coefsaux[0,:cnum],  W.T)
    V_cl[0,cnum] = np.matmul(Vaux_cl[0,:cnum], (W**2).T)

    Nh[:,cnum] = np.nansum(Nh[:,:cnum],axis=1)

    CI[:,cnum] = B[0,cnum] + np.sqrt(V[0,cnum])*norm.ppf(1-(1-level/100)/2)*np.array([-1,1])
    Pv[0,cnum] = 2*(1-norm.cdf(np.abs(B[0,cnum]/np.sqrt(V[0,cnum]))))
    CI_cl[:,cnum] = Coefs[0,cnum] + np.sqrt(V_cl[0,cnum])*norm.ppf(1-(1-level/100)/2)*np.array([-1,1]) 
    Pv_cl[0,cnum] = 2*(1-norm.cdf(np.abs(Coefs[0,cnum]/np.sqrt(V_cl[0,cnum]))))


    #################################################################
    # Display results
    #################################################################

    if verbose==True: print(rdr)
    
    if conventional==False:
        print('')
        print('Cutoff-specific RD estimation with robust bias-corrected inference')
        print('='*90)
        print('Cutoff'.ljust(11),
              'Coef.'.ljust(8),
              'P-value'.ljust(16),
              '95% CI'.ljust(16),
              'hl'.ljust(9),
              'hr'.ljust(9),
              'Nh'.ljust(5), 
              'Weight'.ljust(5))
        print('='*90)

        for k in range(cnum):
            print("{:4.3f}".format(clist[k]).ljust(11),
                    "{:7.3f}".format(Coefs[0,k]).ljust(9),
                    "{:1.3f}".format(Pv[0,k]).ljust(9),
                    "{:4.3f}".format(CI[0,k]).ljust(10),
                    "{:4.3f}".format(CI[1,k]).ljust(10),
                    "{:4.3f}".format(H[0,k]).ljust(9),
                    "{:4.3f}".format(H[1,k]).ljust(9),
                    "{:4.0f}".format(Nh[0,k]+Nh[1,k]).ljust(8),
                    "{:1.3f}".format(W[0,k]).ljust(5))
            
        print('='*90)
        print('Weighted'.ljust(11),
                "{:7.3f}".format(Coefs[0,cnum]).ljust(9),
                "{:1.3f}".format(Pv[0,cnum]).ljust(9),
                "{:4.3f}".format(CI[0,cnum]).ljust(10),
                "{:4.3f}".format(CI[1,cnum]).ljust(10),
                '  .'.ljust(9),
                '  .'.ljust(9),
                "{:4.0f}".format(Nh[0,cnum]+Nh[1,cnum]).ljust(8),
                '  .'.ljust(5))
        
        print('Pooled'.ljust(11),
                "{:7.3f}".format(Coefs[0,cnum+1]).ljust(9),
                "{:1.3f}".format(Pv[0,cnum+1]).ljust(9),
                "{:4.3f}".format(CI[0,cnum+1]).ljust(10),
                "{:4.3f}".format(CI[1,cnum+1]).ljust(10),
                "{:4.3f}".format(H[0,cnum+1]).ljust(9),
                "{:4.3f}".format(H[1,cnum+1]).ljust(9),
                "{:4.0f}".format(Nh[0,cnum+1]+Nh[1,cnum+1]).ljust(8),
                '  .'.ljust(5))
        print('='*90)

    else:
        print('')
        print('Cutoff-specific RD estimation with conventional inference')
        print('='*90)
        print('Cutoff'.ljust(11),
              'Coef.'.ljust(8),
              'P-value'.ljust(16),
              '95% CI'.ljust(16),
              'hl'.ljust(9),
              'hr'.ljust(9),
              'Nh'.ljust(5), 
              'Weight'.ljust(5))
        print('='*90)

        for k in range(cnum):
            print("{:4.3f}".format(clist[k]).ljust(11),
                    "{:7.3f}".format(Coefs[0,k]).ljust(9),
                    "{:1.3f}".format(Pv_cl[0,k]).ljust(9),
                    "{:4.3f}".format(CI_cl[0,k]).ljust(10),
                    "{:4.3f}".format(CI_cl[1,k]).ljust(10),
                    "{:4.3f}".format(H[0,k]).ljust(9),
                    "{:4.3f}".format(H[1,k]).ljust(9),
                    "{:4.0f}".format(Nh[0,k]+Nh[1,k]).ljust(8),
                    "{:1.3f}".format(W[0,k]).ljust(5))
            
        print('='*90)
        print('Weighted'.ljust(11),
                "{:7.3f}".format(Coefs[0,cnum]).ljust(9),
                "{:1.3f}".format(Pv_cl[0,cnum]).ljust(9),
                "{:4.3f}".format(CI_cl[0,cnum]).ljust(10),
                "{:4.3f}".format(CI_cl[1,cnum]).ljust(10),
                '  .'.ljust(9),
                '  .'.ljust(9),
                "{:4.0f}".format(Nh[0,cnum]+Nh[1,cnum]).ljust(8),
                '  .'.ljust(5))
        
        print('Pooled'.ljust(11),
                "{:7.3f}".format(Coefs[0,cnum+1]).ljust(9),
                "{:1.3f}".format(Pv_cl[0,cnum+1]).ljust(9),
                "{:4.3f}".format(CI_cl[0,cnum+1]).ljust(10),
                "{:4.3f}".format(CI_cl[1,cnum+1]).ljust(10),
                "{:4.3f}".format(H[0,cnum+1]).ljust(9),
                "{:4.3f}".format(H[1,cnum+1]).ljust(9),
                "{:4.0f}".format(Nh[0,cnum+1]+Nh[1,cnum+1]).ljust(8),
                '  .'.ljust(5))
        print('='*90)

    if count_fail>0:
        warnings.warn("rdrobust() could not run in one or more cutoffs.")

    #################################################################
    # Plots
    #################################################################

    if plot==True:

        ylim = 1.3*np.array([np.nanmin(CI), np.nanmax(CI)])
        xlim= np.array([np.nanmin(clist), np.nanmax(clist)])

        plt.subplot(1, 2, 1)

        for k in range(cnum):
            if k==0:
                plt.plot(clist[k], Coefs[0,k],'o', c='blue',label ='Estimate')
                if conventional==False:
                    plt.plot(np.repeat(clist[k],2), CI[:,k], c = 'blue', label = str(level)+'% CI')
                else:
                    plt.plot(np.repeat(clist[k],2), CI_cl[:,k], c = 'blue', label = str(level)+'% CI')
            else:
                plt.plot(clist[k], Coefs[0,k],'o', c='blue')
                if conventional==False:
                    plt.plot(np.repeat(clist[k],2), CI[:,k], c = 'blue')
                else:
                    plt.plot(np.repeat(clist[k],2), CI_cl[:,k], c = 'blue')

        plt.axhline(y = Coefs[0,cnum+1], color = 'gray', label ='Pooled estimate')
        plt.axhline(y = Coefs[0,cnum], color = 'red', label ='Weighted estimate')
        plt.axhline(y = 0, color = 'black',linestyle = 'dotted')
        plt.fill_between(xlim, CI[0,cnum], CI[1,cnum], color='red', alpha=.1)
        plt.fill_between(xlim, CI[0,cnum+1], CI[1,cnum+1], color='gray', alpha=.2)
        plt.xlabel("Cutoff")
        plt.ylabel("Treatment Effect")
        plt.legend()
        plt.subplots_adjust(wspace=0.4)

        plt.subplot(1, 2, 2)
        plt.bar(np.array(clist, dtype = str),height = W[0], color = 'gray')
        plt.xlabel("Cutoff")
        plt.ylabel("Weight")
        plt.show()
            

    #################################################################
    # Return values
    #################################################################

    class rdmc_output:
        def __init__(self, tau, se_rb, pv_rb, ci_rb_l, ci_rb_r, hl, hr, Nhl, Nhr, B, V, Coefs,
                     V_cl, W, Nh, CI, CI_cl, H, Bbw, Pv, Pv_cl, rdrobust_results, cfail):

            self.tau = tau,
            self.se_rb = se_rb,
            self.pv_rb = pv_rb,
            self.ci_rb_l = ci_rb_l,
            self.ci_rb_r = ci_rb_r,
            self.hl = hl,
            self.hr = hr,
            self.Nhl = Nhl,
            self.Nhr = Nhr,
            self.B = B,
            self.V = V,
            self.Coefs = Coefs,
            self.V_cl = V_cl,
            self.W = W,
            self.Nh = Nh,
            self.CI = CI,
            self.CI_cl = CI_cl,
            self.H = H,
            self.Bbw = Bbw,
            self.Pv = Pv,
            self.Pv_cl = Pv_cl,
            self.rdrobust_results = rdrobust_results,
            self.cfail = cfail

    colnames = list((np.arange(cnum)+1).astype('str')) + ["weighted","pooled"]
    rownames = pd.Index(["left","right"])

    output = rdmc_output(tau = rdr.Estimate['tau.us'],
                        se_rb = rdr.se.iloc[2],
                        pv_rb = rdr.pv.iloc[2],
                        ci_rb_l = rdr.ci.iloc[2,0],
                        ci_rb_r = rdr.ci.iloc[2,1],
                        hl = rdr.bws.iloc[0,0],
                        hr = rdr.bws.iloc[0,1],
                        Nhl = rdr.N_h[0],
                        Nhr = rdr.N_h[1],
                        B = pd.DataFrame(B, columns = colnames),
                        V = pd.DataFrame(V, columns = colnames),
                        Coefs = pd.DataFrame(Coefs, columns = colnames),
                        V_cl = pd.DataFrame(V_cl, columns = colnames),
                        W = W,
                        Nh = pd.DataFrame(Nh, columns = colnames, index = rownames),
                        CI = pd.DataFrame(CI, columns = colnames),
                        CI_cl = pd.DataFrame(CI_cl, columns = colnames),
                        H = pd.DataFrame(H, columns = colnames, index = rownames),
                        Bbw = pd.DataFrame(Bbw, columns = colnames, index = rownames),
                        Pv = pd.DataFrame(Pv, columns = colnames),
                        Pv_cl = pd.DataFrame(Pv_cl, columns = colnames),
                        rdrobust_results = rdr,
                        cfail = Cfail)
    
    return output