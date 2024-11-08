###################################################################
# rdmulti: analysis of RD designs with multiple cutoffs or scores
# Illustration file
# Authors: Matias Cattaneo, Ricardo Masini, Rocio Titiunik, 
#          Gonzalo Vazquez-Bare
###################################################################

import numpy as np
import pandas as pd
from rdmulti import rdmc, rdms, rdmcplot

########## To run from the source ##################################################
# (in terminal) pip uninstall binsreg
# import sys
# import numpy as np
# sys.path.insert(0, '/Dropbox/rdmulti/python/rdmulti/scr/rdmulti')
# from rdmc import rdmc
# from rdms import rdms
# from rdmcplot import rdmcplot
###################################################################################

data = pd.read_csv('simdata_multic.csv')
Y = data.y
X = data.x
C = data.c

###################################################################
# Multiple cutoffs
###################################################################

aux = rdmc(Y,X,C)

aux = rdmc(Y,X,C, plot = True)

aux = rdmc(Y, X, C, pooled_opt='h=20,p=2', verbose=True)

aux = rdmc(Y, X, C, hmat=[11, 10])

aux = rdmc(Y, X, C, bwselectvec=['msetwo', 'certwo'])

# Add four covariates
Z = np.random.randn(len(Y), 4)

# Including all covariates in each cutoff
aux = rdmc(Y, X, C, covs_mat=Z)

# Use covariates Z1 and Z2 in cutoff 1, all four covariates in cutoff 2
covlist = [[0, 1], list(range(4))]
aux = rdmc(Y, X, C, covs_mat=Z, covs_list=covlist)

# Add weights
wvar = np.random.uniform(0.8, 1.2, len(Y)*len(np.unique(C))).reshape(len(Y),-1)
aux = rdmc(Y, X, C, weightsvec=wvar)

###################################################################
# Plots
###################################################################

aux = rdmcplot(Y, X, C)

aux = rdmcplot(Y, X, C, nobins=True)

aux = rdmcplot(Y, X, C, hmat=[11, 12], pvec=[1, 1])

aux = rdmcplot(Y, X, C, covs_mat=Z, covs_list=covlist)

##################################################################
# Cumulative cutoffs
##################################################################

data = pd.read_csv('simdata_cumul.csv')
Y = data['y']
X = data['x']
cvec = np.array([data['c'][0], data['c'][1]])

aux = rdms(Y, X, cvec)

aux = rdms(Y, X, cvec, hmat=(11, 8), kernelvec=('uniform', 'triangular'))

rangemat = np.array([[0, 35.0], [60.5, 100]])
aux = rdms(Y, X, cvec, rangemat=rangemat)

cutoff = cvec[0] * (X <= 49.5) + cvec[1] * (X > 49.5)
aux = rdmc(Y, X, cutoff)

aux = rdmcplot(Y, X, cutoff)

# Add four covariates
Z = np.random.randn(len(Y), 4)

# Including all covariates in each cutoff
aux = rdms(Y, X, cvec, covs_mat=Z)

# Use covariates Z1 and Z2 in cutoff 1, all four covariates in cutoff 2
covlist = [[0, 1], list(range(4))]
aux = rdms(Y, X, cvec, covs_mat=Z, covs_list=covlist)

# Bivariate score
data = pd.read_csv('simdata_multis.csv')
Y = data['y']
X1 = data['x1']
X2 = data['x2']
zvar = data['t']
cvec = [data['c1'][0], data['c1'][1], data['c1'][2]]
cvec2 = [data['c2'][0], data['c2'][1], data['c2'][2]]

aux = rdms(Y, X1, cvec, X2, zvar, cvec2)

aux = rdms(Y, X1, cvec, X2, zvar, cvec2, hmat=[15, 13, 17])

xnorm = np.min(np.column_stack((np.abs(50 - X1), np.abs(50 - X2))), axis=1) * (2 * zvar - 1)
aux = rdms(Y, X1, cvec, X2, zvar, cvec2, xnorm=xnorm)

# Add four covariates
Z = np.random.randn(len(Y), 4)

# Including all covariates in each cutoff
aux = rdms(Y, X1, cvec, X2, zvar, cvec2, covs_mat=Z)

# Use covariates Z1 and Z2 in cutoff 1, all four covariates in cutoff 2, covariates Z1 and Z3 in cutoff 3
covlist = [[0, 1], list(range(4)), [0, 2]]
aux = rdms(Y, X1, cvec, X2, zvar, cvec2, covs_mat=Z, covs_list=covlist)