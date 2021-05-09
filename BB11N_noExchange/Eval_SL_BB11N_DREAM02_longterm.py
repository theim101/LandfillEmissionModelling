"""
Created on Sat Oct 20 2018
Modules for extracting data from the Chronos database

@author: T.J. Heimovaara
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from collections import namedtuple
import seaborn as sns
import Initialize_BB11N_DREAM_LT as iBB
#import Initialize_WM_DREAM01 as iBB
import scipy.stats as scpstats
from context import mymod
import pickle as pkl
from pydream.convergence import Gelman_Rubin
#import inspect

# set seaborn graphics style
sns.set()

def likelihood(par):
    # Run Forward model using parameters in par
    logp, emOut, lFPar, emLF \
        = mymod.likelihood_Emission_LF(par, iBB.tdseries,iBB.lF, 
                                       iBB.meteo_data,  
                                       iBB.meas_data_flow, iBB.lab_data, 
                                       iBB.sdData)

    return logp, emOut, lFPar, emLF


d_iter = 5000
nchains = 5
d_max = 35000
npar = iBB.par_df.shape[1]
parfn = './Results/BB11N_SL_'
logpsfn = './Results/BB11N_SL_logps_'

sampled_params = np.zeros([nchains,d_max,npar])
log_ps = np.zeros([nchains,d_max])

#GR = Gelman_Rubin(sampled_params)

#Loop through all files
for chain in np.arange(nchains):
    for nsamp in np.arange(d_iter,d_max+d_iter,d_iter):
        fn = parfn + str(chain) + '_' + str(nsamp) + '.npy'
        fnlogps = logpsfn + str(chain) + '_' + str(nsamp) + '.npy'
        sampled_params[chain][nsamp-d_iter:nsamp] = np.load(fn)
        log_ps[chain][nsamp-d_iter:nsamp] = np.load(fnlogps).squeeze()
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]

maxidx = np.unravel_index(log_ps.argmax(), log_ps.shape)
par = sampled_params[maxidx]


# save starts for possible restart of dream optimization
# save starts to pickle...
fo = open('starts01.pkl','wb')
pkl.dump(starts,fo)
fo.close()



# fo = open('starts01.pkl','rb')
# starts = pkl.load(fo)
# par = starts[0]
#par = iBB.par_df.iloc[2].values

logp_tot = []
cumqDr_tot = []
cDr_tot = []
qDr_tot = []
qInf_tot = []

thCL_tot = []
mCL_tot = []
cInf_CL = []
mbalerr_CL = []
qrunoff_CL = []

vW_tot = []
mW_tot = []
totF_sim = []
lFPar_tot = []
drainLev_tot = []
drainLev_sim = []
drainOutF_tot = []
emLF_tot = []

par_set = starts

par_setdf = pd.DataFrame(par_set,columns=iBB.par_df.columns)
for ii in np.arange(len(par_set)):
    
    logp, emOut, lFPar, emLF \
        = likelihood(par_set[ii])

    cumqDr, cDr, qDr, qInf, thCL, mCL, cInfCL, mbalerrCL, qrunoffCL, vW, mW \
        = emOut

    logp_tot.append(logp)
    cumqDr_tot.append(cumqDr)
    cDr_tot.append(cDr)
    qDr_tot.append(qDr)
    qInf_tot.append(qInf)
    
    thCL_tot.append(thCL)
    mCL_tot.append(mCL)
    cInf_CL.append(cInfCL)
    mbalerr_CL.append(mbalerrCL)
    qrunoff_CL.append(qrunoffCL)
    vW_tot.append(vW)
    mW_tot.append(mW)
    lFPar_tot.append(lFPar)
    emLF_tot.append(emLF)
    
    totF_sim.append(cumqDr_tot[ii][iBB.tmeas_ind] - cumqDr_tot[ii][iBB.tmeas_ind[0]])



