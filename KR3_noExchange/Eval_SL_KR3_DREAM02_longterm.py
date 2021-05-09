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
import Initialize_KR3_DREAM_LT as iKR
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
          = mymod.likelihood_Emission_Recirculate_LF (par,iKR.tdseries, iKR.lF, iKR.meteo_data,
                                 iKR.infilF, iKR.meas_data_flow, iKR.lab_data,
                                 iKR.sdData)

    #return logp
    return logp, emOut, lFPar,emLF


d_iter = 1000
nchains = 5
d_max = 63000
npar = iKR.par_df.shape[1]
parfn = './Results/KR3_SL_'
logpsfn = './Results/KR3_SL_logps_'

sampled_params = np.zeros([nchains,d_max,npar])
log_ps = np.zeros([nchains,d_max])


# Loop through all files
for chain in np.arange(nchains):
    for nsamp in np.arange(d_iter,d_max+d_iter,d_iter):
        fn = parfn + str(chain) + '_' + str(nsamp) + '.npy'
        fnlogps = logpsfn + str(chain) + '_' + str(nsamp) + '.npy'
        sampled_params[chain][nsamp-d_iter:nsamp] = np.load(fn, allow_pickle=False)
        log_ps[chain][nsamp-d_iter:nsamp] = np.load(fnlogps).squeeze()
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]

maxidx = np.unravel_index(log_ps.argmax(), log_ps.shape)
par = sampled_params[maxidx]
GR = Gelman_Rubin(sampled_params)


# save starts for possible restart of dream optimization
#save starts to pickle...
fo = open('starts01.pkl','wb')
pkl.dump(starts,fo)
fo.close()



# fo = open('starts01.pkl','rb')
# starts = pkl.load(fo)
# par = starts[0]

# par = iKR.par_df.iloc[2].values
# starts = [par]

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
drainVolN_tot = []
emLF_tot = []

par_set = starts
par_setdf = pd.DataFrame(par_set,columns=iKR.par_df.columns)

for ii in np.arange(len(par_set)):
    # Run model for parameters in chain
    logp, emOut, lFPar, emLF \
        = likelihood(par_set[ii])

    cumqDr, cDr, qDr, qInf, thCL, mCL, cInfCL, mbalerrCL, qrunoffCL, vW, mW, drainLev, \
        drainOutF, drainVolN = emOut

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
    drainLev_tot.append(drainLev)
    drainOutF_tot.append(drainOutF)
    drainVolN_tot.append(drainVolN)
    lFPar_tot.append(lFPar)
    emLF_tot.append(emLF)

    totF_sim.append(cumqDr_tot[ii][iKR.tmeas_ind] - cumqDr_tot[ii][iKR.tmeas_ind[0]])
    drainLev_sim.append(drainLev[iKR.tmeas_ind])
    
#import make_plotsBB
#cLOut = [qinf, qrf, qev, seff, th, mbalerr, cInf, m]
