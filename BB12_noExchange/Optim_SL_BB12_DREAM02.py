"""
Created on Sat Oct 20 2018
Modules for extracting data from the Chronos database

@author: T.J. Heimovaara
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import Initialize_BB12_DREAM01 as iBB
from context import mymod 
import scipy.stats as stats

from pydream.core import run_dream
from pydream.convergence import Gelman_Rubin
from pydream.parameters import SampledParam

import pickle as pkl

sns.set()

def likelihood(par):
    # Run Forward model using parameters in par

    logp, emOut, lFPar, emLF \
          = mymod.likelihood_Emission_LF (par,iBB.tdseries, iBB.lF, iBB.meteo_data,
                                 iBB.meas_data_flow, iBB.lab_data,
                                 iBB.sdData)

    return logp
    #return logp, emOut, lFPar,emLF
    
    
parfn = './Results/BB12_SL_'
logpsfn = './Results/BB12_SL_logps_'

niterations = 5000
converged = False
nchains = 5
d_iter = niterations
d_max = 25000 # only relevant if restart is true... (last file to be read)
restart = False

npar = iBB.par_df.shape[1]
lower_limits = iBB.par_df.iloc[0].values
range_vals = iBB.par_df.iloc[1].values - iBB.par_df.iloc[0].values

parameters_to_sample = SampledParam(stats.uniform, loc=lower_limits, scale=range_vals)
sampled_parameter_names = [parameters_to_sample]


if restart:
    sampled_params = np.zeros([nchains,d_max,npar])
    log_ps = np.zeros([nchains,d_max])
    start_random = False

    #Loop through all files
    for chain in np.arange(nchains):
        for nsamp in np.arange(d_iter,d_max+d_iter,d_iter):
            fn = parfn+str(chain) + '_' + str(nsamp) + '.npy'
            fnlogps = logpsfn + str(chain) + '_' + str(nsamp) + '.npy'
            sampled_params[chain][nsamp-d_iter:nsamp] = np.load(fn)
            log_ps[chain][nsamp-d_iter:nsamp] = np.load(fnlogps).squeeze()

    old_samples = sampled_params
    starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
    total_iterations = d_max + niterations
else:
    # start_random = False
    # fo = open('starts01.pkl','rb')
    # starts = pkl.load(fo)

    start_random = True
    starts = None
    total_iterations = niterations


sampled_params, log_ps = run_dream(sampled_parameter_names,
                                   likelihood=likelihood,
                                   niterations=niterations, nchains=nchains,
                                   restart=restart, start=starts, start_random=start_random,
                                   save_history=True, history_thin=10,
                                   adapt_gamma=True, adapt_crossover=True,
                                   model_name=parfn,
                                   verbose=True)

# Save sampling output (sampled parameter values and their corresponding logps).

for chain in range(len(sampled_params)):
    np.save(parfn + str(chain)+'_'+str(total_iterations), sampled_params[chain])
    np.save(logpsfn + str(chain)+'_'+str(total_iterations), log_ps[chain])

# Check convergence and continue sampling if not converged

GR = Gelman_Rubin(sampled_params)
print('At iteration: ', total_iterations, ' GR = ', GR)
np.savetxt(parfn + 'GR_iteration_' + str(total_iterations) + '.txt', GR)

old_samples = sampled_params
if np.any(GR > 1.2) or total_iterations < 25000:
    while not converged:
        total_iterations += niterations
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]

        sampled_params, log_ps = run_dream(sampled_parameter_names,
                                           likelihood=likelihood,
                                           niterations=niterations,nchains=nchains,
                                           restart=True, start=starts, start_random=False,
                                           save_history=True, history_thin=10,
                                           adapt_gamma=True, adapt_crossover=True,
                                           model_name=parfn,
                                           verbose=True)

        for chain in range(len(sampled_params)):
            np.save(parfn + str(chain) + '_' + str(total_iterations),
                        sampled_params[chain])
            np.save(logpsfn + str(chain) + '_' + str(total_iterations),
                        log_ps[chain])

        old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
        GR = Gelman_Rubin(old_samples)
        print('At iteration: ', total_iterations, ' GR = ', GR)
        np.savetxt(parfn + 'GR_iteration_' + str(total_iterations)+'.txt', GR)

        if np.all(GR < 1.2):
            converged = True
            print('Optimization has converged, all done!')

# try:
    # # Plot output
    # total_iterations = len(old_samples[0])
    # burnin = total_iterations // 2
    # samples = old_samples
    # samples = np.concatenate((old_samples[0][burnin:, :], old_samples[1][burnin:, :],
    #                           old_samples[2][burnin:, :], old_samples[3][burnin:, :],
    #                           old_samples[4][burnin:, :]))

    # ndims = len(old_samples[0][0])
    # colors = sns.color_palette(n_colors=ndims)
    # for dim in range(ndims):
    #     fig = plt.figure()
    #     sns.distplot(samples[:, dim], color=colors[dim])
        #fig.savefig('WM01_SL_' + str(dim))

# except ImportError:
    # pass

# else:
# run_kwargs = {'parameters':iBB.sampled_parameter_names, 'likelihood':iBB.likelihood,
#               'niterations':10000, 'nchains':nchains, 'multitry':False,
#               'gamma_levels':4, 'adapt_gamma':True,
#               'history_thin':1,
#               'model_name':parfn, 'verbose':True}