#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:59:18 2020

@author: theimovaara
"""

#Plot measurements in plots
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Eval_SL_KR3_DREAM02_longterm import *

loc_name = 'Kragge 3'

plt.close('all')

fig0 = plt.figure(0)
plt.plot(log_ps.T,'.')
plt.title('log_ps')


fig1 = plt.figure(1)
plt.clf()
plt.plot(iKR.meas_data_flow.index,iKR.meas_data_flow['level_data'].values,'.')
for ii in range(len(par_set)):
    plt.plot(iKR.trange[:], drainLev_tot[ii][1:])
plt.xlabel('year')
plt.ylabel(r'level in drainage system [$m$]')
plt.title(loc_name)


fig2 = plt.figure(2)
plt.clf()
for ii in range(len(par_set)):
    aap = np.cumsum(drainOutF_tot[ii])
    plt.plot(iKR.meteo_data.index,aap - aap['drainRate'].iloc[iKR.tmeas_ind[0]])

iKR.meas_data_flow['totF'].plot()

plt.xlabel('year')
plt.ylabel(r'cumulative Leachate extracted $[m^3]$')
#plt.legend(['measured','simulated'])
plt.title(loc_name)


fig3 = plt.figure(3)
plt.plot(iKR.lab_data.index,iKR.lab_data['Cl'],'g.')
for ii in range(len(par_set)):
    plt.plot(iKR.trange[:], cDr_tot[ii])
plt.xlabel('year')
plt.ylabel(r'$[{Cl}^-]$ [mg/L]')
plt.title(loc_name)


fig4 = plt.figure(4)
plt.clf()
for ii in range(len(par_set)):
    totF = (cumqDr_tot[ii][:] - cumqDr_tot[ii][iKR.tmeas_ind[0]])
    plt.plot(iKR.trange[:],np.diff(totF))

plt.title('long term, water flow')
plt.xlabel('year')
plt.ylabel('leachate pumprate [$m^3/d$]')

fig5 = plt.figure(5)
for ii in range(len(par_set)):
    plt.plot(iKR.trange[:],cDr_tot[ii])
plt.title('long term, outflow concentration')
plt.xlabel('year')
plt.ylabel(r'$[{Cl}^-]$ [mg/L]')

fig6 = plt.figure(6)
[plt.plot(iKR.trange[:],mW_tot[ii][:,-1]/vW_tot[ii][:,-1]) for ii in range(len(par_set))]
plt.title('long term, bulk concentration')
plt.xlabel('year')
plt.ylabel(r'$[{Cl}^-]$ [mg/L]')

fig7 = plt.figure(7)
[plt.plot(iKR.trange[:],mW_tot[ii][:,-1]) for ii in range(len(par_set))]
plt.title('long term, bulk Mass')
plt.xlabel('year')
plt.ylabel(r'$M_{{Cl}^-}$ [kg/m2]')

fig8 = plt.figure(8)
plt.clf()
[plt.plot(iKR.trange[:],np.sum(mW_tot[ii][:,:],axis=1)) for ii in range(len(par_set))]
plt.title('long term, total mass in wastebody')
plt.xlabel('year')
plt.ylabel(r'$M_{{Cl}^-} [kg/m^2]$')

fig9 = plt.figure(9)
plt.clf()
[plt.plot(iKR.trange[:],vW_tot[ii][:,-1]) for ii in range(len(par_set))]
plt.title('long term, bulk storage')
plt.xlabel('year')
plt.ylabel('Vw [m]')


fig10 =  plt.figure(10)
plt.clf()

for ii in range(len(par_set)):
    vWtmp = np.arange(0,9,0.01)
    bFlow = emLF_tot[ii].BaseFlow (vWtmp)
    plt.plot(vWtmp,(bFlow))
plt.title('long term, base flow')
plt.xlabel('vW [m]')
plt.ylabel('base flow [m/d]')

best_par = pd.Series(par,index=iKR.par_df.columns)
allpardf = par_setdf.append(best_par,ignore_index=True)
print(allpardf.T)



pset = pd.DataFrame()
for ii in range(len(par_set)):
    pset[ii] = pd.Series(par_set[ii],index=iKR.par_df.columns)

plot_chains = False
if plot_chains:
    plt.close('all')
    for ii in np.arange(npar):
        plt.figure()
        plt.clf()
        plt.plot(sampled_params[:,:,ii].T,'.')
        plt.title(iKR.par_df.columns[ii])

    #plt.close(fig4)


plot_dists = False
if plot_dists:
    ndims = len(sampled_params[0][0])
    nsamples = len(sampled_params[0])
    nchains = len(sampled_params)
    #nburnin = nsamples//2
    nburnin = 20000
    plt.close('all')
    colors = sns.color_palette(n_colors=ndims)
    samp_par_reshaped = sampled_params[:, nburnin:, :].reshape((nsamples-nburnin)*nchains,ndims)
    par_mean = samp_par_reshaped.mean(axis=0)
    par_std = samp_par_reshaped.std(axis=0)

    par_mean = pd.Series(par_mean,index=iKR.par_df.columns)
    par_std = pd.Series(par_std,index=iKR.par_df.columns)
    par_stats = pd.DataFrame(par_mean,columns=['mean'])
    par_stats['std'] = par_std
    par_stats['best']= best_par

    #figdist, axdist = plt.subplots(1, 1, figsize=(5, 5))
    for dim in range(ndims):
        #fig = plt.figure()
        sns.displot(samp_par_reshaped[:,dim], color=colors[dim],stat='probability',bins=25)
        plt.title(iKR.par_df.columns[dim])
        #fig.savefig('WM01_SL_' + str(dim))

# calculate mean and stddev of parameters


plot_pairs = False
if plot_pairs:
    plt.close('all')
    ndims = len(sampled_params[0][0])
    nsamples = len(sampled_params[0])
    nchains = len(sampled_params)
    nburnin = nsamples//2
    nel = nsamples-nburnin
    plt.close('all')
    colors = sns.color_palette(n_colors=ndims)
    lpars = sampled_params[:, nburnin:, :].reshape((nsamples-nburnin)*nchains,ndims)
    df_par = pd.DataFrame(lpars)
    df_par.columns = iKR.par_df.columns
    #fig = sns.pairplot(df_par,markers='.')
    #fig.savefig('WM01_SL_PairPlot')
    ccoef = np.corrcoef(lpars.T)
    hicor = np.where((np.abs(ccoef) > 0.85) & (np.abs(ccoef) < 0.999))
    for ii in range(len(hicor[0])):
        plt.figure(ii)
        plt.plot(sampled_params[:,-nel:,hicor[0][ii]].T,
             sampled_params[:,-nel:,hicor[1][ii]].T,'.' )
        plt.xlabel(iKR.par_df.columns[hicor[0][ii]])
        plt.ylabel(iKR.par_df.columns[hicor[1][ii]])

# Interprete the results with respect to the baseline measurements
fDM = 0.773 # fraction of dry matter for WM
# dry weight per m2 landfill [kg/m2]
dryWeightpm2 = iKR.lF.wetWeight * fDM / iKR.lF.surfArea

# dry density landfill
dryDensity = dryWeightpm2 / iKR.lF.height

print('dryWeightpm2 =', dryWeightpm2)
print('dryDensity = ', dryDensity)

# plot development in chloride mass per kg dry weight
fig11 = plt.figure(11)
plt.clf()
[plt.plot(iKR.trange[:],np.sum(mW_tot[ii][:,:],axis=1)/dryWeightpm2) for ii in range(len(par_set))]
plt.title('long term, total mass in wastebody')
plt.xlabel('year')
plt.ylabel(r'$M_{{Cl}^-} [g/kg]$')

fig12, ax12 = plt.subplots(1, 1, figsize=(5, 5))
mlev = iKR.meas_data_flow['level_data'].values[1:].T
ndims = len(par_set)
colors = sns.color_palette(n_colors=ndims)
stats = np.zeros([len(par_set),2])
for ii in range(len(par_set)):
    slev = drainLev_tot[ii][iKR.tmeas_ind][1:]
    sns.histplot(mlev-slev,kde=True,ax=ax12,color=colors[ii])
    stats[ii, 0] = np.mean(mlev-slev)
    stats[ii, 1] = np.std(mlev-slev)
plt.xlabel('delta')
plt.ylabel(r'residuals drainage level$]')
print(stats)

fig13, ax13 = plt.subplots(1, 1, figsize=(5, 5))
mconc = iKR.lab_data['Cl'][1:].T.values
stats_c = np.zeros([len(par_set),2])
for ii in range(len(par_set)):
    sconc = cDr_tot[ii][iKR.tmeaslab_ind][1:]
    sns.histplot(mconc-sconc,color=colors[ii],kde=True,ax=ax13)
    stats[ii, 0] = np.mean(mconc-sconc)
    stats[ii, 1] = np.std(mconc-sconc)
plt.xlabel('delta')
plt.ylabel(r'residuals leachate concentration [$mg/L$]')
print(stats)

#plt.close(fig14)
# fig14,ax14 = plt.subplots(3,1,sharex=True)
# # plot qinf
# [ax14[0].plot(iKR.trange[1:],qInf_tot[ii][0]) for ii in range(len(par_set))]
# [ax14[1].plot(iKR.trange[1:],cLOut_tot[ii][3][1:]) for ii in range(len(par_set))]
# [ax14[1].plot(iKR.trange[1:],cLOut_tot[ii][9][:]) for ii in range(len(par_set))]

# ax14[0].set_ylabel('qInf')
# ax14[1].set_ylabel('Seff')
# ax14[2].set_ylabel('qrunoff')
# ax14[2].set_xlabel('time')
#plt.close(fig14)
fig14,ax14 = plt.subplots(2,1,sharex=True)
# plot qinf
[ax14[0].plot(iKR.trange[1:],qInf_tot[ii]) for ii in range(len(par_set))]
     
[ax14[1].plot(iKR.trange[:], np.maximum(0,(thCL_tot[ii]-emLF_tot[ii].cL.thRes)
                             / (emLF_tot[ii].cL.Por - emLF_tot[ii].cL.thRes)))
              for ii in range(len(par_set))]
ax14[0].set_ylabel('qInf')
ax14[1].set_ylabel('Seff')
ax14[1].set_xlabel('time')

# Plot forward travel time pdf
tage = np.arange(0.0001,25*365,0.1)
rv1 = np.zeros([len(par_set),len(tage)])
rv2 = np.zeros([len(par_set),len(tage)])
rv3 = np.zeros([len(par_set),len(tage)])
rv4 = np.zeros([len(par_set),len(tage)])
cdf1 = np.zeros([len(par_set),len(tage)])
pdf1 = np.zeros([len(par_set),len(tage)])
cdf2 = np.zeros([len(par_set),len(tage)])
pdf2 = np.zeros([len(par_set),len(tage)])

for ii in range(len(par_set)):
    sig1 = 10**allpardf['wBsig1'][ii]
    sig2 = 10**allpardf['wBsig2'][ii]
    sig3 = 10**allpardf['wBsig3'][ii]
    sig4 = 10**allpardf['wBsig4'][ii]
    tau1 = allpardf['wBtau1'][ii]
    tau2 = allpardf['wBtau2'][ii]
    tau3 = allpardf['wBtau3'][ii]
    tau4 = allpardf['wBtau4'][ii]
    rv1 = scpstats.lognorm(sig1,loc=0,scale=tau1)
    rv2 = scpstats.lognorm(sig2,loc=0,scale=tau2)
    rv3 = scpstats.lognorm(sig3,loc=0,scale=tau3)
    rv4 = scpstats.lognorm(sig4,loc=0,scale=tau4)
    cdf1[ii] = allpardf['wBmfrac'][ii]*rv1.cdf(tage) + \
              (1-allpardf['wBmfrac'][ii])*rv2.cdf(tage)
    pdf1[ii] = allpardf['wBmfrac'][ii]*rv1.pdf(tage) + \
              (1-allpardf['wBmfrac'][ii])*rv2.pdf(tage)
    cdf2[ii] = allpardf['wBmfrac2'][ii]*rv3.cdf(tage) + \
              (1-allpardf['wBmfrac2'][ii])*rv4.cdf(tage)
    pdf2[ii] = allpardf['wBmfrac2'][ii]*rv3.pdf(tage) + \
              (1-allpardf['wBmfrac2'][ii])*rv4.pdf(tage)

fig15,ax15 = plt.subplots(2,1, sharex=True)
# plot qinf
[ax15[0].plot(tage,pdf1[ii]) for ii in range(len(par_set))]
[ax15[1].plot(tage,pdf2[ii]) for ii in range(len(par_set))]

[ax15[ii].set_ylabel('pdf' + str(ii+1)) for ii in range(2)]
ax15[1].set_xlabel('tage')

fig16,ax16 = plt.subplots(2,1,sharex=True)
# plot qinf
[ax16[0].plot(tage,cdf1[ii]) for ii in range(len(par_set))]
[ax16[1].plot(tage,cdf2[ii]) for ii in range(len(par_set))]

[ax16[ii].set_ylabel('cdf' + str(ii+1)) for ii in range(2)]
ax16[1].set_xlabel('tage')

fig17 = plt.figure(17)
plt.clf()
[plt.plot(iKR.trange[:],np.sum(vW_tot[ii][:,:],axis=1)) for ii in range(len(par_set))]
plt.title('long term, total storage')
plt.xlabel('year')
plt.ylabel('Vw [m]')


plot_cl = False
if plot_cl:
    plt.close('all')
    #qinf, qrf, qev, seff, th, mbalerr, cInf, m = cLOut_tot

    fig1,ax1 = plt.subplots(2,1,sharex=True)
    # plot qinf
    [ax1[0].plot(iKR.trange[1:],cLOut_tot[ii][0]) for ii in range(len(par_set))]
    [ax1[1].plot(iKR.trange[1:],cLOut_tot[ii][3][1:]) for ii in range(len(par_set))]
    ax1[0].set_ylabel('qInf')
    ax1[1].set_ylabel('Seff')
    ax1[1].set_xlabel('time')
    

fig18 = plt.figure(18)
plt.clf()
#plt.plot(iKR.meas_data_flow.index,iKR.meas_data_flow['level_data'].values,'.')
for ii in range(len(par_set)):
    plt.plot(drainVolN_tot[ii][iKR.tdseries.tmeas_ind], 
             drainVolN_tot[ii][iKR.tdseries.tmeas_ind] /
             iKR.meas_data_flow['level_data'].values  
             ,'.')
    shape_fact = emLF_tot[ii].wB.drApar * (drainVolN_tot[ii] >= emLF_tot[ii].wB.drCpar) \
        + (emLF_tot[ii].wB.drBpar + (emLF_tot[ii].wB.drApar-emLF_tot[ii].wB.drBpar) \
           /emLF_tot[ii].wB.drCpar * drainVolN_tot[ii]) * (drainVolN_tot[ii] < emLF_tot[ii].wB.drCpar) 
    
    plt.plot(drainVolN_tot[ii], 
             shape_fact  
             ,'-')
plt.xlabel('drainVolN')
plt.ylabel(r'drainVolN/drainLev]')
plt.title(loc_name)

fig19 = plt.figure(19)
plt.clf()
#plt.plot(iKR.meas_data_flow.index,iKR.meas_data_flow['level_data'].values,'.')
for ii in range(len(par_set)):
    plt.plot(drainVolN_tot[ii][iKR.tdseries.tmeas_ind], 
             iKR.meas_data_flow['level_data'].values  
             ,'.')
plt.xlabel('drainVolN')
plt.ylabel(r'levDrain]')
plt.title(loc_name)


fig20,ax20 = plt.subplots(3,1,sharex=True)
ax20[0].plot(iKR.meas_data_flow.index,iKR.meas_data_flow['level_data'].values,'.')
ax20[2].plot(iKR.meas_data_flow.index,iKR.meas_data_flow['totF'].diff(21).values/21,'.')

for ii in range(len(par_set)):
    ax20[1].plot(iKR.meas_data_flow.index,drainVolN_tot[ii][iKR.tdseries.tmeas_ind] 
             ,'.')
ax20[1].set_xlabel('time')
ax20[0].set_ylabel('drainLev')
ax20[1].set_ylabel('drainVolN]')
ax20[2].set_ylabel('flow')
ax20[0].set_title(loc_name)

fig21,ax21 = plt.subplots(2,1,sharex=True)
ax21[0].plot(iKR.meas_data_flow.index,iKR.meas_data_flow['level_data'].values,'.')
ax21[1].plot(iKR.meas_data_flow.index,iKR.meas_data_flow['totF'].diff(21).values/21,'.')

# for ii in range(len(par_set)):
#     ax20[1].plot(iKR.meas_data_flow.index,drainVolN_tot[ii][iKR.tdseries.tmeas_ind] 
#              ,'.')
ax21[1].set_xlabel('time')
ax21[0].set_ylabel('drainLev')
# ax20[1].set_ylabel('drainVolN]')
ax21[1].set_ylabel('flow')
ax21[0].set_title(loc_name)