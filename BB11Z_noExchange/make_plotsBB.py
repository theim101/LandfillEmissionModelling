#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:59:18 2020

@author: theimovaara
"""

#Plot measurements in plots
import matplotlib.pyplot as plt
import seaborn as sns
from Eval_SL_BB11Z_DREAM02_longterm import *
loc_name = 'Braambergen 11Z'



loc_name = 'Braambergen 11Z'

plt.close('all')

fig0 = plt.figure(0)
plt.plot(log_ps.T,'.')
plt.title('log_ps')


fig1 = plt.figure(1)
plt.plot(iBB.meas_data_flow.index,iBB.meas_data_flow['totF'].values,'o')
for ii in range(len(par_set)):
    plt.plot(iBB.tmeas, totF_sim[ii])
plt.xlabel('year')
plt.ylabel(r'cumulative leachate [$m^3$]')
plt.title(loc_name)


fig2 = plt.figure(2)
plt.clf()
plt.plot(iBB.tmeas,iBB.meas_data_flow.diff()/iBB.measFreq, 'g.')
for ii in range(len(par_set)):
    plt.plot(iBB.tmeas[1:], np.diff(totF_sim[ii])/iBB.measFreq)
plt.xlabel('year')
plt.ylabel(r'leachate pump rate [$m^3/day]$')
plt.title(loc_name)

fig3 = plt.figure(3)
plt.plot(iBB.lab_data.index,iBB.lab_data['Cl'],'g.')
for ii in range(len(par_set)):
    plt.plot(iBB.trange[iBB.tmeaslab_ind], cDr_tot[ii][iBB.tmeaslab_ind])
plt.xlabel('year')
plt.ylabel(r'$[{Cl}^-]$ [mg/L]')
plt.title(loc_name)

n2006 = np.where(iBB.trange=='2006-01-01')[0][0]
fig4 = plt.figure(4)
plt.clf()
for ii in range(len(par_set)):
    totF = (cumqDr_tot[ii][n2006:] - cumqDr_tot[ii][iBB.tmeas_ind[0]])
    plt.plot(iBB.trange[n2006:],np.diff(totF))
             
plt.title('long term, water flow')
plt.xlabel('year')
plt.ylabel('leachate pumprate [$m^3/d$]')

fig5 = plt.figure(5)
for ii in range(len(par_set)):
    plt.plot(iBB.trange[n2006:],cDr_tot[ii][n2006:])
plt.title('long term, outflow concentration')
plt.xlabel('year')
plt.ylabel(r'$[{Cl}^-]$ [mg/L]')

fig6 = plt.figure(6)
[plt.plot(iBB.trange[n2006:],mW_tot[ii][n2006:,-1]/vW_tot[ii][n2006:,-1]) for ii in range(len(par_set))]
plt.title('long term, bulk concentration')
plt.xlabel('year')
plt.ylabel(r'$[{Cl}^-]$ [mg/L]')

fig7 = plt.figure(7)
[plt.plot(iBB.trange[n2006:],mW_tot[ii][n2006:,-1]) for ii in range(len(par_set))]
plt.title('long term, bulk Mass')
plt.xlabel('year')
plt.ylabel(r'$M_{{Cl}^-}$ [kg/m2]')

fig8 = plt.figure(8)
plt.clf()
[plt.plot(iBB.trange[n2006:],np.sum(mW_tot[ii][n2006:,:],axis=1)) for ii in range(len(par_set))]
plt.title('long term, total mass in wastebody')
plt.xlabel('year')
plt.ylabel(r'$M_{{Cl}^-} [kg/m^2]$')

fig9 = plt.figure(9)
plt.clf()
[plt.plot(iBB.trange[n2006:],vW_tot[ii][n2006:,-1]) for ii in range(len(par_set))]
plt.title('long term, bulk storage')
plt.xlabel('year')
plt.ylabel('Vw [m]')


fig10 =  plt.figure(10)
for ii in range(len(par_set)):
    vWtmp = np.arange(0,9,0.01)
    bFlow = emLF_tot[ii].BaseFlow (vWtmp)
    plt.plot(vWtmp,(bFlow))
plt.title('long term, base flow')
plt.xlabel('vW [m]')
plt.ylabel('base flow [m/d]')


best_par = pd.Series(par,index=iBB.par_df.columns)
allpardf = par_setdf.append(best_par,ignore_index=True)
print(allpardf.T)



pset = pd.DataFrame()
for ii in range(len(par_set)):
    pset[ii] = pd.Series(par_set[ii],index=iBB.par_df.columns) 

plot_chains = False
if plot_chains:
    plt.close('all')
    for ii in np.arange(npar):
        plt.figure()
        plt.clf()
        plt.plot(sampled_params[:,:,ii].T,'.')
        plt.title(iBB.par_df.columns[ii])
      
    #plt.close(fig4)


plot_dists = False
if plot_dists:
    ndims = len(sampled_params[0][0])
    nsamples = len(sampled_params[0])
    nchains = len(sampled_params)
    nburnin = nsamples//2
    plt.close('all')
    colors = sns.color_palette(n_colors=ndims)
    samp_par_reshaped = sampled_params[:, nburnin:, :].reshape((nsamples-nburnin)*nchains,ndims)
    par_mean = samp_par_reshaped.mean(axis=0)
    par_std = samp_par_reshaped.std(axis=0)
    
    par_mean = pd.Series(par_mean,index=iBB.par_df.columns)
    par_std = pd.Series(par_std,index=iBB.par_df.columns)
    par_stats = pd.DataFrame(par_mean,columns=['mean'])
    par_stats['std'] = par_std
    par_stats['best']= best_par
    
    #figdist, axdist = plt.subplots(1, 1, figsize=(5, 5))
    for dim in range(ndims):
        #fig = plt.figure()
        sns.displot(samp_par_reshaped[:,dim], color=colors[dim],stat='probability',bins=25)
        plt.title(iBB.par_df.columns[dim])
        #fig.savefig('WM01_SL_' + str(dim))

# calculate mean and stddev of parameters


plot_pairs = False
if plot_pairs:
    plt.close('all')
    ndims = len(sampled_params[0][0])
    nel = 4000
    lpars = sampled_params[:, -nel:, :].reshape(5*nel,ndims).T
    df_par = pd.DataFrame(lpars).T
    df_par.columns = iBB.par_df.columns
    #fig = sns.pairplot(df_par,markers='.')
    #fig.savefig('WM01_SL_PairPlot')
    ccoef = np.corrcoef(lpars)
    hicor = np.where((np.abs(ccoef) > 0.85) & (np.abs(ccoef) < 0.999))
    for ii in range(len(hicor[0])):
        plt.figure(ii)
        plt.plot(sampled_params[:,-nel:,hicor[0][ii]].T,
             sampled_params[:,-nel:,hicor[1][ii]].T,'.' )
        plt.xlabel(iBB.par_df.columns[hicor[0][ii]])
        plt.ylabel(iBB.par_df.columns[hicor[1][ii]])

# Interprete the results with respect to the baseline measurements
fDM = 0.773 # fraction of dry matter for WM
# dry weight per m2 landfill [kg/m2]
dryWeightpm2 = iBB.lF.wetWeight * fDM / iBB.lF.surfArea

# dry density landfill
dryDensity = dryWeightpm2 / iBB.lF.height

print('dryWeightpm2 =', dryWeightpm2)
print('dryDensity = ', dryDensity)

# plot development in chloride mass per kg dry weight
fig11 = plt.figure(11)
plt.clf()
[plt.plot(iBB.trange[n2006:],np.sum(mW_tot[ii][n2006:,:],axis=1)/dryWeightpm2) for ii in range(len(par_set))]
plt.title('long term, total mass in wastebody')
plt.xlabel('year')
plt.ylabel(r'$M_{{Cl}^-} [g/kg]$')

fig12, ax12 = plt.subplots(1, 1, figsize=(5, 5))
mflow = iBB.meas_data_flow.diff()[1:].T.values[0]
ndims = len(par_set)
colors = sns.color_palette(n_colors=ndims)
stats = np.zeros([len(par_set),2])
for ii in range(len(par_set)):
    sflow = np.diff(totF_sim[ii])
    sns.histplot(mflow-sflow,kde=True,ax=ax12,color=colors[ii])
    stats[ii, 0] = np.mean(mflow-sflow)
    stats[ii, 1] = np.std(mflow-sflow)
plt.xlabel('delta')
plt.ylabel(r'residuals leachate pump rate [$m^3/(7 days)$]')
print(stats)

fig13, ax13 = plt.subplots(1, 1, figsize=(5, 5))
mconc = iBB.lab_data['Cl'][1:].T.values
stats_c = np.zeros([len(par_set),2])
for ii in range(len(par_set)):
    sconc = cDr_tot[ii][iBB.tmeaslab_ind][1:]
    sns.histplot(mconc-sconc,color=colors[ii],kde=True,ax=ax13)
    stats[ii, 0] = np.mean(mconc-sconc)
    stats[ii, 1] = np.std(mconc-sconc)
plt.xlabel('delta')
plt.ylabel(r'residuals leachate concentration [$mg/L$]')
print(stats)

#plt.close(fig14)
fig14,ax14 = plt.subplots(2,1,sharex=True)
# plot qinf
[ax14[0].plot(iBB.trange[n2006+1:],qInf_tot[ii][n2006:]) for ii in range(len(par_set))]
     
[ax14[1].plot(iBB.trange[n2006:], np.maximum(0,(thCL_tot[ii][n2006:]-emLF_tot[ii].cL.thRes)
                             / (emLF_tot[ii].cL.Por - emLF_tot[ii].cL.thRes)))
              for ii in range(len(par_set))]
ax14[0].set_ylabel('qInf')
ax14[1].set_ylabel('Seff')
ax14[1].set_xlabel('time')


# Plot forward travel time pdf
tage = np.arange(0.0001,5*365,0.1)
rv1 = np.zeros([len(par_set),len(tage)])
rv2 = np.zeros([len(par_set),len(tage)])
cdf = np.zeros([len(par_set),len(tage)])
pdf = np.zeros([len(par_set),len(tage)])
for ii in range(len(par_set)):
    sig1 = 10**allpardf['wBsig1'][ii]
    sig2 = 10**allpardf['wBsig2'][ii]
    tau1 = allpardf['wBtau1'][ii]
    tau2 = allpardf['wBtau2'][ii]
    rv1 = scpstats.lognorm(sig1,loc=0,scale=tau1)
    rv2 = scpstats.lognorm(sig2,loc=0,scale=tau2)
    cdf[ii] = allpardf['wBmfrac'][ii]*rv1.cdf(tage) + \
              (1-allpardf['wBmfrac'][ii])*rv2.cdf(tage)
    pdf[ii] = allpardf['wBmfrac'][ii]*rv1.pdf(tage) + \
              (1-allpardf['wBmfrac'][ii])*rv2.pdf(tage)

ffig15,ax15 = plt.subplots(1,1)
# plot qinf
[ax15.plot(tage,pdf[ii]) for ii in range(len(par_set))]
              
ax15.set_ylabel('pdf')
ax15.set_xlabel('tage')

fig16,ax16 = plt.subplots(1,1)
# plot qinf
[ax16.plot(tage,cdf[ii]) for ii in range(len(par_set))]
              
ax16.set_ylabel('cdf')
ax16.set_xlabel('tage')

fig17 = plt.figure(17)
plt.clf()
[plt.plot(iBB.trange[n2006:],np.sum(vW_tot[ii][n2006:,:],axis=1)) for ii in range(len(par_set))]
plt.title('long term, total storage')
plt.xlabel('year')
plt.ylabel('Vw [m]')