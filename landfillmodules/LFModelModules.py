#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 09:46:35 2018

@author: theimovaara
"""

import numpy as np
#import scipy.integrate as spi
import scipy.stats as stats
import pandas as pd
import scipy.special as spspec

import matplotlib.pyplot as plt
import seaborn as sns

class Emission_LF:
    """ Class for simulating water and mass balance in landfill waste bodies.
        Model allows of recirculation of leachate.
        Boudary conditions are controlled by rainfall, potential evapotranspiration,
        leachate infiltrated, concentration in infiltrated leachate, leachate extracted.
    """
    
    #TO DO: implement concentration in infiltrating leachate (recirculation)
    #       mass in coverlayer increases due to infiltrating leachate
    #       mass in cover layer decreases with runoff and infiltration to waste
    #       body
    #       Assumption: All flows leaving cover layer have same concentration!
    
    def __init__(self, tRange, iniSt, cL, wB, lF, meteo_data, infilF):
    
        # Initialize emission models
        self.tRange = tRange
        self.dt = np.diff(tRange)
        iniSt = iniSt
        self.cL = cL
        self.wB = wB
        self.lF = lF
        meteo_data = meteo_data
        self.infilF = infilF
        
        self.ndays = np.int(tRange[-1]-tRange[0])
        ndays = self.ndays
        self.dt = np.diff(tRange)
        
        ## Variables for Cover Layer model
        self.th = np.zeros([ndays+1])
        self.seff = np.zeros([ndays+1])
        self.cInf = np.zeros([ndays+1])
        self.m = np.zeros([ndays+1])
        self.qInf = np.zeros([ndays])
        self.qrf = np.zeros([ndays])
        self.qev = np.zeros([ndays])
        self.qinfil = np.zeros([ndays])
        self.cinfil = np.zeros([ndays+1])
        self.qrunoff = np.zeros([ndays])
        self.mbalerr = np.zeros([ndays+1])

        self.thMin = 0.001
        #save initial states
        self.th[0] = iniSt.thIni
        self.m[0] = iniSt.cLcIni*self.th[0]
        
        # select rainfall and evaporation data from meteo_data
        meteo_time = meteo_data.index.astype(np.int64)/(1e9 * 3600 * 24)
        xy, md_ind, t_ind = np.intersect1d(meteo_time, np.ceil(tRange), return_indices=True)

        self.qrf = -meteo_data['rain_station'].iloc[md_ind].values
        self.qev = meteo_data['pEV'].iloc[md_ind].values * cL.cropFact
        
   
        # expand inifilF into the time series...
        if len(infilF) > 0:
            qInf2 = pd.DataFrame(tRange,columns=['DateTime'])
            qInf2.set_index('DateTime',inplace=True)
            qInf2['infilRate'] = 0

            tValsInf = infilF.index.astype(np.int64)/(1e9 * 3600 * 24)

            idx1 = np.where((qInf2.index >= tValsInf.min()) & (qInf2.index <= tValsInf.max()))
            infilF_resample = infilF.resample('D').sum()
            infilF_resample = infilF_resample.interpolate()
            infilRate = infilF_resample.diff()
            qInf2['infilRate'].iloc[idx1[0][1:]] = infilRate.iloc[1:].values/lF.surfArea
        else:
            qInf2 = pd.DataFrame(tRange,columns=['DateTime'])
            qInf2.set_index('DateTime',inplace=True)
            qInf2['infilRate'] = 0

    
        self.qinfil = -qInf2['infilRate'].iloc[md_ind].values
    
        
        ## Initialize Wastebody Emission Model
        # number of mobile cells...
        ndEx = int(wB.ndEx)
        tage = np.arange(ndEx)
        
        rv1 = stats.lognorm(wB.sig1,loc=0,scale=wB.tau1)
        rv2 = stats.lognorm(wB.sig2,loc=0,scale=wB.tau2)
        
        cdf1 = rv1.cdf(tage)
        cdf2 = rv2.cdf(tage)
        self.cdf = wB.mfrac*cdf1 + (1-wB.mfrac)*cdf2
        self.pdf = np.concatenate(([0],np.diff(self.cdf)))
    
        rv3 = stats.lognorm(wB.sig3,loc=0,scale=wB.tau3)
        rv4 = stats.lognorm(wB.sig4,loc=0,scale=wB.tau4)
        cdf3 = rv3.cdf(tage)
        cdf4 = rv4.cdf(tage)
        self.cdf_infil = wB.mfrac2*cdf3 + (1-wB.mfrac2)*cdf4
        self.pdf_infil = np.concatenate(([0],np.diff(self.cdf_infil)))
    
        self.vW = np.zeros([ndays+1,ndEx+1])
        self.mW = np.zeros([ndays+1,ndEx+1])

        # Initialize water volume and solute mass present in waste body...
        self.vW[0,0:ndEx] = wB.bFlow #iniSt.wBvWIni * wB.mfrac / ndEx
        self.vW[0,ndEx] = iniSt.wBvWIni #np.max([0,iniSt.wBvWIni - sum(vW[0,0:ndEx])])
        
        self.mW[0,0:ndEx] = self.vW[0,0:ndEx] * iniSt.wBcIni
        self.mW[0,ndEx] = self.vW[0,ndEx] * iniSt.wBcIni


    def BaseFlow (self, vWBulk):
        """ Calculate base flow a function of bulk water storage...
        """
        f = (self.wB.bFlow*0.5*(1.0+spspec.erf(self.wB.alphavW*(vWBulk-self.wB.vWK))))*(vWBulk>self.wB.vWRes)
        return f


    def ExRates(t, mW, vW, kEx, ndEx, bFlow):
        """ Calculate the solute exchange between cell and bulk...
        """
        totM = np.sum(mW)
        totV = np.sum(vW)
        dMdt = np.zeros(ndEx+1)
        jj = np.arange(0, ndEx)
        dMdt[jj] = kEx*(vW[jj]*totM/totV - mW[jj].squeeze())
        dMdt[ndEx] = -np.sum(dMdt[jj])
        return dMdt


    def Emission_CL(self, th, m, qrf, qev, dt):
        """ Here we define the doc_string
        Function for calulating water balance of cover layer
        For every # Function for calculating water flows in coverlayer
        th: volumetric water content in layer [-]
        m: mass in layer kg/m2
        qrf: rainfall rate m3/(m2 day)
        qev: potential evapo transpiration rate m3/(m2 day)
        qinfil: leachate infiltration rate m3/(m2 day)
        dt: time step (day)
        """
        thMin = 0.001
        qrunoff = 0
        seff = np.maximum(0,(th-self.cL.thRes)/(self.cL.Por-self.cL.thRes))
        qInf = -self.cL.KSat * seff**self.cL.bPar
        thEst = th -(qrf - qInf + qev + qrunoff) * dt / self.cL.dzCL
        if thEst > self.cL.Por:
            #qInf too small, profile is too wet, excess rainfall should be routed
            #to infiltration flux
            thExcess = thEst-self.cL.Por
            
            qrunoff = self.cL.alphaRO * thExcess * self.cL.dzCL / dt
            qInfFast = (1.0 - self.cL.alphaRO) * thExcess * self.cL.dzCL / dt
            qInf = qInf - qInfFast
            #self.qinf[ii] = self.qinf[ii] 
            thEst = th - (qrf - qInf + qev + qrunoff) \
                * dt / self.cL.dzCL
        else:
            if (thEst < thMin) and (th > self.cL.thRes):
                #qInf is too large and possible qEvap as well We start by reducing
                #qInf
                qInf = -(th-self.cL.thRes)*self.cL.dzCL / dt
                #restimate thEst
                thEst = th - (qrf - qInf + qev) * dt / self.cL.dzCL

            if thEst < thMin:
                #only qEv is too large and needs to be reduce
                qev = qev + thEst * self.cL.dzCL / dt
                thEst = th-(qrf  - qInf + qev) * dt / self.cL.dzCL
                
        thnew = np.maximum(thEst, thMin)
        mbalerr = thEst - (th - (qrf  - qInf + qev \
                                 + qrunoff) * dt / self.cL.dzCL)

        cInf = m/th
        mnew = m - (qrf * self.cL.cRain - qInf *cInf \
                    + qrunoff*cInf)* dt / self.cL.dzCL
    
        #seff[ii+1] = np.maximum(0,(th[ii+1]-self.cL.thRes)/(self.cL.Por-self.cL.thRes))
        #qinf[ii+1] = -cL.KSat * seff[ii+1]**cL.bPar
    
        cLOut = [qInf, thnew, mnew, cInf, mbalerr, qrunoff]
        return cLOut


    def Emission_WB (self, vW, mW, qInf, cInf, qinfil, cinfil):
            """ 
            Function for calulating water and mass fluxes of Waste body for an
            individual time step
            waste body discretized in ndEx cells (where last cell is the bulk)
            Each cell contains water jj days from being removed from the waste body 
            as pumped leachate.
            
            vW[0:ndEx]: volume of water stored in waste body m3/m2 
            mW 0:ndEx]: mass stored in cell: kg/m2
            qInf infiltration rate m3/(m2 day)
            cInf: concentration in infiltrating water (kg/m3)
            dt: time step (day)
            """
            ndEx = np.int(self.wB.ndEx)
            #tage = np.arange(ndEx)

            
            vWNew = np.zeros(ndEx+1)
            mWNew = np.zeros(ndEx+1)
            jj = np.arange(0,ndEx-1)
    
            # solve water balance (water of ttravel time jj+1 + new version of rain)
    
            # shift all cells by one day, extract bF from bulk
            vWNew[jj] = vW[jj+1]
    
            bF = self.BaseFlow(vW[ndEx])
            vWNew[ndEx-1] = bF
            vWNew[ndEx] = vW[ndEx]-bF
    
            # distribute rain over cells
            vWNew[jj] = vWNew[jj] + qInf * self.pdf[jj+1] \
                + qinfil * self.pdf_infil[jj+1]
            
            vWNew[ndEx] = vWNew[ndEx] + qInf * (1 - self.cdf[-1]) \
                + qinfil * (1 - self.cdf_infil[-1])
    
            # distribute rain over cells
            #vW[ii+1,jj] = vW[ii+1,jj] + qInf2['infilRate'].iloc[ii]*pdf3[jj+1]
            #vW[ii+1,ndEx] = vW[ii+1,ndEx] + (1-cdf3[-1])*qInf2['infilRate'].iloc[ii]
    
    
            # solute balance rated dependent exchange
            # def intFun(t, y):
            #     outVals = exRates(t, y, vW[ii,:], wB.kEx, ndEx, wB.bFlow)
            #     return outVals
    
            # mWeX = spi.solve_ivp(intFun, [0, 1], mW[ii,:],
            #               method='RK45',vectorized=True)
            # mWtmp = mWeX.y[:,-1]
            # mWtmp[mWtmp<0] = 0
            # solute balance convective flow
            # mW[ii+1,jj] = mWtmp[jj+1] + qInf[ii]*pdf[jj+1]*cInf[ii] # + exch
            # solute balance convective flow
    
            mWNew[jj] = mW[jj+1] + qInf * self.pdf[jj+1] * cInf \
                + qinfil * self.pdf_infil[jj+1] * cinfil# + exch
            
            mWNew[ndEx] = mWNew[ndEx] + qInf * (1 - self.cdf[-1]) * cInf \
                + qinfil * (1 - self.cdf_infil[-1]) * cinfil
                
            # mass flowing from bulk to mobile
            if bF > 0:
                mWNew[ndEx-1] = bF * mW[ndEx] / vW[ndEx]
                # mW[ii+1,ndEx] = mW[ii,ndEx] - bF * mW[ii,ndEx]/vW[ii,ndEx]
    
            else:
                mWNew[ndEx-1] = 0
    
            #Update mass in bulk
            #Ensure that mass in bulk storage is never less than zero...
            mWNew[ndEx] = np.amax([mW[ndEx] - mWNew[ndEx-1],0])
            #mWNew[ndEx] = mW[ndEx] - mWNew[ndEx-1]
            
            return vWNew, mWNew
    


    def SimulateEmissionLF (self):
        
        for ii in np.arange(0,self.ndays):
            
            cLOut = self.Emission_CL(self.th[ii], self.m[ii], self.qrf[ii], 
                            self.qev[ii],  
                            self.dt[ii])
            
            self.qInf[ii], self.th[ii+1], self.m[ii+1], self.cInf[ii], \
                self.mbalerr[ii], self.qrunoff[ii] = cLOut 
            
            self.vW[ii+1], self.mW[ii+1] = self.Emission_WB (self.vW[ii], 
                self.mW[ii], -self.qInf[ii], self.cInf[ii], 
                -self.qinfil[ii], self.cinfil[ii])
            
            self.cinfil[ii+1] = self.mW[ii+1,0] / (self.vW[ii+1,0] + 1e-14)
            
            
    
        qDr = self.vW[:,0] * self.lF.surfArea
        cDr = self.mW[:,0] / (self.vW[:,0] + 1e-14)
    
        cumqDr = np.concatenate(([0],qDr),axis=0).cumsum()
    
        emOut = [cumqDr, cDr, qDr,  self.qInf, self.th, self.m, self.cInf,
                                       self.mbalerr, self.qrunoff, self.vW, 
                                       self.mW]
        
        return emOut
        

def likelihood_Emission_Recirculate_LF (par, tdseries, lF, meteo_data, infilF, meas_data_flow,
                      lab_data, sdData):
    # Run Forward model using parameters in par
    dzCL = par[0]
    cropFact = par[1]
    cLPor = par[2]
    cLthRes = par[3]*cLPor
    cLKSat = 10**par[4]
    cLbPar = par[5]

    wBtau1 = par[6]
    wBsig1 = 10**par[7]
    wBtau2 = wBtau1 + par[8]
    wBsig2 = 10**par[9]
    wBmfrac = par[10]
    wBvWIni = par[11] * (lF.height - dzCL)
    wBvWRes = par[12] * wBvWIni
    wBbFlow = 10**par[13]
    cLcIni = 10**par[14]
    wBcIni = 10**par[15]
    wBvWK = par[16]*(lF.height-dzCL)
    wBalphavW = 10**par[17]
    wBndEx = par[18]
    wBdrVIni = 10**par[19]
    wBdrApar = 10**par[20]
    wBdrBpar = 10**par[21]
    wBdrCpar = 10**par[22]
    cLalphaRO = par[23]
    wBtau3 = par[24]
    wBsig3 = 10**par[25]
    wBtau4 = wBtau3 + par[26]
    wBsig4 = 10**par[27]    
    wBmfrac2 = par[28]
    
    cLthIni = cLPor
    stdMeas_levData = sdData[0]
    stdMeas_cDr = sdData[1]

    coverLayer = {'dzCL': dzCL,
                  'cropFact': cropFact,
                  'Por': cLPor,
                  'thRes': cLthRes,
                  'KSat': cLKSat,
                  'bPar': cLbPar,
                  'cRain': 0.0,
                  'alphaRO': cLalphaRO}
    
    cL = pd.Series(coverLayer)


    wasteBody = {'tau1': wBtau1,
              'sig1': wBsig1,
              'tau2': wBtau2,
              'sig2': wBsig2,
              'mfrac': wBmfrac,
              'bFlow': wBbFlow,
              'vWRes': wBvWRes,
              'vWK': wBvWK,
              'alphavW': wBalphavW,
              'ndEx': wBndEx,
              'drApar': wBdrApar,
              'drBpar': wBdrBpar,
              'drCpar': wBdrCpar,
              'tau3': wBtau3,
              'sig3': wBsig3,
              'tau4': wBtau4,
              'sig4': wBsig4,
              'mfrac2': wBmfrac2
              }

    wB = pd.Series(wasteBody)


    #Calculate initial volume in drainage system. Please note that the bottom liner is installed at a slope

    initialStates = {'thIni': cLthIni,
                     'cumQIni': meas_data_flow['totF'].iloc[0],
                     'cLcIni': cLcIni,
                     'wBvWIni': wBvWIni,
                     'wBcIni': wBcIni,
                     'drVini': wBdrVIni
                     }

    iniSt = pd.Series(initialStates)

    trange_sim = tdseries.trange.astype(np.int64)/(1e9 * 3600 * 24)

    

    # Initialize the Emission Model Class
    
    emLF = Emission_LF(trange_sim.values, iniSt, cL, wB, lF, meteo_data,
                         infilF)
    # Call model: this function solves the problem...
    emOut = emLF.SimulateEmissionLF()
    
    cumqDr, cDr, qDr, qInf, th, m, cInf, mbalerr, qrunoff, vW, mW = emOut
    
    
    # For the Kragge we need to estimate the level in the drainage system.

    # first expand measured outflow to time series...
    drainOutF = pd.DataFrame(trange_sim,columns=['DateTime'])
    drainOutF.set_index('DateTime',inplace=True)

    #initially we assume that there is no storage in drainagesystem
    drainOutF['drainRate'] = qDr #np.concatenate([[0],qDr])

    tVals = meas_data_flow.index.astype(np.int64)/(1e9 * 3600 * 24)

    # outflow data is very noisy, better to take 14 day cumulative and work with thes
    idx = np.where((drainOutF.index > tVals.min()) & (drainOutF.index <= tVals.max()))
    outF_resampled1 = meas_data_flow['totF'].resample('1D').asfreq()
    outF_resampled2 = outF_resampled1.resample('1D').asfreq()
    outF_resampled = outF_resampled2.interpolate()

    outFRate = outF_resampled.diff()
    drainOutF['drainRate'].iloc[idx[0]] = outFRate.values[1:]

    # we use a shape factor to relate drainlevel to volume    
    # drainLev = volume/shapefact 
    # where shapefact = a Volume + b
    #
    # therefore Vol = b*drainLev/(1-a*drainLev)
    
    drainVolIni = wBdrVIni
    
    #Volume change in drainage system
    dQ = (qDr-drainOutF['drainRate'].values)
        
    drainVolN = np.cumsum(np.concatenate([[drainVolIni],dQ]))
    shape_fact = wB.drApar * (drainVolN >= wB.drCpar) + (wB.drBpar + (wB.drApar-wB.drBpar)/wB.drCpar * drainVolN) * (drainVolN < wB.drCpar) 
        
    drainLev = drainVolN/shape_fact

    # Add drainLev data to output list
    emOut = cumqDr, cDr, qDr, qInf, th, m, cInf, mbalerr, qrunoff, vW, mW, \
        drainLev, drainOutF, drainVolN

    lev_data_meas = meas_data_flow['level_data'].iloc[tdseries.tmeascal_ind].values
    like_lev_data = stats.norm(lev_data_meas, scale=stdMeas_levData)
    
    # Calculate likelihood: drainage level and outflow concentration
    cDr_meas = lab_data['Cl'].iloc[tdseries.tcalmeas_lab_ind].values
    like_cDr_data = stats.norm(loc=cDr_meas, scale=stdMeas_cDr)

    drainLev_cal = drainLev[tdseries.tcalib_ind[tdseries.tmeascal_ind]]
    cDr_cal = cDr[tdseries.tcalib_ind[tdseries.tmeascal_lab_ind]]

    est1 = like_lev_data.logpdf(drainLev_cal)
    est2 = like_cDr_data.logpdf(cDr_cal)

    est1[np.isnan(est1)] = 0
    est2[np.isnan(est2)] = 0

    #const = -0.5*np.log(2*np.pi*stdMeas_outF**2)
    #tmp = -(np.diff(outF)-np.diff(outF_sim))**2 /(2*stdMeas_outF**2)
    #est2 = const + tmp
    #lnp2 = np.sum(est2)

    logp_totFOut = np.sum(est1) + np.sum(est2)
    logp = logp_totFOut

    if np.isnan(logp):
        logp = -np.inf

    lFPar = {'cL': cL,
             'wB': wB,
             'lF': lF,
             'iniSt': initialStates}

    return logp, emOut, lFPar, emLF


def likelihood_Emission_LF (par, tdseries, lF, meteo_data,
                            meas_data_flow, lab_data, sdData):
    # Using parameter values in par we initialize the model class
    dzCL = par[0]
    cropFact = par[1]
    cLPor = par[2]
    cLthRes = par[3]*cLPor
    cLKSat = 10**par[4]
    cLbPar = par[5]

    wBtau1 = par[6]
    wBsig1 = 10**par[7]
    wBtau2 = wBtau1 + par[8]
    wBsig2 = 10**par[9]
    wBmfrac = par[10]
    wBvWIni = par[11] * (lF.height - dzCL)
    wBvWRes = par[12] * wBvWIni
    wBbFlow = 10**par[13]
    cLcIni = 10**par[14]
    wBcIni = 10**par[15]
    wBvWK = par[16]*(lF.height-dzCL)
    wBalphavW = 10**par[17]
    wBndEx = par[18]
    
    cLthIni = cLPor
    stdMeas_outF = sdData[0]
    stdMeas_cDr = sdData[1]

    coverLayer = {'dzCL': dzCL,
                  'cropFact': cropFact,
                  'Por': cLPor,
                  'thRes': cLthRes,
                  'KSat': cLKSat,
                  'bPar': cLbPar,
                  'cRain': 0.0,
                  'alphaRO': 0}
    
    cL = pd.Series(coverLayer)


    # wasteBody = {'tau1': wBtau1,
    #           'sig1': wBsig1,
    #           'tau2': wBtau2,
    #           'sig2': wBsig2,
    #           'mfrac': wBmfrac,
    #           'bFlow': wBbFlow,
    #           'vWRes': wBvWRes,
    #           'vWK': wBvWK,
    #           'alphavW': wBalphavW,
    #           'ndEx': wBndEx
    #           }
    wasteBody = {'tau1': wBtau1,
              'sig1': wBsig1,
              'tau2': wBtau2,
              'sig2': wBsig2,
              'mfrac': wBmfrac,
              'bFlow': wBbFlow,
              'vWRes': wBvWRes,
              'vWK': wBvWK,
              'alphavW': wBalphavW,
              'ndEx': wBndEx,
              'drApar': 0,
              'drBpar': 0,
              'drCpar': 0,
              'tau3': 1000,
              'sig3': 1,
              'tau4': 1000,
              'sig4': 1,
              'mfrac2': 0.5
              }
    wB = pd.Series(wasteBody)


    #Calculate initial volume in drainage system. Please note that the bottom liner is installed at a slope

    initialStates = {'thIni': cLthIni,
                     'cumQIni': meas_data_flow['totF'].iloc[0],
                     'cLcIni': cLcIni,
                     'wBvWIni': wBvWIni,
                     'wBcIni': wBcIni
                     }

    iniSt = pd.Series(initialStates)

    trange_sim = tdseries.trange.astype(np.int64)/(1e9 * 3600 * 24)

    # No infiltration (we pass an empty list)   
    infilF = []
    
    emLF = Emission_LF(trange_sim.values, iniSt, cL, wB, lF, meteo_data,
                         infilF)
    # Call model: this function solves the problem...
    emOut = emLF.SimulateEmissionLF()
    
    cumqDr, cDr, qDr, qInf, th, m, cInf, mbalerr, qrunoff, vW, mW = emOut
   
    # Calculate likelihood for simulated data matching measured data
    # We want to fit the data to the weekly outflow data, therefore we take
    # the difference (in the input we selected a 7 day cumulative total)

    outF = meas_data_flow['totF'].iloc[tdseries.tcalmeas_ind].values
    like_outF_data = stats.norm(loc=np.diff(outF), scale=stdMeas_outF)

    outF_sim = cumqDr[tdseries.tcalib_ind[tdseries.tmeascal_ind]]
    outF_sim = outF_sim - outF_sim[0]


    cDr_meas = lab_data['Cl'].iloc[tdseries.tcalmeas_lab_ind].values
    like_cDr_data = stats.norm(loc=cDr_meas, scale=stdMeas_cDr)

    cDr_sim = cDr[tdseries.tcalib_ind[tdseries.tmeascal_lab_ind]]

    est1 = like_outF_data.logpdf(np.diff(outF_sim))
    est2 = like_cDr_data.logpdf(cDr_sim)

    est1[np.isnan(est1)] = 0
    est2[np.isnan(est2)] = 0

    #const = -0.5*np.log(2*np.pi*stdMeas_outF**2)
    #tmp = -(np.diff(outF)-np.diff(outF_sim))**2 /(2*stdMeas_outF**2)
    #est2 = const + tmp
    #lnp2 = np.sum(est2)

    logp_totFOut = np.sum(est1) + np.sum(est2)
    logp = logp_totFOut

    if np.isnan(logp):
        logp = -np.inf

    lFPar = {'cL': cL,
             'wB': wB,
             'lF': lF,
             'iniSt': initialStates}

    return logp, emOut, lFPar, emLF