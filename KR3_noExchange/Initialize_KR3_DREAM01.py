"""
Created on Sat Oct 20 2018
Modules for extracting data from the Chronos database

@author: T.J. Heimovaara
"""

import numpy as np
import scipy as sp

import pandas as pd
import matplotlib.pyplot as plt
from context import dbl
# import DataBaseLibrary as dbl
#from pydream.core import run_dream
#from pydream.parameters import SampledParam
#from pydream.convergence import Gelman_Rubin
#import inspect

# Meteorological data will be obtained from two sources:
# 1: a close by weather station (for WM Berkhout, for BB: Lelystad)
#    we will use the evapotranspiration data obtained from the weather station...
# 2: rainfall from the 1km gridded radar corrected interpolated rainfall data obtained
#    from climate4impact...

# surface areas of Kragge compartment 3 and 4
topArea = np.array([57593,61012]) #m2
baseArea = np.array([56137, 58206]) # m2
baseLevel = np.array([5.04, 5.01]) # m
fCrop = 1 # optimized crop factor by Loys Vermeijden (2018)

weather_station = '350'
pklfile = 'meteoGilzeRijen.pkl'
t_range = ['20030101','20210101']


pklfile = './DataFiles/meteoKR.bz2'
#path = './MeteoData/WM_Rain_2008-2019.bz2'

# Read data from close by meteo station
meteo_data_stat = dbl.download_meteoKNMI (t_range, weather_station, pklfile)
meteo_data_stat = meteo_data_stat.rename(columns={'rain':'rain_station'})

# Read data from gridded radar corrected interpolated rainfall data
#ain_radar = pd.read_pickle(fpath,compression='infer')
# transform rain values from kg/m2 (mm) to m water column
#ain_radar['rain'] = rain_radar['rain']/1e3
# Merge the station data and the interpolated rain data in to a single dataframe
meteo_data = meteo_data_stat
# meteo_data is top boundary condition. We run the model from 2003 onward
meteo_data = meteo_data[slice('2003-01-01','2021-01-01')]

pump_code = 'CF013'
pklfile = './DataFiles/flowdata_CF013.bz2'

measFreq = 1
trange = pd.date_range(start='2003-01-01',end = '2020-12-31',freq='D')
tmeas = pd.date_range(start='2012-06-12',end = '2020-12-31',freq='D')


tcalib1 = pd.date_range(start='2012-06-12',end = '2017-03-03',freq='D')
tcalib2 = pd.date_range(start='2018-05-21',end = '2020-07-01',freq='D')
tcalib = tcalib1.append(tcalib2)
#tcalib = tcalib1

sensData = dbl.download_sens_data_Kragge (pump_code, tmeas, pklfile)

# We create a pivot table based on column cname (component names)
#inline_par = pd.pivot_table(df_inline, values='rawvalue', index=['datetime'],
#                      columns=['cname'], aggfunc=np.sum)
totF = sensData['totalFlow']
levelD = sensData['levelD']
infilF = sensData['totInfilF']

## Correction of the data due to datagap from sept 2017 to January 2018
# Assumption: leachate was trucked off site during 2017. From January 2018 
# leachate was buffered until treatement system started.
# volume trucked of site is estimated from historical trend from 2014-12-30 to
# 2017-09-05 which gives 30 m3/day removed leachate. This was over 119 days.
# The difference in drainage level between 2017-09-05 and 2018-05-22 leads to 
# a buffered volume of about 8043 m3. This was over 259 days.
#  Leachate infiltration started before level measurements, the volume 
# infilrated is 3179 m3

# volume estimated from change in level during measurement gap
bufferedlVol = 8043 #m3 leachate
# volume infiltrated in infiltration system until level measurements restarted
infiltratedVol = 3179 # m3 

#on average 30 m3/day
#259 days until level measurements restart
#119 dats until start measurement system 20180102

extra_qOut = 30 * 119 # extracted leachate
extra_storage = 30 * (259-119) + 3179

# add extra_qOut to totF after 2018-01-02
totF['2018-01-02':] = totF['2018-01-02':] + totF['2017-09-05'] + extra_qOut
totF = totF.interpolate()

levelD = levelD.interpolate()
#totF = dbl.remove_outliers_inline(inline_par)
#totF = dbl.remove_outliers_inlineBB(totF0)

## Download laboratory data for pump pit
pklfile = './DataFiles/labdata_CF013.bz2'
df_lab = dbl.download_lab(pump_code, pklfile)

# We create a pivot table based on column cname (component names)
lab_data_all = pd.pivot_table(df_lab, values='value', index=['date'],
                          columns=['cname'], aggfunc=np.sum)

lab_data_all = lab_data_all.rename(index=str, columns={'Ammonium (als N)': 'NH4',
                                   'Bromide': 'Br',
                                   'Chloride': 'Cl',
                                   'Totaal organisch koolstof (TOC)': 'TOC'})
lab_data = pd.DataFrame()
lab_data['Cl'] = lab_data_all['Cl'].dropna()

lab_data.index = pd.to_datetime(lab_data.index)



# Select measurements, should fall within trange.
# tmeas contains times where measurements are available!
# can have multiple tmeas vectors for different types of measurements
# totF contains measured data from mid 2012. We choose to start on the 2012-07-01
# Because the outflow is influenced by operator decisions we choose to select weekly
# cumulative totals...


totF2 = totF-totF.iloc[0]
meas_data_flow = pd.DataFrame(totF2).rename(columns={'totalFlow':'totF'})
meas_data_flow['level_data'] = levelD
# Define calibration time range. This will be used by DREAM to compare
# simulated values with calibration set...
# Data set to be matched by modifying parameters...


# In order to facilitate quick and easy comparison of simulation with data
# we need to define the overlapping indices:
# tmeas_ind: trange[tmeas_ind] = tmeas
# tcalib_ind: trange[tcalib_ind]=tcalib
# tcalmeas_ind, tmeascal_ind: tmeas[tcalmeas_ind]=tcalib[tmeascal_ind]


xy, ind1, tmeas_ind = np.intersect1d(tmeas,  trange,
                                     return_indices=True)
xy, ind1, tcalib_ind = np.intersect1d(tcalib, trange,
                                      return_indices=True)
xy, tmeascal_ind, tcalmeas_ind = np.intersect1d(tcalib, tmeas,
                                                return_indices=True)

xy, tlabmeas_ind, tmeaslab_ind = np.intersect1d(lab_data.index,  trange,
                                         return_indices=True)
xy, tmeascal_lab_ind, tcalmeas_lab_ind = np.intersect1d(tcalib, lab_data.index,
                                                        return_indices=True)


tdata = {'trange':trange,
         'tmeas':tmeas,
         'tcalib':tcalib,
         'tmeas_ind':tmeas_ind,
         'tcalib_ind':tcalib_ind,
         'tcalmeas_ind':tcalmeas_ind,
         'tmeascal_ind':tmeascal_ind,
         'tlabmeas_ind': tlabmeas_ind,
         'tmeaslab_ind': tmeaslab_ind,
         'tmeascal_lab_ind': tmeascal_lab_ind,
         'tcalmeas_lab_ind': tcalmeas_lab_ind}

tdseries = pd.Series(tdata)


# Obtain landfill specific properties
cellIdx = 0 # Cell 3 = 0, Cell 4 = 1
lF = dbl.wastebodyPropertiesKR(cellIdx) #m2


## Prepare DREAM model...
# Model parameters which are required to calculate fluxes etc. (often need to
# optimized).

# List of parameters

par_d = {'dzCL': [0.5, 1.9, 0.9267], #%m
         'cropFact': [0.75, 1.5, 0.9584], #[-]
         'cLPor': [0.15, 0.60, 0.1722], #[-]
         'cLthRes': [0.001, 0.9, 0.0293], # percentage of thTot
         'cLKSat': [-5, 1, -1.1936], # 10 log
         'cLbpar': [0, 8, 2.8698], #%10log! [m/d]
         'wBtau1': [1, 200, 37.5672],
         'wBsig1': [-5, 1, -0.4642],
         'wBtau2': [0, 5*365, 170.3241],
         'wBsig2': [-5, 1, 0.3246],
         'wBmfrac': [0.01, 1.0, 0.2338],
         'wBthIni': [0.05, 0.95, 0.5516],
         'wBfRes' : [0, 0.9, 0.0851],
         'wBbFlow': [-7,-2,-4.0935], #10log!
         'cLcIni': [-4, 3, -1.5343], #10log!
         'wBcIni': [2, 6, 3.1901],
         'wBvWK': [0, 0.9, 0.1546],
         'wBalphavW': [-9,2,-4.3908],#, #10log!
         'wBndEx': [365, 750, 500],
         'wBdrVIni':[0, 6, 5], #10log! initial leachate level in drainage system
         'wBdrApar': [0, 5, -0.3], # 10log! surface area for shape factor (when volume surpasses drainage volume)
         'wBdrBpar': [0, 5, 3.5 ], # 10log! surface factor (empirical fit) should optimize to be larger than A 
         'wBdrCpar': [0, 6, 3], # Max Volume Leachate in drainage layer
         'cLalphaRO': [0, 1, 0.5], # parameter indicating fraction of runoff when thExcess > 0
         'wBtau3' : [0, 10*365, 10],
         'wBsig3' : [-5, 2, 1],
         'wBtau4' : [0, 500*365, 3100],
         'wBsig4' : [-5, 2, 1],
         'wBmfrac2': [0.01, 1.0, 0.2338],
         }

par_df = pd.DataFrame(data=par_d)

# lower_limits = par_df.iloc[0].values
# range_vals = par_df.iloc[1].values - par_df.iloc[0].values
#start_pos = np.array([-3.4, 0.2, 0.8, 1.5, 0, (0.3*1+0.5*lF.height)*lF.surfArea,
#                      0, 5, 10])

sdData = np.array([0.37,700])

