"""
Created on Sat Oct 20 2018
Modules for extracting data from the Chronos database

@author: T.J. Heimovaara
"""

import numpy as np
import scipy as sp
import scipy.stats as stats
import pandas as pd

from context import dbl
# import DataBaseLibrary as dbl
#from pydream.core import run_dream
from pydream.parameters import SampledParam
#from pydream.convergence import Gelman_Rubin
#import inspect

# Meteorological data will be obtained from two sources:
# 1: a close by weather station (for WM Berkhout, for BB: Lelystad)
#    we will use the evapotranspiration data obtained from the weather station...
# 2: rainfall from the 1km gridded radar corrected interpolated rainfall data obtained
#    from climate4impact...

weather_station = '269' #Lelystadut
t_range = ['20030101','20210301']

pklfile = './DataFiles/meteoLS.gz'
#path = './MeteoData/WM_Rain_2008-2019.bz2'

inpfile = 'etmgeg_269.txt'
# Read data from close by meteo station
meteo_data_stat = dbl.download_meteoKNMI_etmgeg (t_range, weather_station, pklfile, inpfile)

#meteo_data_stat = dbl.download_meteoKNMI (t_range, weather_station, pklfile)
meteo_data_stat = meteo_data_stat.rename(columns={'rain':'rain_station'})

# Read data from gridded radar corrected interpolated rainfall data
#ain_radar = pd.read_pickle(fpath,compression='infer')
# transform rain values from kg/m2 (mm) to m water column
#ain_radar['rain'] = rain_radar['rain']/1e3
# Merge the station data and the interpolated rain data in to a single dataframe
meteo_data = meteo_data_stat
# meteo_data is top boundary condition. We run the model from 2003 onward
meteo_data = meteo_data[slice('2003-01-01','2021-03-01')]

#eteo_data.rain.loc[meteo_data['rain'].isnull()] = \
#   meteo_data.rain_station.loc[meteo_data['rain'].isnull()]

## Download flow and level data from CHRONOS
pump_code = 'PP-11Z'
pklfile = './DataFiles/flowdata_PP-11Z.gz'

df_inline = dbl.download_flow_level (pump_code, pklfile)
# We create a pivot table based on column cname (component names)
#inline_par = pd.pivot_table(df_inline, values='rawvalue', index=['datetime'],
#                      columns=['cname'], aggfunc=np.sum)
totF0 = pd.pivot_table(df_inline, values='rawvalue', index=['datetime'],
                      columns=['cname'], aggfunc=np.sum)
#totF = dbl.remove_outliers_inline(inline_par)
totF = dbl.remove_outliers_inlineBB(totF0)

# as the model allows for leachate recirculation it expects a totIniflF dataset
# For a situation where no leachate is recirculated, we set the totIniflF to zero

totF0['totInfilF'] = 0

levelD = totF0['level']
infilF = totF0['totInfilF']

#sensData = dbl.download_sens_data_Kragge (pump_code, tmeas, pklfile)

# We create a pivot table based on column cname (component names)
#inline_par = pd.pivot_table(df_inline, values='rawvalue', index=['datetime'],
#                      columns=['cname'], aggfunc=np.sum)
#totF = sensData['totalFlow']
#levelD = sensData['levelD']
#infilF = sensData['totInfilF']

# Download laboratory data for pump pit
pklfile = './DataFiles/labdata_PP-11Z.gz'

df_lab = dbl.download_lab(pump_code, pklfile)

# We create a pivot table based on column cname (component names)
lab_data = pd.pivot_table(df_lab, values='value', index=['date'],
                          columns=['cname'], aggfunc=np.sum)

lab_data = lab_data.rename(index=str, columns={'Ammonium (als N)': 'NH4',
                                   'Bromide': 'Br',
                                   'Chloride': 'Cl',
                                   'Totaal organisch koolstof (TOC)': 'TOC'})

lab_data.index = pd.to_datetime(lab_data.index)

# meteo_data is top boundary condition. We run the model over 10 years
meteo_data = meteo_data[slice('2003-01-01','2021-03-01')]


# Define simulation time range (trange)
trange = pd.date_range(start='2003-01-01',end = '2021-03-01',freq='D')

# Select measurements, should fall within trange.
# tmeas contains times where measurements are available!
# can have multiple tmeas vectors for different types of measurements
# totF contains measured data from mid 2012. We choose to start on the 2012-07-01
# Because the outflow is influenced by operator decisions we choose to select weekly
# cumulative totals...
measFreq = 7
tmeas = pd.date_range(start='2012-06-14',end = '2021-03-01',freq='7D')
finter = sp.interpolate.interp1d(totF.index.astype(int),totF.values)
totF_val = finter(tmeas.astype(int))
totF2 = pd.DataFrame(data = totF_val, index=tmeas)
totF2 = totF2-totF2.iloc[0]
meas_data_flow = totF2.rename(columns = {0:'totF'})

# Define calibration time range. This will be used by DREAM to compare
# simulated values with calibration set...
# Data set to be matched by modifying parameters...
tcalib = pd.date_range(start='2012-06-14',end = '2020-01-01',freq='D')

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
cellIdx = 1 # 11Z = 0, 11Z = 1, 12 = 2
lF = dbl.wastebodyPropertiesBB(cellIdx) #m2

## Run model wil optimal parameter set...

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
         'wBndEx': [200, 600, 500]
         }

par_df = pd.DataFrame(data=par_d)

# lower_limits = par_df.iloc[0].values
# range_vals = par_df.iloc[1].values - par_df.iloc[0].values
# #start_pos = np.array([-3.4, 0.2, 0.8, 1.5, 0, (0.3*1+0.5*lF.height)*lF.surfArea,
# #                      0, 5, 10])

# parameters_to_sample = SampledParam(stats.uniform, loc=lower_limits, scale=range_vals)
# sampled_parameter_names = [parameters_to_sample]

sdData = np.array([1,25])

