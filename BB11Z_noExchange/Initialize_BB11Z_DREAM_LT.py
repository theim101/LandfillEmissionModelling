"""
Created on Sat Oct 20 2018
Modules for extracting data from the Chronos database

@author: T.J. Heimovaara
"""

import numpy as np
import pandas as pd

from Initialize_BB11Z_DREAM01 import *

# meteo_data is top boundary condition. We run the model over 10 years
meteo_data = meteo_data[slice('2003-01-01','2021-03-01')]
meteo_data_tmp = meteo_data[slice('2003-03-02','2021-03-01')]
md2 = meteo_data.append(meteo_data_tmp,ignore_index=True)
md2 = md2.append(meteo_data_tmp,ignore_index=True)
#md2 = md2.append(meteo_data_tmp,ignore_index=True)

trange = pd.date_range('2003-01-01',periods=md2.shape[0],freq='D')
md2 = md2.set_index(trange)
meteo_data = md2

# Define simulation time range (trange)
#trange = pd.date_range(start='2008-01-01',end = '2019-06-30',freq='D')


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


