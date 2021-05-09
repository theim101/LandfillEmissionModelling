#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 2018
Modules for extracting data from the Chronos database

@author: T.J. Heimovaara
"""

import sqlalchemy as sa
import pandas as pd
import numpy as np
import scipy as sp
import os
import subprocess, shlex
from collections import namedtuple
import seaborn as sns
import matplotlib.pyplot as plt

conn_str = 'mysql+pymysql://Guest_01:TUDelft01@localhost:3307/database_iDS03'

def wastebodyPropertiesWM(cellIdx):

    lFCode = np.array(['VP-06','VP-05a'])
    baseArea = np.array([28355,21343]) #m2
    topArea = np.array([9100,16000]) #m2
    slopeArea = baseArea-topArea
    slopeWidth = [38,42]
    slopeLength = slopeArea/slopeWidth
    lFHeight = 12+1.5
    surfArea = slopeLength *slopeWidth + topArea
    lFVolume = lFHeight * (0.5*slopeLength *slopeWidth + topArea)
    totVolume = sum(lFVolume)
    totWetWeight = [281083e3, 268407e3] #kg
    wetWeight = totWetWeight*lFVolume/totVolume;
    wetDensity = wetWeight/lFVolume;

    lFProp = namedtuple('lFProp', ['code', 'volume', 'wetWeight', 'wetDensity',
                                 'surfArea', 'baseArea', 'height'])
    lFProp = {'code': lFCode[cellIdx],
              'volume': lFVolume[cellIdx],
              'wetWeight': wetWeight[cellIdx],
              'wetDensity': wetDensity[cellIdx],
              'surfArea': surfArea[cellIdx],
              'baseArea': baseArea[cellIdx],
              'height': lFHeight}
    lF = pd.Series(lFProp)
    return lF

def wastebodyPropertiesBB(cellIdx):

    lFCode = np.array(['PP-11N','PP-11Z','PP-12'])
    baseArea = np.array([42084,39374,37963]) #m2
    topArea = np.array([20778,24365,16269]) #m2
    slopeArea = baseArea-topArea
    slopeWidth = 45
    slopeLength = slopeArea/slopeWidth
    lFHeight = 15+1.5
    surfArea = slopeLength *slopeWidth + topArea
    lFVolume = lFHeight * (0.5*slopeLength *slopeWidth + topArea)
    totVolume = sum(lFVolume)
    totWetWeight = 1216723e3 #kg
    wetWeight = totWetWeight*lFVolume/totVolume
    wetDensity = wetWeight/lFVolume

    lFProp = {'code': lFCode[cellIdx],
          'volume': lFVolume[cellIdx],
          'wetWeight': wetWeight[cellIdx],
          'wetDensity': wetDensity[cellIdx],
          'surfArea': surfArea[cellIdx],
          'baseArea': baseArea[cellIdx],
          'height': lFHeight
        }
    lF = pd.Series(lFProp)
    return lF

def wastebodyPropertiesKR(cellIdx):

    lFCode = np.array(['KR3','Kr4'])
    topArea = np.array([49700,59750]) #m2
    baseArea = np.array([56750, 59750]) #m2
    #baseLevel = np.array([5.04, 5.01]) # m
    wasteVolume = np.array([851561,622421]) #m3 obtained from Attero

    #slopeArea = baseArea-topArea
    #slopeWidth = 60
    #slopeLength = slopeArea/slopeWidth
    #lFHeight = wasteVolume/baseArea + 0.5
    #surfArea = slopeLength *slopeWidth + topArea
    lFHeight = wasteVolume / baseArea
    lFVolume = wasteVolume + 0.5*topArea  # assume a coverlayer of 50 cm

    wetWeight = wasteVolume*1300 + 0.5*topArea*1800 #kg wet waste density of 1300 kg/m3 and wet soil density of 1800 (20 % moisture)
    wetDensity = wetWeight/lFVolume

    drainWidth = [159, 147] # [m]width of base liner in order to calculate the storage volume in the drainage system
    drainLength = np.array([0,24,82,96,215,344]) # cumulative length of drain sections in drainage system
    drainHeight = np.array([5.15, 5.85, 6.20, 6.35, 6.85, 6.85]) # cumulative height of drain sections in drainage system (with respect to NAP)
    drainPor = 0.3 # porosity in drainage system
    # include rectangles for earlier sections...
   


    lFProp = {'code': lFCode[cellIdx],
          'volume': lFVolume[cellIdx],
          'wetWeight': wetWeight[cellIdx],
          'wetDensity': wetDensity[cellIdx],
          'surfArea': topArea[cellIdx],
          'baseArea': baseArea[cellIdx],
          'height': lFHeight[cellIdx],
          'drainWidth': drainWidth[cellIdx],
          'drainLength': drainLength,
          'drainHeight': drainHeight,
          'drainPor': drainPor,
          }
    lF = pd.Series(lFProp)
    return lF


def download_flow_level (pump_code,pklfile):
    # function checks if pklfile[0] exists, if not it will run a query on Chronos and download the relevant data and
    # save the data to the pklfile (in the format defined by the calling function)


    if pump_code =='PP-11N':
        sensor_lst = '2, 11'
    elif pump_code =='PP-11Z':
        sensor_lst = '3, 12'
    elif pump_code =='PP-12':
        sensor_lst = '1, 10'
    elif pump_code in ['VP-06','IDS-put 6']:
        sensor_lst = '4, 13'
    elif pump_code =='CF013':
        sensor_lst = '95, 96'
    elif pump_code =='CF014':
        sensor_lst = '100, 108'


    # check if pickel file exists
    if not (os.path.exists(pklfile)):
        # Need to query the data base for online sensor and manual gas measurements
        # Open a connection
        db_engine = sa.create_engine(conn_str)

        ## Read the inline measurements from the database
        # define the query, we need to extract data from the tables: fieldsensor,
        # sensormeas, component and unit. We explicitly select the sensor numbers
        # of the inline sensors. Please note the list of sensor-
        # numbers can be found by querying the database...
        sql_query = ("select fs.sensorno, fs.sensorname, sm.datetime, sm.componentno, " \
                     "c.name_nl as cname, sm.rawvalue, sm.value, u.name_nl as uname " \
                     "from fieldsensor as fs, sensormeas as sm, component as c, "
                     "unit as u where fs.sensorno = sm.sensorno and " \
                     "sm.unitno = u.unitno and sm.componentno=c.componentno " \
                     "and fs.sensorno in (" + sensor_lst + ")")

        # import the data in to the dataframe df_inline
        df_inline = pd.read_sql_query(sql_query, db_engine)
        # save df_inline to a pickle file
        df_inline.to_pickle(pklfile, compression='infer')
    else:
        # if the first file exists, we assume the second does as well
        df_inline = pd.read_pickle(pklfile)

    return df_inline

def download_sens_data_Kragge (pump_code,tmeas,pklfile):
    # function checks if pklfile[0] exists, if not it will run a query on Chronos and download the relevant data and
    # save the data to the pklfile (in the format defined by the calling function)

     # check if pickel file exists
    if not (os.path.exists(pklfile)):
        db_engine = sa.create_engine(conn_str)

        if pump_code == 'CF013':
            flow_sensnm = "'CF013','CF_013_o'"
            level_sensnm = "'CL012','CL_012_o'"
        elif pump_code == 'CF014':
            flow_sensnm = "'CF014','CF_014_o'"
            level_sensnm = "'CL014','CL_014_o'"
        elif pump_code == 'CF015':
            flow_sensnm = 'CF015'
            level_sensnm = 'CL008'



        infil_sensnm = 'CF001'

        # select the data for the different sensors per group
        # Cumulative Flow data
        sql_query = ("select fs.sensorno, fs.sensorname, sm.datetime, sm.componentno, "\
                         "c.name_nl as cname, sm.value, u.name_nl as uname "\
                         "from fieldsensor as fs, sensormeas as sm, component as c, "\
                         "unit as u where fs.sensorno = sm.sensorno and "\
                         "sm.unitno = u.unitno and sm.componentno=c.componentno "\
                         "and fs.sensorname in (" + flow_sensnm + ") and sm.componentno = 2701 "\
                         "and sm.datetime >= '2012-01-01 00:00:00'")

        #import the data in to a dataframe
        df_cumflow = pd.read_sql_query(sql_query,db_engine)

        #correct the values of the old time series as they are in 0.01*m3
        df_cumflow['value'].loc[(df_cumflow['sensorname']=='CF_013_o')] *= 0.01
        df_cumflow['value'].loc[(df_cumflow['sensorname']=='CF_014_o')] *= 0.01

        # Instantaneous Flow data
        # sql_query = ("select fs.sensorno, fs.sensorname, sm.datetime, sm.componentno, "\
        #                 "c.name_nl as cname, sm.value, u.name_nl as uname "\
        #                 "from fieldsensor as fs, sensormeas as sm, component as c, "\
        #                 "unit as u where fs.sensorno = sm.sensorno and "\
        #                 "sm.unitno = u.unitno and sm.componentno=c.componentno "\
        #                 "and fs.sensorname in (" + flow_sensnm + ") and sm.componentno = 2728 "\
        #                 "and sm.datetime >= '2018-01-01 00:00:00'")

        # #import the data in to a dataframe
        # df_flow = pd.read_sql_query(sql_query,db_engine)

        # Level data
        sql_query = ("select fs.sensorno, fs.sensorname, sm.datetime, sm.componentno, "\
                        "c.name_nl as cname, sm.value, u.name_nl as uname "\
                        "from fieldsensor as fs, sensormeas as sm, component as c, "\
                        "unit as u where fs.sensorno = sm.sensorno and "\
                        "sm.unitno = u.unitno and sm.componentno=c.componentno "\
                        "and fs.sensorname in (" + level_sensnm + ") "\
                        "and sm.datetime >= '2012-01-01 00:00:00'")

        df_level = pd.read_sql_query(sql_query,db_engine)


        sql_query = ("select fs.sensorno, fs.sensorname, sm.datetime, sm.componentno, "\
                        "c.name_nl as cname, sm.value, u.name_nl as uname "\
                        "from fieldsensor as fs, sensormeas as sm, component as c, "\
                        "unit as u where fs.sensorno = sm.sensorno and "\
                        "sm.unitno = u.unitno and sm.componentno=c.componentno "\
                        "and fs.sensorname in ('" + infil_sensnm + "') and sm.componentno = 2701 "\
                        "and sm.datetime >= '2012-01-01 00:00:00'")
        #import the data in to the dataframe df_inline
        df_cum_infil = pd.read_sql_query(sql_query,db_engine)
        ## create pivot tables for easy handling

        totF = df_cumflow[['datetime','value']].set_index('datetime')

        # flow_data = pd.pivot_table(df_flow, values='value', index=['datetime'],
        #                     columns=['sensorname'], aggfunc=np.sum)

        #clean up data (there are gaps in the data, the data is at a much higher frequency than required for our analysis...)

        level_data = df_level[['datetime','value']].set_index('datetime')

        cum_infil_data = df_cum_infil[['datetime','value']].set_index('datetime')


        finter = sp.interpolate.interp1d(totF.index.astype(np.int),totF.value,fill_value='extrapolate')
        totF_val = finter(tmeas.astype(np.int))
        totalFlow = pd.DataFrame(data = totF_val, index=tmeas)

        finter = sp.interpolate.interp1d(level_data.index.astype(np.int),level_data.value,fill_value='extrapolate')
        level_val = finter(tmeas.astype(np.int))
        levelD = pd.DataFrame(data = level_val, index=tmeas)

        finter = sp.interpolate.interp1d(cum_infil_data.index.astype(np.int),cum_infil_data.value,fill_value='extrapolate')
        cinfil_val = finter(tmeas.astype(np.int))
        totInfilF = pd.DataFrame(data = cinfil_val, index=tmeas)

        df_inline = totalFlow.rename(columns={0:'totalFlow'})
        df_inline['levelD'] = levelD
        df_inline['totInfilF'] = totInfilF

        #remove missing values from data
        df_inline['levelD'].loc['2017-09-12':'2018-05-16']=np.nan
        df_inline['totalFlow'].loc['2017-09-12':'2018-01-01']=np.nan

        # save df_inline to a pickle file
        df_inline.to_pickle(pklfile, compression='infer')
    else:
        # if the first file exists, we assume the second does as well
        df_inline = pd.read_pickle(pklfile)

    return df_inline




def download_lab (pump_code,pklfile):
    # function checks if pklfile[0] exists, if not it will run a query on Chronos and download the relevant data and
    # save the data to the pklfile (in the format defined by the calling function)

    # components: Cl:508, NH4: 290, Br:391, TOC: 1574
    clist = '508, 290, 391, 1574'

    if pump_code =='PP-11N':
        measpointno = '3033'
    elif pump_code =='PP-11Z':
        measpointno = '3034'
    elif pump_code =='PP-12':
        measpointno = '3032'
    elif pump_code == 'VP-06':
        measpointno = '9638'
    elif pump_code == 'CF013':
        measpointno = '9726'
    elif pump_code == 'CF014':
        measpointno = '9728'

    # check if pickel file exists
    if not (os.path.exists(pklfile)):
        # Need to query the data base for online sensor and manual gas measurements
        # Open a connection
        db_engine = sa.create_engine(conn_str)

        ## Read the inline measurements from the database
        # define the query, we need to extract data from the tables: fieldsensor,
        # sensormeas, component and unit. We explicitly select the sensor numbers
        # of the inline sensors. Please note the list of sensor-
        # numbers can be found by querying the database...
        sql_query = ("select s.date, s.measpointno, s.measpointname, "\
            "a.componentno, a.value, c.name_nl as cname, u.name_nl as uname "\
            "from sample as s, analysis as a, component as c, unit as u "\
            "where s.sampleno = a.sampleno and a.componentno = c.componentno "\
            "and a.unitno = u.unitno "\
            "and s.measpointno = " + measpointno \
            + " and a.componentno in (" + clist +");")

        # import the data in to the dataframe df_inline
        df_lab = pd.read_sql_query(sql_query, db_engine)
        # save df_inline to a pickle file
        df_lab.to_pickle(pklfile, compression='infer')
    else:
        # if the first file exists, we assume the second does as well
        df_lab = pd.read_pickle(pklfile)

    return df_lab


def download_meteoKNMI (t_range, weather_station, pklfile):
    # This function runs a command line argument to download data from
    # the KNMI website. The meteo station code is given in weather_station
    # the time range for the data is given in the tuple t_range

    # two meteo stations:
    # 249:         4.979       52.644      -2.50  BERKHOUT
    # 269:         5.526       52.458      -4.00  LELYSTAD

    # 340: Woensdrecht
    # 350: Gilze-Rijen
    # download rainfall data in t_range to 0DailyDataRLT.txt
    #http://projects.knmi.nl/klimatologie/daggegevens/getdata_dag.cgi
    #!wget -O 0DailyDataRLT.txt --post-data="stns=249:269&vars=ALL&start=20000101&end20180915" http://projects.knmi.nl/klimatologie/daggegevens/getdata_dag.cgi
    # check if pickel file exists
    if not (os.path.exists(pklfile)):
        my_command = 'wget -O 0DailyDataRLT.txt --post-data="stns=' + weather_station \
            + '&vars=ALL&start=' + t_range[0] + '&end=' + t_range[1] \
            + '"  http://projects.knmi.nl/klimatologie/daggegevens/getdata_dag.cgi'

        my_args = shlex.split(my_command)

        subprocess.run(my_args)  # doesn't capture output

        m_dat = pd.read_csv('0DailyDataRLT.txt', sep=',', header=None, engine='c',
                   comment='#', parse_dates=True, infer_datetime_format=True,
                   names=['STN', 'YYYYMMDD', 'DDVEC', 'FHVEC', 'FG', 'FHX', 'FHXH',
                          'FHN', 'FHNH', 'FXX','FXXH', 'TG', 'TN', 'TNH', 'TX',
                          'TXH', 'T10N', 'T10NH', 'SQ', 'SP', 'Q', 'DR', 'RH',
                          'RHX', 'RHXH', 'EV24', 'PG', 'PX', 'PXH', 'PN', 'PNH',
                          'VVN', 'VVNH', 'VVX', 'VVXH', 'NG', 'UG','UX','UXH',
                          'UN', 'UNH'])

        m_out = pd.DataFrame()
        m_out['datetime'] = pd.to_datetime(m_dat['YYYYMMDD'],format='%Y%m%d')
        m_out['rain'] = m_dat['RH']
        m_out['rain'] = m_out['rain'].replace(-1,0.25)
        m_out['rain'] = m_out['rain']/1e4
        m_out['pEV'] = m_dat['EV24']/1e4
        m_out['temp'] = m_dat['TG']/10
        m_out.set_index('datetime',inplace=True)

        # we only require rainfall data and potential evaporation
        # we also want change the YYYYMMDD column to a pandas datetime

        m_out.to_pickle(pklfile,compression='infer')
    else:
        m_out = pd.read_pickle(pklfile)
    return(m_out)

def remove_outliers_inline(inline_par):
    #inline_par.index = pd.to_datetime(inline_par.index)
    # resample from totalflow to daily data
    tmpF = inline_par.dropna(subset=['totalflow']).copy()
    #delt = np.diff(tmpF.index).astype('float64')/(1e9 * 3600)
    tmpF['diffFlow']=tmpF['totalflow'].diff()
    tmpF['flowR'] = tmpF['diffFlow']
    tmpF['flowR'].iloc[0] = 0
    #tmpF[['totalflow','flowR','level']].plot(subplots=True,style='.-',grid=True)

    max_flow_ini = np.amax(tmpF['flowR'].iloc[0:20000])
    # Identify outliers using z-score
    zscore = (tmpF['flowR'] - tmpF['flowR'].mean())/tmpF['flowR'].std(ddof=0)
    tmpF['outliers_zscore'] = np.abs(zscore) > 3

    outlieridx = np.where(tmpF.outliers_zscore==True)[0]
    rngidx = np.where(tmpF['flowR'].iloc[outlieridx[1]:outlieridx[1]+10000]> 5*max_flow_ini)
    lastidx = outlieridx[1] + rngidx[0][-1]

    # slope of flowR between outlieridx[0] and outlieridx[1] is 50 times to high
    # slope of flowR between outlieridx[1] and lastidx is 10 times to high

    tmpF['cFlowR'] = tmpF['flowR']
    tmpF['cFlowR'].iloc[outlieridx[0]:outlieridx[1]] = tmpF['cFlowR'].iloc[outlieridx[0]:outlieridx[1]]/50
    tmpF['cFlowR'].iloc[outlieridx[1]:lastidx] = tmpF['cFlowR'].iloc[outlieridx[1]:lastidx]/10
    tmpF['cFlowR'].iloc[outlieridx[0]] = tmpF['cFlowR'].iloc[outlieridx[0]-1]
    tmpF['cFlowR'].iloc[outlieridx[1]] = tmpF['cFlowR'].iloc[outlieridx[1]-1]
    tmpF['totF'] = tmpF['cFlowR'].cumsum()

    totF = tmpF['totF']

    return totF

def remove_outliers_inlineBB(totF0):

    # %% Correct totalflow (remove drops in data set):
    # replace nan in totflow column by interpolated values
    totF1 = totF0[slice('2012-06-01','2021-03-01')].interpolate()
    #totF1 = totF0
    # difference in level
    dlev = totF1['level'].diff()

    # cumulative drop in level...
    dlev2 = np.cumsum(dlev*(dlev<0))

    # Step 1
    # Find the first negative flow,
    # then correct the value based on the cumulative level change
    # which has been calibrated using the level change from a
    # time period of 4 days before the point of correction.
    # We calibrate the level change by relating the total level change
    # over 4 days with the cumulative pumped volume in these four days...
    # correction is carried out on all subsequent data...

    totF2 = totF1.copy()
    # find the negative values in the gradient of the cumulative flow data set
    idx1 = np.where(totF1['totalflow'].diff()<0)[0]

    # create time series of four days before idx1
    for ii in np.arange(len(idx1)):
        irange = np.arange(idx1[ii]-24*4*4,idx1[ii])
        cumQ = totF2['totalflow'].iloc[irange[-20]]-totF2['totalflow'].iloc[irange[0]]
        dLevtot = dlev2[irange[-20]]-dlev2[irange[0]]
        effSurf = -cumQ/dLevtot

        # coordinate with negative diff value: idx1[0]
        # replace value with
        dQFlow = totF2['totalflow'].iloc[idx1[ii]-1] - \
            totF2['totalflow'].iloc[idx1[ii]] + \
            (dlev2[idx1[ii]-1]-dlev2[idx1[ii]])*effSurf

        # add correction from dlev2: (dlev2[idx1[0]]-dlev2[idx1[0]])*effSurf

        totF2['totalflow'].iloc[idx1[ii]:] += dQFlow

    # Step 2

    totF3 = totF2.copy()
    # find the negative values in the gradient of the cumulative flow data set
    idx2 = np.where(totF3['totalflow'].diff()>1000*totF3['totalflow'].diff().mean())[0]

    for ii in np.arange(len(idx2)):
        irange = np.arange(idx2[ii]-24*4*4,idx2[ii])
        cumQ = totF3['totalflow'].iloc[irange[-20]]-totF3['totalflow'].iloc[irange[0]]
        dLevtot = dlev2[irange[-20]]-dlev2[irange[0]]
        effSurf = -cumQ/dLevtot

        # coordinate with negative diff value: idx1[0]
        # replace value with
        dQFlow = totF3['totalflow'].iloc[idx2[ii]-1] - \
            totF3['totalflow'].iloc[idx2[ii]] + \
            (dlev2[idx2[ii]-1]-dlev2[idx2[ii]])*effSurf

        # add correction from dlev2: (dlev2[idx1[0]]-dlev2[idx1[0]])*effSurf

        totF3['totalflow'].iloc[idx2[ii]:] += dQFlow


    return totF3['totalflow']

def download_meteoKNMI_etmgeg (t_range, weather_station, pklfile, inpfile):
    # This function runs a command line argument to download data from
    # the KNMI website. The meteo station code is given in weather_station
    # the time range for the data is given in the tuple t_range

    # two meteo stations:
    # 249:         4.979       52.644      -2.50  BERKHOUT
    # 269:         5.526       52.458      -4.00  LELYSTAD

    # 340: Woensdrecht
    # 350: Gilze-Rijen
    # download rainfall data in t_range to 0DailyDataRLT.txt
    #http://projects.knmi.nl/klimatologie/daggegevens/getdata_dag.cgi
    #!wget -O 0DailyDataRLT.txt --post-data="stns=249:269&vars=ALL&start=20000101&end20180915" http://projects.knmi.nl/klimatologie/daggegevens/getdata_dag.cgi
    # check if pickel file exists
    if not (os.path.exists(pklfile)):
        m_dat = pd.read_csv(inpfile, sep=',', header=None, engine='c',
               na_values='     ', comment='#', parse_dates=True, 
               infer_datetime_format=True, skiprows=48,
               names=['STN', 'YYYYMMDD', 'DDVEC', 'FHVEC', 'FG', 'FHX', 'FHXH',  
                      'FHN', 'FHNH', 'FXX', 'FXXH', 'TG', 'TN', 'TNH', 'TX', 'TXH', 
                      'T10N', 'T10NH', 'SQ', 'SP', 'Q', 'DR', 'RH', 'RHX', 'RHXH', 
                      'PG', 'PX', 'PXH', 'PN', 'PNH', 'VVN', 'VVNH', 'VVX', 'VVXH',
                      'NG', 'UG', 'UX', 'UXH', 'UN', 'UNH', 'EV24'])

        m_out = pd.DataFrame()
        m_out['datetime'] = pd.to_datetime(m_dat['YYYYMMDD'],format='%Y%m%d')
        m_out['rain'] = m_dat['RH']
        m_out['rain'] = m_out['rain'].replace(-1,0.25)
        m_out['rain'] = m_out['rain']/1e4
        m_out['pEV'] = m_dat['EV24']/1e4
        m_out['temp'] = m_dat['TG']/10
        m_out.set_index('datetime',inplace=True)

        # we only require rainfall data and potential evaporation
        # we also want change the YYYYMMDD column to a pandas datetime

        m_out.to_pickle(pklfile,compression='infer')
    else:
        m_out = pd.read_pickle(pklfile)
    return(m_out)