#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################################################/
 # @section LICENSE
 # The reservoir daily evaporation model (Version 1.0)
 # Copyright (C) 2020 The Surface Hydrology Group, Department of Civil
 # and Environmental Engineering, Texas A&M University.
 #
 # @section REFERENCES
 # If you use the code or a modified version of it in a publication, please cite -
 # [1] Gang Zhao, and Huilin Gao. "Estimating reservoir evaporation losses for the
 # United States: Fusing remote sensing and modeling approaches." Remote Sensing of
 # Environment 226 (2019):109-124.
 # [2] “Monitoring Daily Reservoir Evaporation Losses for the State of Texas”, Bing
 # jie Zhao, et al… (tentative authors from TAMU, DRI, TWDB, COE, and LCRA), in
 # preparation.
 #

 # @section NOTIFICATIONS
 # 1) These codes will be further examined and may be updated/revised at a later stage.
 # 2) Permission from Dr. Huilin Gao (hgao@tamu.edu) is needed for using the codes   # beyond the ongoing Texas Gridded Daily Reservoir project.
 #  3) Commercial usage is not allowed.
 #
 # @section DISCLAIMER
 # Please refer to the “Readme” file for the disclaimer of these codes.
 #############################################################################/

from equilibrium import equilibrim
import pandas as pd

def read_meteo(filename="powell_meteo.csv"):
    """
    Input the meteorological forcings file
    
    *make sure the unit of air temperature is degree Celsius
    *make sure the wind is measured at 2 meters
    
    Output the forcings dataframe
    """
    df_meteo=pd.read_csv(filename)
    df_meteo.index=pd.to_datetime(df_meteo.date)
    df_meteo.tem=df_meteo.tem-273.15
    ##select period
    mask=(df_meteo.index >= '2016-01-01') & (df_meteo.index <= '2021-12-31')
    df_meteo=df_meteo[mask]
    return df_meteo

def read_attributes(area_file="powell_area_depth.csv",fetch_file="powell_fetch.csv"):
    """
    Input the reservoir attribute file
    
    Output the reservoir attribute dataframe
    """
    area_df=pd.read_csv(area_file)
    area_df['date'] = pd.to_datetime(area_df['date'], format='%Y%m%d')
    area_df.set_index('date', inplace=True)
    ##select period
    mask=(area_df.index >= '2016-01-01') & (area_df.index <= '2021-12-31')
    area_df=area_df[mask]
    
    fetch_df=pd.read_csv(fetch_file)
    fetch_list=list(fetch_df.iloc[0])[1:]
    return area_df,fetch_list

def run_model():
    ##read_inputs
    meteo_d=read_meteo()
    area_d,fetch_l=read_attributes()
    
    ##calculate the daily fetch
    wind_dir=meteo_d['dir']
    wind_dir=list(((wind_dir+15)/30).astype('int'))
    fetch=[]
    for wd in (wind_dir):
        fetch.append(fetch_l[wd])
    
    ##set some constant, no necessary change
    longwave=-9999
    lat=36
    elev=1600
    
    Date=meteo_d.date
    month=area_d.index.month
    ##calculation daily evaporation
    tw=0
    ##Output file name
    f=open('powell_DLEM_output1.csv','w')
    f.write('date'+','+'evap_hs'+','+'error'+'\n')
    for i in range(len(Date)):
        if i ==0:
            tw0=meteo_d.tem[0]
        else:
            tw0=tw 
        tw,evap_hs,evap_nohs,error_=equilibrim(lat,area_d.depth_m[i],area_d.area_m2[i],elev,meteo_d.srd[i],longwave,meteo_d.tem[i],meteo_d.vpd[i],meteo_d.wind[i],tw0,fetch[i],month[i])
        f.write(Date[i]+','+str(evap_hs)+','+str(error_)+'\n')
    return None

def main():
    run_model()
    return None