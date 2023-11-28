 # ###########################################################################/
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

from math import log, exp, sqrt, pow
from math import sin, cos, tan, atan, pi
from constants import T_abs

def delcalc(ta):
    # FUNCTION TO CALCULATE THE SLOPE OF THE VAPOUR PRESSURE CURVE 
    # INPUT    TA - AIR TEMPERATURE (deg. C)
    # OUTPUT    DELCALC - SLOPE OF THE VAPOUR PRESSURE CURVE (kPa deg. C-1)
    
    ea = 0.6108*exp(17.27*ta/(ta+237.3))
    delcalc = 4098*ea/pow((ta+237.3),2.)
    return delcalc

def alambdat(ta):
    # FUNCTION TO CORRECT THE LATENT HEAT OF VAPORISATION FOR TEMPERATURE
    # INPUT:     T = TEMPERATURE (deg. C)
    # OUTPUT:    ALAMBDAT = LATENT HEAT OF VAPORISATION (MJ kg-1)
    return 2.501-ta*2.361e-3


def psyconst(p, alambda):   ## S2.9 
    # INPUT:    P = ATMOSPHERIC PRESSURE (kPa)    ALAMBDA = LATENT HEAT OF VAPORISATION (MJ kg-1)
    # OUTPUT:    PSYCONST = PSYCHROMETRIC CONStanT (kPa deg. C-1) 
    return 0.00163*p/alambda


def v_sat(ta):
    # ta in celcius, v_sat in kPa 
    return 0.6108*exp(17.27*ta/(ta+237.3))


def tvpd2wbt(ta, vpd): ##S2.2 S2.3
    # ta is celsius, vpd in kPa 
    # t_wb, vpsat, t_d #dew point temperature
    vpsat = v_sat(ta)
    t_d = (116.9+237.3*log(vpsat-vpd))/(16.78-log(vpsat-vpd))
    t_wb = (0.00066*100.*ta + 4098.*(vpsat-vpd)/pow(t_d+237.3,2)*t_d)/(0.00066*100.+4098.*(vpsat-vpd)/pow(t_d+237.3,2.))
    # t_wb in celsius 
    return t_wb


def wspd2m(wspd, hgt):
    # wspd is wind speed at height hgt (m)
    # http:##edis.ifas.ufl.edu/ae459 
    wpd2m = wspd*4.87/log(67.8*hgt-5.42)
    # wpd2m in m/s 
    return wpd2m


def airdens(ta, elev):
    # airds is air density (kg/m3) at elev (m) and ta (C) 
    # https:##en.wikipedia.org/wiki/Density_of_air#Altitude 
    # p, p0     sea level standard atmospheric pressure, 101.325 kPa 
    # g     earth-surface gravitational acceleration, 9.80665 m/s2 
    # r     ideal (universal) gas constant, 8.31447 J/(mol K) 
    # m     molar mass of dry air, 0.0289644 kg/mol 
    # t0, l     temperature lapse rate, 0.0065 K/m 
    p0 = 101.325
    g = 9.80665
    r = 8.31447
    m = 0.0289644
    l = 0.0065
    
    ta = ta+273.15     ## celsius to kelvin 
    t0 = ta+l*elev
    p = p0*pow(1-l*elev/t0,g*m/r/l)*1000.0     # in Pa 
    airds = p*m/r/ta
    airds = 1.225 if airds>1.225 else airds
    return airds

def atm_p(ta, elev):
    # ta (C) and elev (m) 
    atmp = 101.3*pow((T_abs + ta -0.0065*elev)/(T_abs + ta),5.26)  ## S2.10 kPa 
    return atmp

def air_emissivity(ta, ea, fcd):
    # ta (C) and ea (kPa) and fcd 
    # Satterlund 1979 and ...
    em_a = 1.08*(1.-exp(-pow(ea*10.,(ta+T_abs)/2016.)))*(1+0.22*pow(fcd,2.75))
    return em_a

def cloud_factor(K, M, lat, elev):
    # K is incoming shortwave radiation (MJ/m2/d) and M is the month 
    # J, omega, delta, dr
    # Kso, Ket, Kr, fcd
    # lat_r
    
    J = (int)(30.4*M-15)
    delta = 0.409*sin(2.*pi*J/365. - 1.39)
    lat_r = lat/180.*pi
    omega = pi/2. - atan(-tan(lat_r)*tan(delta)/sqrt(1-tan(lat_r)*tan(lat_r)*tan(delta)*tan(delta)))
    dr = 1.+0.033*cos(2.*pi/365.*J)
    
    Ket = 24./pi*4.92*dr*(omega*sin(lat_r)*sin(delta) + cos(lat_r)*cos(delta)*sin(omega))
    Kso = (0.75+2e-5*elev)*Ket
    Kr = K/Kso
    
    Kr=1. if Kr > 1. else Kr
    Kr=0. if Kr < 0. else Kr
    
    fcd = 1. - Kr
    return fcd
    
