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

from functions import *
import constants
from math import log, exp, sqrt, pow

def equilibrim(lat, depth, area, elev, solrad, longrad, ta, vpd, ut, tw0, fch, mth):

    ##################### Inputs #####################
    # lat - latitude of the lake (degree)
    # depth - depth of the water body (m)
    # area - area of the water body (km2)         (not used, use fetch instead)
    # elev - elevation of the water body (m)
    
    # solrad - incoming surface solar radiation (W m-2 per day)
    # longrad - imcoming long wave radiation (W m-2 per day)
    # ta - air temperature (deg. C)
    # vpd - vapor pressure deficit (kPa)
    # ut - wind speed at 2m (m s-1)
    # tw0 - temperature of the water on the previous time step (deg.C)
    # fch - fetch length of the water body for each month (m)

    ##################### Constants #####################
    # constants.waterds - density of water (kg m-3)
    # constants.cw - specific heat of water (MJ kg-1 deg.V-1)
    # constants.ca - specific heat of air (kJ kg-1 deg.C-1)
    # constants.sigma - stefan-boltzmann constant (MJ m-2 deg.C-4 d-1)
    # constants.T_abs - difference between degrees kelvin and degrees celsius
    # constants.alb - albedo of the water body
    # constants.tstep - the time step for the model to use (days)

    ##################### Variables #####################
    # alambda - latent heat of vaporisation (MJ kg-1)
    # gamma - pschrometric constant (kPa deg.C-1)
    # airds - density of air (kg m-3)

    # deltaa - slope of the temperature-satuartion water vapour curve at air temperature (kPa deg C-1)
    # deltawb - slope of the temperature-satuartion water vapour curve at wet bulb temperature (kPa deg C-1)
    # deltaw - slope of the temperature-satuartion water vapour curve at water temperature (kPa deg C-1)
    
    # atmp - atmospheric pressure at specific elevation (kPa)
    # windf - wind function 
    
    # rn - net radiation (MJ m-2 per day)
    # heat_stg - change in heat storage (MJ m-2 per day)
    # tau - time constant of the water body (days)
    # te - equilibrium temperature (deg. C)
    # tw - temperature of the water at the end of the time period (deg.C)
    # le - latent heat flux (w m-2 per day)
    # evap - evaporation calculated using the penman-monteith formula w/o heat storage(mm per day)
    
    # ierr - error flag
        # 0 = ok
        # 1 = albedo  = < 0 or  = > 1
        # 2 = depth  = < 0
        # 3 = air temperature < wet bulb temperature
        # 4 = downwelling solar radiation  = < 0
        # 5 = wind speed < 0.01 m/s
        # 6 = vpd  = < 0
        # 7 = vpd > vsat

    # initialize constants 
    waterds = constants.waterds
    cw = constants.cw
    ca = constants.ca
    sigma = constants.sigma
    T_abs = constants.T_abs
    alb = constants.alb
    tstep = constants.tstep
    
    ierr = 0
    # check input data errors
    ierr = 1 if alb <= 0.0 or alb >= 1.0 else ierr
    ierr = 2 if depth <= 0 else ierr
    depth = 20.0 if depth>20.0 else depth
    ierr = 3 if solrad <= 0. else ierr
    [ierr, ut] = [4, 0.01] if ut <= 0.01 else [ierr, ut]
    [ierr, vpd] = [5, 0.0001] if vpd <= 0.0 else [ierr, vpd]

##############################################################################################
########################### Calculate water equilibrium temperature ##########################
##############################################################################################

## some variables
    alambda = alambdat(ta)
    gamma = psyconst(atm_p(ta, elev), alambda)  ##S2.10 kPa
    airds = airdens(ta, elev)

## slope of the saturation water vapour curve at the air temperature (kPa deg C-1) 
    deltaa = delcalc(ta)

    
## actual vapor pressure and saturated vapor pressure (kPa) 
    es = v_sat(ta)
    ea = es - vpd  
    ea = 0.01 if ea < 0 else ea
    [ierr, vpd] = [6, es*0.99] if es < vpd else [ierr, vpd]

    
## slope of the saturation water vapour curve at the wet bulb temperature (kPa deg C-1) */
    twb = tvpd2wbt(ta, vpd)

    
    [ierr, twb] = [7, ta] if twb > ta else [ierr, twb]
    deltawb = delcalc(twb)
    

## Emissvity of air and water (unitless)
    sradj = solrad*0.0864       ## convert from W m-2 to MJ m-2 d-1

    fcd = cloud_factor(sradj, mth, lat, elev)
    em_a = 1.08*(1.-exp(-pow(ea*10.,(ta+T_abs)/2016.)))*(1+0.22*pow(fcd,2.75))
    em_w = 0.97
    longrad = -9999
    lradj = em_a*sigma*pow((ta+T_abs),4.) if longrad == -9999 else longrad*0.0864

## wind function using the method of McJannet, 2012 (MJ m-2 d-1 kPa-1) */    
    windf = (2.33+1.65*ut)*pow(fch, -0.1)*alambda

## calculate regression coefficients 
## Stream temperature/air temperature relationship: a physical interpretation (O. Mohseni , H.G. Stefan) */
    dcu = windf*(deltaa+gamma)                                  ## (MJ/m2/d)
    B0 = (0.46*em_a+dcu)/(0.46*em_w+dcu)                        ## unitless 
    B1 = ((1.-alb)*sradj-28.38*(em_w-em_a))/(0.46*em_w+dcu)     ## deg. C 
    B2 = windf/(em_w*0.46+dcu)

    te = B0*ta+B1-B2*(es-ea)

###############################################################################################
############################## Calculate water column temperature #############################
###############################################################################################

## time constant (d) 
    tau = (waterds*cw*depth)/(4.0*sigma*pow((twb+T_abs),3.)+windf*(deltawb+gamma))
## water column temperature (deg. C) 
    tw = te+(tw0-te)*exp(-tstep/tau)
    tw = 0. if tw<0. else tw
## change in heat storage (MJ m-2 d-1)
    heat_stg = waterds*cw*depth*(tw-tw0)/tstep

################################################################################################
################################### Calculate the evaporation ##################################
################################################################################################
    #print (sradj,alb,lradj,em_w,ta,T_abs)
## calculate the Penman evaporation
    rn = sradj*(1.-alb)+lradj-em_w*(sigma*pow((ta+T_abs),4.))

    le = (deltaa*(rn-heat_stg)+gamma*windf*vpd)/(deltaa+gamma)
    evap_hs = le/alambda
    
    #evap_hs = 0 if evap_hs < 0 else evap_hs

    le = (deltaa*(rn)+gamma*windf*vpd)/(deltaa+gamma)
    evap_nohs = le/alambda
    #evap_nohs = 0 if evap_nohs < 0 else evap_nohs
    return tw, evap_hs, evap_nohs,ierr

