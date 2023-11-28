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

## constants
waterds = 1000.0
cw = 0.0042
ca = 1.013
sigma = 4.9e-9
T_abs = 273.15
alb = 0.1
tstep = 1
