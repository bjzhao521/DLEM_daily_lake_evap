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


from osgeo import ogr
import numpy as np
import math
import shapefile

COEF=111139 ##convert degree to meters

def Nrotate(angle,valuex,valuey,pointx,pointy):
    valuex = np.array(valuex)
    valuey = np.array(valuey)
    nRotatex = (valuex-pointx)*math.cos(angle) - (valuey-pointy)*math.sin(angle) + pointx
    nRotatey = (valuex-pointx)*math.sin(angle) + (valuey-pointy)*math.cos(angle) + pointy
    return nRotatex, nRotatey

def getArea(shapename,ith):
    ##get the ith shapefile's area
    path_to_shp_data=shapename
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(path_to_shp_data, 0)
    layer = dataSource.GetLayer()
    feature=layer[ith]
    geom = feature.GetGeometryRef()
    area = geom.GetArea() 
    return area

def calFetch(shapename,ith):
    sf=shapefile.Reader(shapename)
    Shapes=sf.shapes()
    in_vector=Shapes[ith].points
    area=getArea(shapename,ith)
    fetch_list=[]
    #fetch_list.append(area*COEF*COEF)
    for j in range(6):
        deg=30*j
        out_vector=np.zeros((len(in_vector),2))
        for i in range(len(in_vector)):
            out_vector[i][0],out_vector[i][1]= Nrotate(math.radians(deg),in_vector[i][0],in_vector[i][1],0,0)
        ymax=np.max(out_vector, axis=0)[1]
        ymin=np.min(out_vector, axis=0)[1]
        width=ymax-ymin
        fetch=area/width*COEF
        fetch_list.append(fetch)
    return fetch_list

def Output(shapename):
    sf1=shapefile.Reader(shapename)
    Shapes1=sf1.shapes()
    length=len(Shapes1)
    output_name=shapename.split(".")[0]
    f=open(output_name+'_fetch.csv','w')
    f.write('ID'+','+'0'+','+'30'+','+'60'+','+'90'+','+'120'+','+'150'+','+'180'+','+'210'+','+
            '240'+','+'270'+','+'300'+','+'330'+','+'360'+'\n')  
    for ith in range(length):
        fl=calFetch(shapename,ith)
        ##the angle difference of width and fetch is 90 degree (90=3*30), so started with 3
        f.write(str(ith)+','+str(fl[3])+','+str(fl[4])+','+str(fl[5])+','+str(fl[0])+','+str(fl[1])+','+
                str(fl[2])+','+str(fl[3])+','+str(fl[4])+','+str(fl[5])+','+str(fl[0])+','+str(fl[1])+','+
                str(fl[2])+','+str(fl[3])+'\n')
    return None


def main():
    Output('powell.shp')
    return None
