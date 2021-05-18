#Import necessary libraries

import matplotlib.pyplot as plt
import numpy as np
import requests

# USER MODIFY 
scene_area=(np.sin(20*np.pi/180.)-np.sin(-20*np.pi/180.))/2*5.101e8 #in sq km

def retrieve_size_characteristics(dir_loc,platform_str,date_str):
    MOD_date_size = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_size_',date_str,'.txt']),delimiter=',')
    MOD_date_PRECT = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_PRECT_',date_str,'.txt']),delimiter=',')
    MOD_date_lat = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_lat_',date_str,'.txt']),delimiter=',')
    MOD_date_lon = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_lon_',date_str,'.txt']),delimiter=',')

    MOD_date_size_reshape=np.reshape(MOD_date_size,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))
    MOD_date_size_cldonly=MOD_date_size_reshape[MOD_date_size_reshape>0]

    MOD_date_PRECT_cldonly=np.reshape(MOD_date_PRECT,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    MOD_date_lat_cldonly=np.reshape(MOD_date_lat,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    MOD_date_lon_cldonly=np.reshape(MOD_date_lon,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    return MOD_date_size_cldonly,MOD_date_PRECT_cldonly,MOD_date_lat_cldonly,MOD_date_lon_cldonly

def retrieve_size_other_characteristics(dir_loc,platform_str,date_str):
    MOD_date_size = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_size_',date_str,'.txt']),delimiter=',')
    MOD_date_PRECT = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_PRECT_',date_str,'.txt']),delimiter=',')
    MOD_date_lat = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_lat_',date_str,'.txt']),delimiter=',')
    MOD_date_lon = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_lon_',date_str,'.txt']),delimiter=',')
    MOD_date_LANDFRAC = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_LANDFRAC_',date_str,'.txt']),delimiter=',')
    
    MOD_date_size_reshape=np.reshape(MOD_date_size,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))
    MOD_date_size_cldonly=MOD_date_size_reshape[MOD_date_size_reshape>0]

    MOD_date_PRECT_cldonly=np.reshape(MOD_date_PRECT,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    MOD_date_lat_cldonly=np.reshape(MOD_date_lat,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    MOD_date_lon_cldonly=np.reshape(MOD_date_lon,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    MOD_date_LANDFRAC_cldonly=np.reshape(MOD_date_LANDFRAC,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    return MOD_date_size_cldonly,MOD_date_PRECT_cldonly,MOD_date_lat_cldonly,MOD_date_lon_cldonly,MOD_date_LANDFRAC_cldonly

def chord_histogram(input_array):
    input_chord_cldonly=(input_array/(np.pi))**0.5
    Diameter_length_binedges=2**np.arange(3,10.1,0.4)
    counts2,bin_edgesb = np.histogram(input_chord_cldonly,bins=Diameter_length_binedges)
    bin_width=Diameter_length_binedges[1:]-Diameter_length_binedges[:-1]
    Box_area=scene_area
    counts_norm=counts2/bin_width/Box_area
    counts_errors=(counts2)**0.5/bin_width/Box_area
    return counts_norm,counts_errors

def size_PRECT_histogram(size_array,PRECT_array,radius_edge=None,prect_edge=None):
    if radius_edge is None:
        radius_edge=2**np.arange(3,10.1,0.4)
    if prect_edge is None:
        prect_edge=2**np.arange(0.2,10.1,0.4)
    input_chord_cldonly=(size_array/(np.pi))**0.5
    counts,xedges,yedges = np.histogram2d(input_chord_cldonly,PRECT_array,bins=(radius_edge,prect_edge))
    x_mesh,y_mesh=np.meshgrid(xedges,yedges)
    counts=counts.T
    return counts,x_mesh,y_mesh

def chord_histogram_areaPRECT_weight(input_array,prect_area_array):
    input_chord_cldonly=(input_array/(np.pi))**0.5
    Diameter_length_binedges=2**np.arange(3,11.5,0.4)
    counts_areaweight,bin_edgesb = np.histogram(input_chord_cldonly,bins=Diameter_length_binedges,weights=input_array)
    PRECT_weight=input_array*prect_area_array
    counts_PRECTweight,bin_edgesb = np.histogram(input_chord_cldonly,bins=Diameter_length_binedges,weights=PRECT_weight)
    return counts_areaweight,counts_PRECTweight,bin_edgesb

def size_PRECT_histogram_areaPRECT_weight(size_array,PRECT_array,radius_edge=None,prect_edge=None):
    if radius_edge is None:
        radius_edge=2**np.arange(3,11.5,0.4)
    if prect_edge is None:
        prect_edge=2**np.arange(0.2,10.1,0.4)
    input_chord_cldonly=(size_array/(np.pi))**0.5
    PRECT_weight=size_array*PRECT_array
    counts_area,xedges,yedges = np.histogram2d(input_chord_cldonly,PRECT_array,bins=(radius_edge,prect_edge),weights=size_array)
    counts_PRECT,xedges,yedges = np.histogram2d(input_chord_cldonly,PRECT_array,bins=(radius_edge,prect_edge),weights=PRECT_weight)
    x_mesh,y_mesh=np.meshgrid(xedges,yedges)
    counts_area=counts_area.T
    counts_PRECT=counts_PRECT.T
    return counts_area,counts_PRECT,x_mesh,y_mesh


# Retrieve the data pertaining to clusters from the text files 
dir_loc='/global/cfs/cdirs/e3sm/terai/SCREAM/Analysis/DYAMOND2/RainClusters/'
#SCREAM_native_Feb20_size_cldonly,SCREAM_native_Feb20_PRECT_cldonly,SCREAM_native_Feb20_lat_cldonly,SCREAM_native_Feb20_lon_cldonly=retrieve_size_characteristics(dir_loc,'SCREAM_native','Feb20')
SCREAM_Feb20_size_cldonly,SCREAM_Feb20_PRECT_cldonly,SCREAM_Feb20_lat_cldonly,SCREAM_Feb20_lon_cldonly,SCREAM_Feb20_LANDFRAC_cldonly=retrieve_size_other_characteristics(dir_loc,'SCREAM','Feb20')
SCREAM_Feb21_size_cldonly,SCREAM_Feb21_PRECT_cldonly,SCREAM_Feb21_lat_cldonly,SCREAM_Feb21_lon_cldonly,SCREAM_Feb21_LANDFRAC_cldonly=retrieve_size_other_characteristics(dir_loc,'SCREAM','Feb21')

GPM_Jan20_size_cldonly,GPM_Jan20_PRECT_cldonly,GPM_Jan20_lat_cldonly,GPM_Jan20_lon_cldonly,GPM_Jan20_LANDFRAC_cldonly=retrieve_size_other_characteristics(dir_loc,'GPM','Jan20')
GPM_Jan21_size_cldonly,GPM_Jan21_PRECT_cldonly,GPM_Jan21_lat_cldonly,GPM_Jan21_lon_cldonly,GPM_Jan21_LANDFRAC_cldonly=retrieve_size_other_characteristics(dir_loc,'GPM','Jan21')


#SCREAM_native_Feb20_counts_normalized,SCREAM_native_Feb20_counts_error=chord_histogram(SCREAM_native_Feb20_size_cldonly)
SCREAM_Feb20_counts_normalized,SCREAM_Feb20_counts_error=chord_histogram(SCREAM_Feb20_size_cldonly)
SCREAM_Feb21_counts_normalized,SCREAM_Feb21_counts_error=chord_histogram(SCREAM_Feb21_size_cldonly)

GPM_Jan20_counts_normalized,GPM_Jan20_counts_error=chord_histogram(GPM_Jan20_size_cldonly)
GPM_Jan21_counts_normalized,GPM_Jan21_counts_error=chord_histogram(GPM_Jan21_size_cldonly)

Diameter_length_binedges=2**np.arange(3,10.1,0.4)

fig=plt.figure(figsize=(12,12))
ax=fig.add_subplot(1,1,1)
#plt.plot((Diameter_length_binedges[1:]*Diameter_length_binedges[:-1])**0.5,SCREAM_native_Feb20_counts_normalized,'*:',color='cyan',label='SCREAMv0native')
plt.plot((Diameter_length_binedges[1:]*Diameter_length_binedges[:-1])**0.5,SCREAM_Feb20_counts_normalized,'o:',color='tab:blue')
plt.plot((Diameter_length_binedges[1:]*Diameter_length_binedges[:-1])**0.5,SCREAM_Feb21_counts_normalized,'o:',color='tab:blue',label='SCREAMv0')
plt.plot((Diameter_length_binedges[1:]*Diameter_length_binedges[:-1])**0.5,GPM_Jan20_counts_normalized,'o:',color='tab:orange',label='GPM')#,label='GPM')
plt.plot((Diameter_length_binedges[1:]*Diameter_length_binedges[:-1])**0.5,GPM_Jan21_counts_normalized,'o:',color='tab:orange')


#plt.plot(lon_cold_location,lat_cold_location,'*',color='tab:gray',markersize=10)
ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_ylim([1e-9, 1e-3])
plt.yticks(size=16)
plt.xticks(size=16)
plt.legend(loc='upper right',fontsize=16)
plt.xlabel('log10(Rain obj size) (km)',fontsize=16)
plt.ylabel('log10(PDF) (km$^{-3}$)',fontsize=16)
plt.savefig('RainClusters_test.png')
