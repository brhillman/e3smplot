#!/usr/bin/env python3

### About the script
"""
# Creating rain clusters from precipitation data
Takes ne1024pg2 precipitation (can be modified for other cloud-related variable) data and then clusters are created based on a precipitation cut-off and DBSCAN (nearest neighbor clustering). Then we retrieve the cluster information, including
* size
* location (lat/lon)
* time
* precipitation rate
* Any others you would want to add
"""

# Load relevant libraries

from glob import glob
import sys, os
import xarray
import numpy as np
from sklearn import cluster

import os.path
import datetime

import matplotlib.pyplot as plt

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler


# ********* Change the orig_data_path to fit your case ***************
orig_data_dir='/global/cscratch1/sd/terai/testing/'
orig_data_path=''.join([orig_data_dir,'GATEIII_20x20.p400.eam.h1.1974-08-30-00000.nc'])

# Directory for output of clustering
dir_loc='/global/cscratch1/sd/terai/testing/Clusters/'
# Figures location
fig_dir='/global/cscratch1/sd/terai/testing/Figures/'
run_name='FirstTest'  #for naming figures

time_stride=1      # assumes you want every snapshot of data
time_slices=450    # of timeslices to count up to
eps_width=0.1
PRECT_coefficient=3600.*24.*1000.
PRECT_cutoff=10 #10 mm/d cutoff
# ********* Change the orig_data_path to fit your case ***************


ds = xarray.open_dataset( orig_data_path )

area=ds['area']
scene_area=np.sum(np.array(area))/1000000.
domain_area_coeff=1/1000000.

ncol=ds['ncol']

# create array with size, location, and mean conditions
# [time steps, 9 scenes * anticipated max # of clusters in each scene]
Feb16_size=np.zeros((time_slices,2000))
Feb16_lat=np.zeros((time_slices,2000))
Feb16_lon=np.zeros((time_slices,2000))
Feb16_PRECT=np.zeros((time_slices,2000))

Total_rain_area=np.zeros((time_slices,9))
Total_area=np.zeros((time_slices,9))

Feb16_size[:,:]=np.nan

for i in np.arange(time_slices):   #loop through time steps
    time_slice=int(i*time_stride)  
    ii=0
    
    #Subset lat, lon, PRECT
    lat_area=ds['crm_grid_y'].isel(time=time_slice)
    lon_area=ds['crm_grid_x'].isel(time=time_slice)

    PRECT_area=ds['PRECL'].isel(time=time_slice)
    PRECT_area=PRECT_area * PRECT_coefficient #from m/s to mm/d 
    lat_array=np.array(lat_area)
    lon_array=np.array(lon_area)


    #Create rain mask
    PRECT_one_snapshot=np.squeeze(np.array(PRECT_area))
    print(PRECT_one_snapshot.shape)
    RAIN_mask=np.zeros(PRECT_one_snapshot.shape)
    RAIN_mask[PRECT_one_snapshot>PRECT_cutoff]=1.

    if np.sum(RAIN_mask)<4:  #skip cases where there are fewer than 4 grid boxes that are raining
        print('skipping')
        continue

    # Figure out the mean and standard deviation of lat and lon to normalize the data
    lon_array_mean=np.mean(lon_array)
    lon_array_std=np.std(lon_array)

    lat_array_mean=np.mean(lat_array)
    lat_array_std=np.std(lat_array)


    lon_1D=lon_array
    lat_1D=lat_array
    lon_rain_only=lon_1D[RAIN_mask>0.3] #Pick out lon locations that pass mask
    lat_rain_only=lat_1D[RAIN_mask>0.3] #Pick out lat locations that pass mask

    # Create array with location info
    X_rains=np.zeros((len(lon_rain_only),2))

    X_rains[:,0]=lon_rain_only
    X_rains[:,1]=lat_rain_only

    # ----------------------------------------
    # Clustering step
    # ----------------------------------------
    # Perform clustering using DBSCAN
    #X = StandardScaler().fit_transform(X_rains)   #instead of using standardization, use the following algorithm to standardize myself
    X = X_rains.copy()
    X[:,0]=(X_rains[:,0]-lon_array_mean)/lon_array_std
    X[:,1]=(X_rains[:,1]-lat_array_mean)/lat_array_std
    db = DBSCAN(eps=eps_width, min_samples=2).fit(X)     #eps is a tuning parameter that is sensitive to the spacing
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print('Estimated number of clusters: %d' % n_clusters_)

    PRECT_1D=PRECT_one_snapshot

    #optional map plotting step to see if the clustering is doing what it's supposed to be doing


    label_colors=np.arange(-1,5000)
    label_colors_20=label_colors-np.floor(label_colors/20.)*20.
    # Black removed and is used for noise instead.
    


    # Longitude and Latitude diagnostics
    diag_lon_1D=lon_array.copy()
    diag_lat_1D=lat_array.copy()

    diag_lon_rain_only=diag_lon_1D[RAIN_mask>0.3]
    diag_lat_rain_only=diag_lat_1D[RAIN_mask>0.3]


    # Precip rate diagnostic
    PRECT_clouds_only=PRECT_1D[RAIN_mask>0.3]


    #Retrieve area data from  the dataset
    try:
        area_area=ds['area']
    except:
        print("need to create array with area pertaining to each grid")

    area_1D=np.array(area_area)
    area_clouds_only=area_1D[RAIN_mask>0.3] #rain mask

    unique_labels = set(labels)

    for k in unique_labels:
        class_member_mask = (labels == k)
        if k == -1:
            Feb16_size[i,ii]=0.
            Feb16_lat[i,ii]=0.
            Feb16_lon[i,ii]=0.
            Feb16_PRECT[i,ii]=0.
        else:
            CLDarea_label=np.array(area_clouds_only[class_member_mask])
            PRECT_label=np.array(PRECT_clouds_only[class_member_mask])
            lat_label=np.array(diag_lat_rain_only[class_member_mask])
            lon_label=np.array(diag_lon_rain_only[class_member_mask])

            # Save values for each cloud object and area weight
            Feb16_size[i,ii]=np.sum(CLDarea_label)*domain_area_coeff
            Feb16_lat[i,ii]=np.sum(CLDarea_label*lat_label)/np.sum(CLDarea_label)
            Feb16_lon[i,ii]=np.sum(CLDarea_label*lon_label)/np.sum(CLDarea_label)
            Feb16_PRECT[i,ii]=np.sum(CLDarea_label*PRECT_label)/np.sum(CLDarea_label)
        ii=ii+1

    PRECT_1D_lowcutoff=PRECT_1D.copy()
    PRECT_1D_lowcutoff[PRECT_1D_lowcutoff<10]=np.nan
    fig_idx_list=[40,100,200,300,400,500]
    if i in fig_idx_list:
        fig=plt.figure(figsize=(18,7))
        ax=fig.add_subplot(1,2,1)
        cm = plt.cm.get_cmap('rainbow')
        sc=plt.scatter(np.array(lon_area)/1000.,np.array(lat_area)/1000.,
                       c=PRECT_1D_lowcutoff,vmin=1.2, vmax=100, s=8,edgecolors='none',cmap=cm)
        colbar=plt.colorbar(sc)
        colbar.set_label(label='mm/d',fontsize=14)
        colbar.ax.tick_params(labelsize=14)
        plt.yticks(size=14)
        plt.xticks(size=14)
        plt.xlim([0,201])
        plt.ylim([0,201])
        plt.xlabel('X (km)',fontsize=16)
        plt.ylabel('Y (km)',fontsize=16)
        plt.title(''.join(['PRECT (Jan 21 step','{:03d}'.format(i),')']),fontsize=16)
        ax=fig.add_subplot(1,2,2)
        unique_labels = set(labels)
        cm = plt.cm.get_cmap('tab20')
        for k in unique_labels:
            class_member_mask = (labels == k)
            if k == -1:
                # Black used for noise.
                xy = X_rains[class_member_mask & core_samples_mask]
                #plt.plot(xy[:, 0], xy[:, 1],'o',markersize=2,markerfacecolor='k',markeredgecolor='k')
                xy = X_rains[class_member_mask & ~core_samples_mask]
                #plt.plot(xy[:, 0], xy[:, 1],'o',markersize=1,markerfacecolor='k',markeredgecolor='k')
            else:
                sc2=plt.scatter(diag_lon_rain_only[class_member_mask]/1000., diag_lat_rain_only[class_member_mask]/1000.,c=label_colors_20[k+1]*np.ones((len(diag_lon_rain_only[class_member_mask]),1)),s=8,vmin=0, vmax=19,edgecolors='none',cmap=cm)#, 
        plt.title('Jan 21: Estimated number of clusters: %d' % n_clusters_,fontsize=16)
        colbar=plt.colorbar(sc2)
        plt.yticks(size=14)
        plt.xticks(size=14)
        plt.xlim([0,201])
        plt.ylim([0,201])
        plt.xlabel('X (km)',fontsize=16)
        plt.ylabel('Y (km)',fontsize=16)
        plt.savefig(''.join([fig_dir,'DP_SCREAM_',run_name,'_Clusters_',"{:02d}".format(time_slice),'.png']))
        plt.close()

    del PRECT_area
    del lat_area
    del lon_area
    del lat_array
    del lon_array
np.savetxt(''.join([dir_loc,'DPSCREAM_RainCluster_size_',run_name,'.txt']),Feb16_size,delimiter=',')
np.savetxt(''.join([dir_loc,'DPSCREAM_RainCluster_y_',run_name,'.txt']),Feb16_lat,delimiter=',')
np.savetxt(''.join([dir_loc,'DPSCREAM_RainCluster_x_',run_name,'.txt']),Feb16_lon,delimiter=',')
np.savetxt(''.join([dir_loc,'DPSCREAM_RainCluster_PRECT_',run_name,'.txt']),Feb16_PRECT,delimiter=',')

print('Completed making clusters!')

# ***********************************
# Making histogram with just one set
# ***********************************



def retrieve_size_characteristics(dir_loc,platform_str,run_name):
    MOD_date_size = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_size_',run_name,'.txt']),delimiter=',')
    MOD_date_PRECT = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_PRECT_',run_name,'.txt']),delimiter=',')
    MOD_date_y = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_y_',run_name,'.txt']),delimiter=',')
    MOD_date_x = np.loadtxt(''.join([dir_loc,platform_str,'_RainCluster_x_',run_name,'.txt']),delimiter=',')

    MOD_date_size_reshape=np.reshape(MOD_date_size,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))
    MOD_date_size_cldonly=MOD_date_size_reshape[MOD_date_size_reshape>0]

    MOD_date_PRECT_cldonly=np.reshape(MOD_date_PRECT,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    MOD_date_y_cldonly=np.reshape(MOD_date_y,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    MOD_date_x_cldonly=np.reshape(MOD_date_x,(MOD_date_size.shape[0]*MOD_date_size.shape[1],1))[MOD_date_size_reshape>0]
    return MOD_date_size_cldonly,MOD_date_PRECT_cldonly,MOD_date_y_cldonly,MOD_date_x_cldonly

def chord_histogram(input_array):
    input_chord_cldonly=(input_array/(np.pi))**0.5
    Radius_length_binedges=2**np.arange(1,10.1,0.4)
    counts2,bin_edgesb = np.histogram(input_chord_cldonly,bins=Radius_length_binedges)
    bin_width=Radius_length_binedges[1:]-Radius_length_binedges[:-1]
    Box_area=scene_area
    counts_norm=counts2/bin_width/Box_area
    counts_errors=(counts2)**0.5/bin_width/Box_area
    return counts_norm,counts_errors

MOD_date_size_cldonly,MOD_date_PRECT_cldonly,MOD_date_y_cldonly,MOD_date_x_cldonly=retrieve_size_characteristics(dir_loc,'DPSCREAM',run_name)

counts_norm,counts_error=chord_histogram(MOD_date_size_cldonly)

Radius_length_binedges=2**np.arange(1,10.1,0.4)

fig=plt.figure(figsize=(12,6))
ax=fig.add_subplot(1,2,1)
plt.plot((Radius_length_binedges[1:]*Radius_length_binedges[:-1])**0.5,counts_norm,'o:',color='tab:blue')
ax.set_xscale('log')
#ax.set_yscale('log')
plt.yticks(size=16)
plt.xticks(size=16)
plt.legend(loc='upper right',fontsize=16)
plt.xlabel('log10(Rain obj radius) (km)',fontsize=16)
plt.ylabel('Histogram (km$^{-3}$)',fontsize=16)
ax=fig.add_subplot(1,2,2)
plt.plot((Radius_length_binedges[1:]*Radius_length_binedges[:-1])**0.5,counts_norm/np.sum(counts_norm),'o:',color='tab:blue')
ax.set_xscale('log')
#ax.set_yscale('log')
plt.yticks(size=16)
plt.xticks(size=16)
plt.legend(loc='upper right',fontsize=16)
plt.xlabel('log10(Rain obj radius) (km)',fontsize=16)
plt.ylabel('PDF (km$^{-3}$)',fontsize=16)
plt.savefig(''.join([fig_dir,'RainClusters_histogram_',run_name,'.png']))
