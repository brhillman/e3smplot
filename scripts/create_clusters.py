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
import ngl
import numpy as np
from sklearn import cluster

import os.path
import datetime

import matplotlib.pyplot as plt

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler


# divide earth into 40 deg x 40 deg scenes
  # Non-overlapping to begin with
Location_quad=np.zeros((9,4))
Location_quad[:,0]=-20.
Location_quad[:,1]=20.
Every_five=np.arange(-180,181,40)
Every_five_inv=Every_five[:]
Location_quad[:,2]=Every_five_inv[:-1]
Location_quad[:,3]=Every_five_inv[1:]


# ------------------------------------------------
# USER SPECIFIED parameters
# ------------------------------------------------

# Specify latitude band to examine (relevant for global data)
lat1=-20
lat2=20

# Specify precipitation rate cutoff (tunable parameter)
PRECT_coefficient=3600. * 1000. * 24.
PRECT_cutoff=1.2

#Specify eps parameter which sets the distance used to do clustering
eps_width=0.01 #0.02 maybe

# Specify coefficient to multiply with area to get actual area in km**2
domain_area_coeff=1./(4.*np.pi)*5.101e8 # area integrates to 4*pi so divide by that and  multiply by total area of Earth

dir_loc='/global/cfs/cdirs/e3sm/terai/SCREAM/Analysis/DYAMOND2/RainClusters/' #For outputting Cluster info (currently outputs into text files but will be changed to DataFrames and csv)

time_slices=48
time_stride=2

days=['02-20']

# --------------------------------------------------------
# --------------------------------------------------------

# Cycles through days and create clusters

for iii in days:
    print (iii)
    
    # input_file_name is the history file with the instantaneous PRECT output
    input_file_name=''.join(['/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127/run/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h1.2020-',iii,'-00000.nc'])
    ds_h1 = xarray.open_dataset( input_file_name )
    ncol=ds_h1['ncol']
    
    # Calculate the size distribution from the data


    # create array with size, location, and mean conditions
    # [time steps, 9 scenes * anticipated max # of clusters in each scene]
    Feb16_size=np.zeros((48,9*6000))
    Feb16_lat=np.zeros((48,9*6000))
    Feb16_lon=np.zeros((48,9*6000))
    Feb16_PRECT=np.zeros((48,9*6000))
    
    Total_rain_area=np.zeros((48,9))
    Total_area=np.zeros((48,9))


    Feb16_size[:,:]=np.nan

    for i in np.arange(time_slices):   #loop through time steps
        time_slice=int(i*time_stride)  
        ii=0
        for j in np.arange(9): #loop through scenes
            #print(j)
            lat1,lat2,lon1,lon2=Location_quad[j,:]
            mask = xarray.DataArray( np.ones([len(ncol)],dtype=bool), coords=[('ncol', ncol)], dims='ncol' )
            if 'lat1' in vars(): mask = mask & (ds_h1['lat']>=lat1)
            if 'lat2' in vars(): mask = mask & (ds_h1['lat']<=lat2)
            if lon1 < 0:
                if 'lon1' in vars(): mask = mask & ((ds_h1['lon']>=lon1+360.) | (ds_h1['lon']<=lon2))
            else:
                if 'lon1' in vars(): mask = mask & (ds_h1['lon']>=lon1)
                if 'lon2' in vars(): mask = mask & (ds_h1['lon']<=lon2)

            #Subset lat, lon, PRECT
            lat_area=ds_h1['lat'].where(mask,drop=True)
            lon_area=ds_h1['lon'].where(mask,drop=True)
            
            PRECT_area=ds_h1['PRECT'].isel(time=time_slice).where(mask,drop=True)
            PRECT_area=PRECT_area * PRECT_coefficient #from m/s to mm/d 
            lat_array=np.array(lat_area)
            lon_array=np.array(lon_area)


            #Create rain mask
            PRECT_one_snapshot=np.squeeze(np.array(PRECT_area))
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

            
            lon_1D=lon_array*np.cos(lat_array/180.*np.pi)  # account for change in lon with latitude
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
            """
            fig=plt.figure(figsize=(12,12))
            ax=fig.add_subplot(1,1,1)
            cm = plt.cm.get_cmap('jet')
            sc=plt.scatter(lon_1D,lat_1D,
                           c=PRECT_1D,vmin=1, vmax=50, s=4,edgecolors='none',cmap=cm)
            colbar=plt.colorbar(sc)
            colbar.set_label(label='mm/d',fontsize=14)
            colbar.ax.tick_params(labelsize=14)
            plt.yticks(size=14)
            plt.xticks(size=14)
            plt.xlabel('Longitude',fontsize=16)
            plt.ylabel('Latitude',fontsize=16)
            plt.title('PRECT (Feb 20 0:00)',fontsize=16)
            plt.savefig('RainClusters_PRECT_Feb16.png')


            label_colors=np.arange(-1,5000)
            label_colors_20=label_colors-np.floor(label_colors/20.)*20.
            # Black removed and is used for noise instead.
            fig=plt.figure(figsize=(12,12))
            ax=fig.add_subplot(1,1,1)
            unique_labels = set(labels)
            cm = plt.cm.get_cmap('tab20')

            for k in unique_labels:
                class_member_mask = (labels == k)

                if k == -1:
                    # Black used for noise.
                    xy = X[class_member_mask & core_samples_mask]
                    #plt.plot(xy[:, 0], xy[:, 1],'o',markersize=2,markerfacecolor='k',markeredgecolor='k')
                    xy = X[class_member_mask & ~core_samples_mask]
                    #plt.plot(xy[:, 0], xy[:, 1],'o',markersize=1,markerfacecolor='k',markeredgecolor='k')
                else:
                    xy = X[class_member_mask & core_samples_mask]
                    plt.scatter(xy[:, 0], xy[:, 1],c=label_colors_20[k+1]*np.ones((xy.shape[0],1)),s=4,vmin=0, vmax=19,edgecolors='none',cmap=cm)#, 
                    xy = X[class_member_mask & ~core_samples_mask]
                    plt.scatter(xy[:, 0], xy[:, 1],c=label_colors_20[k+1]*np.ones((xy.shape[0],1)),s=2,vmin=0, vmax=19,edgecolors='none',cmap=cm)


            fig=plt.figure(figsize=(12,12))
            ax=fig.add_subplot(1,1,1)
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
                    xy = X_rains[class_member_mask & core_samples_mask]
                    plt.scatter(xy[:, 0], xy[:, 1],c=label_colors_20[k+1]*np.ones((xy.shape[0],1)),s=4,vmin=0, vmax=19,edgecolors='none',cmap=cm)#, 
                    xy = X_rains[class_member_mask & ~core_samples_mask]
                    plt.scatter(xy[:, 0], xy[:, 1],c=label_colors_20[k+1]*np.ones((xy.shape[0],1)),s=2,vmin=0, vmax=19,edgecolors='none',cmap=cm)

            plt.title('Feb 20: Estimated number of clusters: %d' % n_clusters_,fontsize=16)
            plt.yticks(size=14)
            plt.xticks(size=14)
            plt.xlabel('Longitude',fontsize=16)
            plt.ylabel('Latitude',fontsize=16)
            plt.savefig('RainClusters_Clusters_Feb16.png')
            """

            
            # Longitude and Latitude diagnostics
            diag_lon_1D=lon_array.copy()
            diag_lat_1D=lat_array.copy()

            diag_lon_rain_only=diag_lon_1D[RAIN_mask>0.3]
            diag_lat_rain_only=diag_lat_1D[RAIN_mask>0.3]


            # Precip rate diagnostic
            PRECT_clouds_only=PRECT_1D[RAIN_mask>0.3]


            #Retrieve area data from  the dataset
            try:
                area_area=ds_h1['area'].where(mask,drop=True)
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
            Total_area[i,j]=np.sum(area_1D)/(4*np.pi)*5.101e8 #in sq km
            Total_rain_area[i,j]=np.sum(area_clouds_only)*domain_area_coeff #in sq km

    np.savetxt(''.join([dir_loc,'SCREAM_native_RainCluster_size_Feb',iii,'.txt']),Feb16_size,delimiter=',')
    np.savetxt(''.join([dir_loc,'SCREAM_native_RainCluster_lat_Feb',iii,'.txt']),Feb16_lat,delimiter=',')
    np.savetxt(''.join([dir_loc,'SCREAM_native_RainCluster_lon_Feb',iii,'.txt']),Feb16_lon,delimiter=',')
    np.savetxt(''.join([dir_loc,'SCREAM_native_RainCluster_PRECT_Feb',iii,'.txt']),Feb16_PRECT,delimiter=',')
    np.savetxt(''.join([dir_loc,'SCREAM_native_RainCluster_RainArea_Feb',iii,'.txt']),Total_rain_area,delimiter=',')
    np.savetxt(''.join([dir_loc,'SCREAM_native_RainCluster_TotalArea_Feb',iii,'.txt']),Total_area,delimiter=',')
    
