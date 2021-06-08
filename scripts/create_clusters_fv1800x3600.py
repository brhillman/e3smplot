#!/usr/bin/env python3

"""
# Creating rain clusters from precipitation data
Takes ne1024pg2 precipitation (can be modified for other cloud-related variable) data regridded to fv1800x3600 and then clusters are created based on a precipitation cut-off and DBSCAN (nearest neighbor clustering). Then we retrieve the cluster information, including
* size
* location (lat/lon)
* time
* precipitation rate
* Any others you would want to add
"""


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

# divide earth into longitude bands of 40 deg
  # Non-overlapping to begin with
Location_quad=np.zeros((9,4))
Location_quad[:,0]=-20.
Location_quad[:,1]=20.
Every_five=np.arange(-180,181,40)
Every_five_inv=Every_five[:]
Location_quad[:,2]=Every_five_inv[:-1]
Location_quad[:,3]=Every_five_inv[1:]

input_file_name='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127/run/regridded/LANDFRAC_1800x3600_USGS-gtopo30_ne1024np4pg2_16xconsistentSGH_20190528_converted.nc'
ds_landfrac = xarray.open_dataset( input_file_name )

days=['20','21']
outputdir_loc='/global/cfs/cdirs/e3sm/terai/SCREAM/Analysis/DYAMOND2/RainClusters/'   #User modify
run_identifier='SCREAMv01SPA'                                                         #User modify

for iii in days:
    print (iii)
    input_file_name=''.join(['/global/cfs/cdirs/e3sm/terai/SCREAM/SCREAMv01/SCREAMv010.SCREAM-DY2.ne1024pg2.may21.n1536a8x16.topo1.ts4.spa.cldfrc/rgr_1800x3600/SCREAMv010.SCREAM-DY2.ne1024pg2.may21.n1536a8x16.topo1.ts4.spa.cldfrc.eam.h1.2020-01-',iii,'-00000.nc'])
    ds_h1 = xarray.open_dataset( input_file_name )
    
    input_file_name=''.join(['/global/cfs/cdirs/e3sm/terai/SCREAM/SCREAMv01/SCREAMv010.SCREAM-DY2.ne1024pg2.may21.n1536a8x16.topo1.ts4.spa.cldfrc/rgr_1800x3600/SCREAMv010.SCREAM-DY2.ne1024pg2.may21.n1536a8x16.topo1.ts4.spa.cldfrc.eam.h0.2020-01-',iii,'-00000.nc'])
    ds_h0 = xarray.open_dataset( input_file_name )
    
    input_file_name=''.join(['/global/cfs/cdirs/e3sm/terai/SCREAM/SCREAMv01/SCREAMv010.SCREAM-DY2.ne1024pg2.may21.n1536a8x16.topo1.ts4.spa.cldfrc/rgr_1800x3600/SCREAMv010.SCREAM-DY2.ne1024pg2.may21.n1536a8x16.topo1.ts4.spa.cldfrc.eam.h2.2020-01-',iii,'-00000.nc'])
    ds_h2 = xarray.open_dataset( input_file_name )
    # Calculate the size distribution from the data
    # This gets at the size distribution at noon around the world - so they might differ

    # Specify latitude band to examine
    lat1=-20
    lat2=20

    # Specify precipitation rate cutoff
    PRECT_coefficient=3600. * 1000. * 24.
    PRECT_cutoff=1.2


    Jan16_size=np.zeros((48,9*4000))
    Jan16_lat=np.zeros((48,9*4000))
    Jan16_lon=np.zeros((48,9*4000))
    Jan16_PRECT=np.zeros((48,9*4000))
    Jan16_LANDFRAC=np.zeros((48,9*4000))
    Jan16_FLNT=np.zeros((48,9*4000))
    Jan16_TMCLDICE=np.zeros((48,9*4000))
    Jan16_TMCLDLIQ=np.zeros((48,9*4000))
    Jan16_TMRAINQM=np.zeros((48,9*4000))
    Jan16_TMQ=np.zeros((48,9*4000))
    
    Total_rain_area=np.zeros((48,9))
    Total_area=np.zeros((48,9))


    Jan16_size[:,:]=np.nan

    lat_all=ds_h1['lat']
    lon_all=ds_h1['lon']

    for i in np.arange(48):
        #print('   ')
        #print(i)
        #print('timestep')
        time_slice=int(i*2)  
        ii=0
        for j in np.arange(9):
            #print(j)
            lat1,lat2,lon1,lon2=Location_quad[j,:]
            PRECT_area=ds_h1['PRECT'].isel(time=time_slice).sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            PRECT_area=PRECT_area * PRECT_coefficient #from m/s to mm/d 
            lat_area=np.array(PRECT_area.lat)
            lon_area=np.array(PRECT_area.lon)
            lon_areaB=lon_area-lon_area[0]

            try:
                area_area=ds_h1['area'].sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            except:
                area_area=np.cos(lat_area*pi/180.)*np.ones((len(lon_areaB)))
            PRECT_one_snapshot=np.squeeze(np.array(PRECT_area))
            RAIN_mask=np.zeros(PRECT_one_snapshot.shape)
            RAIN_mask[PRECT_one_snapshot>PRECT_cutoff]=1.

            if np.sum(RAIN_mask)<4:
                print('skipping')
                continue

            # Figure out lon of values
            lon_array,lat_array=np.meshgrid(lon_areaB,lat_area)
            lon_array=lon_array*np.cos(lat_array/180.*np.pi)
            
            lon_areaB_mean=np.mean(lon_areaB)
            lon_areaB_std=np.std(lon_areaB)
            
            lat_areaB_mean=np.mean(lat_area)
            lat_areaB_std=np.std(lat_area)

            lon_1D=np.reshape(lon_array,(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            lat_1D=np.reshape(lat_array,(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            RAIN_mask_1D=np.reshape(RAIN_mask,(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            lon_rain_only=lon_1D[RAIN_mask_1D>0.3]
            lat_rain_only=lat_1D[RAIN_mask_1D>0.3]
            X_rains=np.zeros((len(lon_rain_only),2))
            X_rains[:,0]=lon_rain_only
            X_rains[:,1]=lat_rain_only

            # Perform clustering using DBSCAN
            #X = StandardScaler().fit_transform(X_rains)   #instead of using standardization, use the following algorithm to standardize myself
            X = X_rains.copy()
            X[:,0]=(X_rains[:,0]-lon_areaB_mean)/lon_areaB_std
            X[:,1]=(X_rains[:,1]-lat_areaB_mean)/lat_areaB_std
            db = DBSCAN(eps=0.01, min_samples=2).fit(X)     #eps is a tuning parameter that is sensitive to the spacing
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_

            # Number of clusters in labels, ignoring noise if present.
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            n_noise_ = list(labels).count(-1)

            #print('Estimated number of clusters: %d' % n_clusters_)

            PRECT_1D=np.reshape(PRECT_one_snapshot,(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
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
            plt.title('PRECT (Jan 20 0:00)',fontsize=16)
            plt.savefig('RainClusters_PRECT_Jan20.png')


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

            plt.title('Jan 20: Estimated number of clusters: %d' % n_clusters_,fontsize=16)
            plt.yticks(size=14)
            plt.xticks(size=14)
            plt.xlabel('Longitude',fontsize=16)
            plt.ylabel('Latitude',fontsize=16)
            plt.savefig('RainClusters_Clusters_Jan20.png')
            """
            # For diagnostics
            diag_lon_array,diag_lat_array=np.meshgrid(lon_area,lat_area)
            diag_lon_1D=np.reshape(diag_lon_array,(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            diag_lat_1D=np.reshape(diag_lat_array,(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))

            diag_lon_rain_only=diag_lon_1D[RAIN_mask_1D>0.3]
            diag_lat_rain_only=diag_lat_1D[RAIN_mask_1D>0.3]


            PRECT_1D=np.reshape(np.array(PRECT_one_snapshot),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            PRECT_clouds_only=PRECT_1D[RAIN_mask_1D>0.3]

            area_1D=np.reshape(np.array(area_area),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            # Figure out the masked areas
            area_clouds_only=area_1D[RAIN_mask_1D>0.3]

            #obtain land fraction 
            landfrac_area=ds_landfrac['LANDFRAC'].sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            landfrac_1D=np.reshape(np.array(landfrac_area),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            landfrac_clouds_only=landfrac_1D[RAIN_mask_1D>0.3]
            
            #obtain other variables (TMQ, FLNT, TMCLDLIQ, TMRAINQM, TMCLDICE)
            FLNT_area=ds_h2['FLNT'].isel(time=time_slice).sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            FLNT_1D=np.reshape(np.array(FLNT_area),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            FLNT_clouds_only=FLNT_1D[RAIN_mask_1D>0.3]
            
            TMQ_area=ds_h0['TMQ'].isel(time=time_slice).sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            TMQ_1D=np.reshape(np.array(TMQ_area),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            TMQ_clouds_only=TMQ_1D[RAIN_mask_1D>0.3]
            
            TMCLDLIQ_area=ds_h0['TMCLDLIQ'].isel(time=time_slice).sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            TMCLDLIQ_1D=np.reshape(np.array(TMCLDLIQ_area),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            TMCLDLIQ_clouds_only=TMCLDLIQ_1D[RAIN_mask_1D>0.3]
            
            TMRAINQM_area=ds_h0['TMRAINQM'].isel(time=time_slice).sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            TMRAINQM_1D=np.reshape(np.array(TMRAINQM_area),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            TMRAINQM_clouds_only=TMRAINQM_1D[RAIN_mask_1D>0.3]
            
            TMCLDICE_area=ds_h0['TMCLDICE'].isel(time=time_slice).sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            TMCLDICE_1D=np.reshape(np.array(TMCLDICE_area),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            TMCLDICE_clouds_only=TMCLDICE_1D[RAIN_mask_1D>0.3]
            
            TMRAINQM_area=ds_h0['TMRAINQM'].isel(time=time_slice).sel(lat=slice(lat1,lat2),lon=slice(lon1,lon2))
            TMRAINQM_1D=np.reshape(np.array(TMRAINQM_area),(PRECT_one_snapshot.shape[0]*PRECT_one_snapshot.shape[1],1))
            TMRAINQM_clouds_only=TMRAINQM_1D[RAIN_mask_1D>0.3]
            


            unique_labels = set(labels)
            
            for k in unique_labels:
                class_member_mask = (labels == k)
                if k == -1:
                    Jan16_size[i,ii]=0.
                    Jan16_lat[i,ii]=0.
                    Jan16_lon[i,ii]=0.
                    Jan16_PRECT[i,ii]=0.
                    Jan16_LANDFRAC[i,ii]=0.

                else:
                    CLDarea_label=np.array(area_clouds_only[class_member_mask])
                    PRECT_label=np.array(PRECT_clouds_only[class_member_mask])
                    lat_label=np.array(diag_lat_rain_only[class_member_mask])
                    lon_label=np.array(diag_lon_rain_only[class_member_mask])
                    landfrac_label=np.array(landfrac_clouds_only[class_member_mask])
                    
                    FLNT_label=np.array(FLNT_clouds_only[class_member_mask])
                    TMQ_label=np.array(TMQ_clouds_only[class_member_mask])
                    TMCLDLIQ_label=np.array(TMCLDLIQ_clouds_only[class_member_mask])
                    TMCLDICE_label=np.array(TMCLDICE_clouds_only[class_member_mask])
                    TMRAINQM_label=np.array(TMRAINQM_clouds_only[class_member_mask])

                    # Save values for each cloud object
                    Jan16_size[i,ii]=np.sum(CLDarea_label)/(4*np.pi)*5.101e8 #in sq km
                    Jan16_lat[i,ii]=np.sum(CLDarea_label*lat_label)/np.sum(CLDarea_label)
                    Jan16_lon[i,ii]=np.sum(CLDarea_label*lon_label)/np.sum(CLDarea_label)
                    Jan16_PRECT[i,ii]=np.sum(CLDarea_label*PRECT_label)/np.sum(CLDarea_label)
                    Jan16_LANDFRAC[i,ii]=np.sum(CLDarea_label*landfrac_label)/np.sum(CLDarea_label)
                    Jan16_FLNT[i,ii]=np.sum(CLDarea_label*FLNT_label)/np.sum(CLDarea_label)
                    Jan16_TMQ[i,ii]=np.sum(CLDarea_label*TMQ_label)/np.sum(CLDarea_label)
                    Jan16_TMCLDLIQ[i,ii]=np.sum(CLDarea_label*TMCLDLIQ_label)/np.sum(CLDarea_label)
                    Jan16_TMCLDICE[i,ii]=np.sum(CLDarea_label*TMCLDICE_label)/np.sum(CLDarea_label)
                    Jan16_TMRAINQM[i,ii]=np.sum(CLDarea_label*TMRAINQM_label)/np.sum(CLDarea_label)
                ii=ii+1
            Total_area[i,j]=np.sum(area_1D)/(4*np.pi)*5.101e8 #in sq km
            Total_rain_area[i,j]=np.sum(area_clouds_only)/(4*np.pi)*5.101e8 #in sq km

    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_size_Jan',iii,'.txt']),Jan16_size,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_lat_Jan',iii,'.txt']),Jan16_lat,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_lon_Jan',iii,'.txt']),Jan16_lon,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_PRECT_Jan',iii,'.txt']),Jan16_PRECT,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_LANDFRAC_Jan',iii,'.txt']),Jan16_LANDFRAC,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_FLNT_Jan',iii,'.txt']),Jan16_FLNT,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_TMQ_Jan',iii,'.txt']),Jan16_TMQ,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_TMCLDLIQ_Jan',iii,'.txt']),Jan16_TMCLDLIQ,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_TMCLDICE_Jan',iii,'.txt']),Jan16_TMCLDICE,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_TMRAINQM_Jan',iii,'.txt']),Jan16_TMRAINQM,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_RainArea_Jan',iii,'.txt']),Total_rain_area,delimiter=',')
    np.savetxt(''.join([outputdir_loc,run_identifier,'_RainCluster_TotalArea_Jan',iii,'.txt']),Total_area,delimiter=',')
    
