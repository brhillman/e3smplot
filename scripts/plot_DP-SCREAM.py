#!/usr/bin/env python3

### About the script
# plots any variables to look at in DPSCREAM

from glob import glob
import sys, os
import xarray
import numpy as np
import os.path

import matplotlib.pyplot as plt

# Load model output file using cdms2 

# ********* Change the orig_data_path to fit your case ***************
orig_data_dir='/global/cscratch1/sd/terai/testing/'
orig_data_path=''.join([orig_data_dir,'GATEIII_20x20.p400.eam.h1.1974-08-30-00000.nc'])

# Figures location
fig_dir='/global/cscratch1/sd/terai/testing/Figures/'
run_name='FirstTest'  #for naming figures

# Choose timestep to look at
time_step=240  # Choose timestep

# Choose variable and any scaling coefficients
var='PRECL'
Coeff=3600.*1000.*24.  # =1. if no scaling necessary
units='mm/d'

# **************************************************************************

ds = xarray.open_dataset( orig_data_path )
x_area=np.array(ds['crm_grid_x'])
y_area=np.array(ds['crm_grid_y'])
var_area=np.array(ds[var]) # Choose variable 

fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(1,1,1)
cm = plt.cm.get_cmap('jet')
sc=plt.scatter(x_area[time_step,:]/1000.,y_area[time_step,:]/1000.,
               c=var_area[time_step,:]*Coeff,vmin=0, vmax=50, s=50,edgecolors='none',cmap=cm)
colbar=plt.colorbar(sc)
colbar.set_label(label=units,fontsize=14)
colbar.ax.tick_params(labelsize=14)
plt.yticks(size=14)
plt.xticks(size=14)
plt.xlabel('y (km)',fontsize=16)
plt.ylabel('x (km)',fontsize=16)
plt.title(''.join([var,': time ',str(time_step)]),fontsize=16)
plt.savefig(''.join([fig_dir,'DPSCREAM_',run_name,'_',var,'_ts',"{:03d}".format(time_step),'.png']))
