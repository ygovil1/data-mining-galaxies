# import basic libraries
import pickle
import yt
import numpy as np
import matplotlib.pyplot as plt
from math import log, log10, pi
from astropy import units as u

# import libraries - not sure what they do
# used to ensure halo catalog loads properly
import tempfile
import shutil
import os

# Create temporary directory for storing files
tmpdir = tempfile.mkdtemp()

# import halo catalogue func
from yt.analysis_modules.halo_analysis.api import *

# load halo dataset
halos_ds = yt.load('../halo_catalogs/catalog/catalog0070_thres160.0.h5')

# load raw dataset
ds = yt.load('~/../../tigress/cen/LG4_2048_40pc/RD0070/redshift0070')
ds.print_stats()

# Instantiate a catalog using those two paramter files
hc = HaloCatalog(halos_ds=halos_ds, output_dir=os.path.join(tmpdir, 'halo_catalog'))
hc.load()

# load python halo list
with open('../halo_list0070', 'rb') as file1:
    halo_list = pickle.load(file1)

# load mass of smallest halo
with open('../finest_particle0070', 'rb') as file2:
    min_mass = pickle.load(file2)
    min_mass = 100 * min_mass * u.g # 100*finest_particle

# load halo count from file
with open('../count_halo0070_160', 'rb') as file3:
    count = pickle.load(file3)
    
# use threshold as specified by Renyue
threshold = 200 * 5.92e-28 * (u.g / (u.cm ** 3))

# min and max bounds for radial profile
# max = 0.5 Mpc comoving
rad_min, rad_max = 2e21, 3.71e23

# specify boundaries of zoom-in box
# scaling factor multiplied by info from text file 
# units in cm
scaling = 2.22535525e+25 # scales dataset coords to cm
xmin = scaling*0.39319589 * u.cm
ymin = scaling*0.42984636 * u.cm
zmin = scaling*0.41706725 * u.cm
xmax = scaling*0.56298484 * u.cm
ymax = scaling*0.55089246 * u.cm
zmax = scaling*0.56698254 * u.cm

# array to store old and new mass and radius of halos
# orig_mass, orig_v_radius, new_mass, new_v_radius
halo_array = []

for i in range(100,200):
    halo = halo_list[i]
    
    # create array to store info of this halo
    halo_info = []
    
    # find parameters of halo
    index = halo[0]
    x = halo[1]
    y = halo[2]
    z = halo[3]
    center = [x.value/scaling, y.value/scaling, z.value/scaling]
    mass = halo[4]
    radius = halo[5]
    
    # store original mass and radius in info list
    halo_info.append(index)
    halo_info.append(mass.value)
    halo_info.append(radius.to('kpc').value)
    
    # check if halo is inside zoom-in box
    if xmin <= x < xmax and ymin <= y < ymax and zmin <= z < zmax:
        pass
    else:
        # add 0's otherwise
        halo_info.append(0)
        halo_info.append(0)
        
        # skip rest of loop
        continue
    
    # check for min mass based on finest particle mass
    if min_mass > mass:
        # add 0's otherwise
        halo_info.append(0)
        halo_info.append(0)
        
        # skip rest of loops
        continue
    
    # create sphere
    sp = ds.sphere(center, (rad_max, 'cm'))
    # create radial density profile
    rp = yt.create_profile(sp, 'radius', 'density', accumulation=True, 
                           units = {'radius': 'cm'}, 
                           logs = {'radius': True, 'density': True}, 
                           n_bins = 64, 
                           extrema = {'radius': (rad_min, rad_max)})
    
    # find radius and density where density > threshold
    bool_mask = rp['density'] > threshold.value
    thresh_rad = rp.x[bool_mask]
    thresh_dens = rp['density'][bool_mask]
    
    # check if density is ever above threshold
    if thresh_rad.size > 1:
        # initialize appended quantities
        new_mass = 0 * u.cm
        new_rad = 0 * u.cm
        
        # find boundary radius and density at that radius, and index of that bin 
        rad1 = thresh_rad[-1] * u.cm
        dens1 = thresh_dens[-1] * u.g / (u.cm**3)
        index1 = np.where(rp['density']==thresh_dens[-1])[0]
        
        # in the case that all density are above threshold
        if index1[0] + 1 == rp.x.size:
            new_rad = thresh_rad[-1] * u.cm
            new_dens = thresh_dens[-1] * u.g / (u.cm**3)
        else:
            rad2 = rp.x[index1 + 1] * u.cm
            dens2 = rp['density'][index1 + 1] * u.g / (u.cm**3)
            
            # use interpolation to find new radius and new density
            new_rad = ((threshold - dens2)*rad1 + (dens1 - threshold)*rad2) / (dens1 - dens2)
            new_dens = threshold
        
        
        # find new mass = dens * vol
        new_mass = new_dens * (4/3 * pi * (new_rad**3))
        
        
        
        '''
        # FOR DEBUGGING
        mass1 = dens1 * (4/3 * pi * (rad1**3))
        print('new_mass, new_rad, new_dens, mass1', (new_mass, new_rad, new_dens, mass1))
        print('threshold, dens2, rad2', (threshold, dens2, rad2))
        '''
        
        
        
        # add new radius and mass to info list
        halo_info.append(new_mass.to('Msun').value)
        halo_info.append(new_rad.to('kpc').value)
    
    else:
        # add 0's otherwise
        halo_info.append(0)
        halo_info.append(0)
    
    # print result of halo
    print(halo_info)
    
    # append halo info list to array
    halo_array.append(halo_info)

# store list to file
with open('calc_list0070_100_200', 'wb') as outfile:
    pickle.dump(halo_array, outfile)