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
halos_ds = yt.load('./halo_catalogs/catalog/catalog0070_thres160.0.h5')

# load raw dataset
ds = yt.load('~/../../tigress/cen/LG4_2048_40pc/RD0070/redshift0070')

# Instantiate a catalog using those two paramter files
hc = HaloCatalog(halos_ds=halos_ds, output_dir=os.path.join(tmpdir, 'halo_catalog'))
hc.load()

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

# load mass of smallest particle
with open('finest_particle0070', 'rb') as file:
    min_mass = pickle.load(file)
    min_mass = 100 * min_mass # 100*finest_particle

# use threshold as specified by Renyue
threshold = 200 * 5.92e-28 * (u.g / (u.cm ** 3))

# min and max bounds for radial profile
# max = 1.5 Mpc comoving
rad_min, rad_max = 2e21, 1.11e24

# array to store old and new mass and radius of halos
# orig_mass, orig_v_radius, new_mass, new_v_radius
calc_array = []


for halo in hc.halo_list:
    # create list to store info of this halo
    halo_info = []
    
    # find parameters of halo
    x = halo.quantities.get('particle_position_x') * u.cm
    y = halo.quantities.get('particle_position_y') * u.cm
    z = halo.quantities.get('particle_position_z') * u.cm
    center = [x.value/scaling, y.value/scaling, z.value/scaling]
    mass = halo.quantities.get('particle_mass').in_units('Msun')
    radius = halo.quantities.get('virial_radius') * u.cm
    
    # store original mass and radius in info list
    halo_info.append(mass.to('Msun'))
    halo_info.append(radius.to('kpc'))
    
    # check if halo is inside zoom-in box
    if xmin <= x < xmax and ymin <= y < ymax and zmin <= z < zmax:
        pass
    else:
        continue
    
    # check for min mass based on finest particle mass
    if min_mass > mass:
        continue
    
    # create sphere
    sp = ds.sphere(center, (100*radius.value, 'cm'))
    # create radial density profile
    rp = yt.create_profile(sp, 'radius', 'density', accumulation=True, 
                           units = {'radius': 'cm'}, 
                           logs = {'radius': True, 'density': True}, 
                           n_bins = 64, 
                           extrema = {'radius': (rad_min, rad_max)})
    
    # find radius and density where density > threshold
    thresh_rad = rp.x[rp['density'] > threshold.value]
    thresh_dens = rp['density'][rp['density'] > threshold.value]
    
    # check if density is ever above threshold
    if thresh_rad.size > 0:
        # find boundary radius and density at that radius, and index of that bin 
        rad1 = thresh_rad[-1] * u.cm
        dens1 = thresh_dens[-1] * u.g / (u.cm**3)
        index1 = np.where(rp['density']==dens1)[0]
        
        rad2 = rp.x[index1 + 1] * u.cm
        dens2 = rp['density'][index1 + 1] * u.g / (u.cm**3)
        
        # use interpolation to find new radius and new density
        new_rad = ((threshold - dens2)*rad1 + (dens1 - threshold)*rad2)/(dens1 - dens2)
        new_dens = threshold
        
        # find new mass = dens * vol
        new_mass = new_dens * (4/3 * pi * (new_rad**3))
        
        # add new radius and mass to info list
        halo_info.append(new_mass.to('Msun'))
        halo_info.append(new_rad.to('kpc'))
    else:
        # add 0's otherwise
        halo_info.append(0)
        halo_info.append(0)
    
    # append halo info list to array
    calc_array.append(halo_info)

# store list to file
with open('calc_array0070', 'wb') as outfile:
    pickle.dump(calc_array, outfile)