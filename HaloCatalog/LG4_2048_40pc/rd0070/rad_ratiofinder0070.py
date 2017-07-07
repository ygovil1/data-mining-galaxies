# import basic libraries
import pickle
import yt
import numpy as np
import matplotlib.pyplot as plt
from math import log, log10
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

# load radius threshold
# 200 times (the mean density of the universe times
# (1-Omega_Baryon/Omega_matter))
with open('rad_threshold0070_2', 'rb') as infile:
    threshold = pickle.load(infile)

# find ratio of computed and catalog radii
radii_ratio = []

for halo in hc.halo_list:    
    # find coord of halo
    x = halo.quantities.get('particle_position_x') * u.cm
    y = halo.quantities.get('particle_position_y') * u.cm
    z = halo.quantities.get('particle_position_z') * u.cm
    center = [x.value/scaling, y.value/scaling, z.value/scaling]
    
    # check if halo is inside zoom-in box
    if xmin <= x < xmax and ymin <= y < ymax and zmin <= z < zmax:
        pass
    else:
        continue
    
    # find values of quantities
    radius = halo.quantities.get('virial_radius') * u.cm
    rad_val = radius.value
    
    # create a sphere around halo
    sp1 = ds.sphere(center, (10*rad_val, 'cm'))
    rad_min, rad_max = 2e21, 10*rad_val
    
    # create radial density profile
    rp = yt.create_profile(sp1, 'radius', 'density', accumulation=True, 
                           units = {'radius': 'cm'}, 
                           logs = {'radius': True, 'density': True}, 
                           n_bins = 64, 
                           extrema = {'radius': (rad_min, rad_max)})
    
    # find max radius where density > threshold
    req_rad = rp.x[rp['density'] > threshold.value]
    if req_rad.size > 0:
        # find ratio and append to list
        calc_rad = req_rad[-1]
        ratio = calc_rad.value / rad_val
        
        radii_ratio.append(ratio)
    else:
        radii_ratio.append(0)

# store list to file
with open('rad_ratiolist0070_2', 'wb') as outfile:
    pickle.dump(radii_ratio, outfile)