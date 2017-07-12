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
    
# load redshift and Omega values from parameter file
with open('../redshift0070', 'rt') as param_file:
    param_contents = param_file.read()
    
    #redshift
    cindex1 = param_contents.find('CosmologyCurrentRedshift')
    cindex_eq = param_contents.find('=', cindex1)
    cindex2 = param_contents.find('\n', cindex_eq, cindex_eq + 100)
    redshift = float(param_contents[cindex_eq+2:cindex2])
    
    # omega_m
    cindex1 = param_contents.find('CosmologyOmegaMatterNow')
    cindex_eq = param_contents.find('=', cindex1)
    cindex2 = param_contents.find('\n', cindex_eq, cindex_eq + 100)
    omega_m = float(param_contents[cindex_eq+2:cindex2])
    
    # hubble const now
    cindex1 = param_contents.find('CosmologyHubbleConstantNow')
    cindex_eq = param_contents.find('=', cindex1)
    cindex2 = param_contents.find('\n', cindex_eq, cindex_eq + 100)
    hubb_now = float(param_contents[cindex_eq+2:cindex2]) * 100 *(u.km / u.s / u.Mpc)
    hubb_now = hubb_now.to('s**-1')
    
    # cosmological constant
    cindex1 = param_contents.find('CosmologyOmegaLambdaNow')
    cindex_eq = param_contents.find('=', cindex1)
    cindex2 = param_contents.find('\n', cindex_eq, cindex_eq + 100)
    cos_const = float(param_contents[cindex_eq+2:cindex2])
    
    # omega baryon as specified by Renyue
    omega_b = 0.048
    
# calculate hubble const for simulation
hubb_z = hubb_now * ((omega_m * (1 + redshift)**3) + (1 - omega_m))**0.5

# calculate crit density and threshold
GRAV_CONST = (6.67408e-11 * u.m**3 /(u.kg * u.s**2)).to('cm^3*g^-1*s^-2')
crit_dens = (3 * hubb_z**2) / (8 * pi * GRAV_CONST)
omegas = (1 - (omega_b / omega_m))
threshold = 200 * omegas * crit_dens

# min and max bounds for radial profile
# min = 1 kpc proper
# max = 0.5 Mpc comoving
# convert to centimeters value (without astropy units)
rad_min = 1 * u.kpc
rad_max = 0.5 * u.Mpc
rad_max = rad_max / (1 + redshift) # convert to physical

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

for i in range(1013,1100):
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
    sp = ds.sphere(center, (rad_max.to('cm').value, 'cm'))
    # create radial density profile
    rp = yt.create_profile(sp, 'radius', 'Dark_Matter_Density', accumulation=True, 
                           units = {'radius': 'kpc', 'Dark_Matter_Density': 'g/cm**3'}, 
                           logs = {'radius': True, 'Dark_Matter_Density': True}, 
                           n_bins = 64, 
                           weight_field='cell_volume', 
                           extrema = {'radius': (rad_min.to('kpc').value, 
                                                 rad_max.to('kpc').value)})
    
    # --find radius and density where density > threshold\
    # create true an false bool masks
    true_mask = rp['Dark_Matter_Density'] > threshold.value
    # create array with true values of rad and dens
    thresh_rad = rp.x[true_mask]
    thresh_dens = rp['Dark_Matter_Density'][true_mask]
    
    crossings = []
    # find crossings
    for i in range(0, true_mask.size - 1):
        if true_mask[i] == True and true_mask[i+1] == False:
            crossings.append(i)
    
    # check if density is ever above threshold
    if thresh_rad.size > 1:
        # initialize isMethod bool
        isMethod1 = False
        
        # if all dens > threshold sent isMethod to True
        if thresh_rad.size == rp.x.size or len(crossings) == 0:
            isMethod1 = True
        
        
        # pick which method to use
        if isMethod1:
            new_rad = thresh_rad[-1] * u.kpc
            new_dens = thresh_dens[-1] * u.g / (u.cm**3)
        else:
            # find index of boundary bin 
            # set index1 to the first crossing in case there are no others
            index1 = crossings[0]
            for ind in crossings:
                rad1 = thresh_rad[ind] * u.kpc
                if rad1 > radius:
                    index1 = ind
                    break # to ensure that only first rad > rad_catalog is set
            
            # find radius and density before and after threshold
            # for interpolation
            rad1 = rp.x[index1] * u.kpc
            dens1 = rp['Dark_Matter_Density'][index1] * u.g / (u.cm**3)
            rad2 = rp.x[index1 + 1] * u.kpc
            dens2 = rp['Dark_Matter_Density'][index1 + 1] * u.g / (u.cm**3)
            
            # use interpolation to find new radius and new density
            new_rad = ((threshold - dens2)*rad1 + (dens1 - threshold)*rad2) / (dens1 - dens2)
            new_dens = threshold
        
        
        # find new mass = dens * vol
        new_mass = new_dens * (4/3 * pi * (new_rad**3))
        
        # scale by omegas
        new_mass = new_mass / omegas
        print(index, new_mass, new_rad, omegas)
        
        # add new radius and mass to info list depending 
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
with open('calc_list0070_test', 'wb') as outfile:
    pickle.dump(halo_array, outfile)
