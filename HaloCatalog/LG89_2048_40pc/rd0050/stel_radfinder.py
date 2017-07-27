#specify input filenames
halo_catalog_filename = './halo_catalogs/catalog/catalog.0.h5'
dataset_filename = '~/../../tigress/cen/LG89_2048_40pc/RD0050/redshift0050'
redshift_filename = 'redshift0050'

# some constants
start_num_entries = 14 # number of entries halos have at start of this program
end_num_entries = 15 # min number of entries halos have at end of this program
# note- start_num_entires is index of first changed/appended value

# import basic libraries
import pickle
import yt
import numpy as np
import matplotlib.pyplot as plt
from math import log, log10, pi
from astropy import units as u
from operator import itemgetter
# import libraries - not sure what they do
# used to ensure halo catalog loads properly
import tempfile
import shutil
import os
# Create temporary directory for storing files
tmpdir = tempfile.mkdtemp()
# import halo catalogue func and cosmology func
from yt.analysis_modules.halo_analysis.api import *
from yt.utilities.cosmology import Cosmology
# load halo dataset
halos_ds = yt.load(halo_catalog_filename)
# load raw dataset
ds = yt.load(dataset_filename)
# Instantiate a catalog using those two paramter files
hc = HaloCatalog(halos_ds=halos_ds, output_dir=os.path.join(tmpdir, 'halo_catalog'))
hc.load()

# load python halo list
with open('calc_list_3000', 'rb') as file1:
    halo_list = pickle.load(file1)
    count = len(halo_list)

# load redshift and Omega values from parameter file
with open(redshift_filename, 'rt') as param_file:
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
    hubb_now = float(param_contents[cindex_eq+2:cindex2]) *(u.km / u.s / u.Mpc) 
    
    # cosmological constant
    cindex1 = param_contents.find('CosmologyOmegaLambdaNow')
    cindex_eq = param_contents.find('=', cindex1)
    cindex2 = param_contents.find('\n', cindex_eq, cindex_eq + 100)
    cos_const = float(param_contents[cindex_eq+2:cindex2])
    
    # omega baryon as specified by Renyue
    omega_b = 0.048
    
    # box size
    cindex1 = param_contents.find('CosmologyComovingBoxSize')
    cindex_eq = param_contents.find('=', cindex1)
    cindex2 = param_contents.find('\n', cindex_eq, cindex_eq + 100)
    box_size = float(param_contents[cindex_eq+2:cindex2]) * u.Mpc
    
    # left edge strings
    cindex1 = param_contents.find('RefineRegionLeftEdge')
    cindex_eq = param_contents.find('=', cindex1)
    cindex2 = param_contents.find('\n', cindex_eq, cindex_eq + 100)
    left_edge_string = param_contents[cindex_eq+2:cindex2]
    left_edges = left_edge_string.split()
    
    # right edge strings
    cindex1 = param_contents.find('RefineRegionRightEdge')
    cindex_eq = param_contents.find('=', cindex1)
    cindex2 = param_contents.find('\n', cindex_eq, cindex_eq + 100)
    right_edge_string = param_contents[cindex_eq+2:cindex2]
    right_edges = right_edge_string.split()
    
# calculate hubble const for simulation
hubb_z = (100* hubb_now.to('s**-1')) * ((omega_m * (1 + redshift)**3) + (1 - omega_m))**0.5

# calculate crit density and threshold
GRAV_CONST = (6.67408e-11 * u.m**3 /(u.kg * u.s**2)).to('cm^3*g^-1*s^-2')
crit_dens = (3 * hubb_z**2) / (8 * pi * GRAV_CONST)
omegas = (1 - (omega_b / omega_m))
threshold = 200 * omegas * crit_dens

# min and max bounds for radial profile
# min = 0.1 kpc proper
# max = 100 kpc proper
# convert to centimeters value (without astropy units)
rad_min = 0.1 * u.kpc
rad_max = 20 * u.kpc

# specify boundaries of zoom-in box
# scaling factor multiplied by info from text file 
# units in cm
scaling =  ((box_size / hubb_now.value) / (1 + redshift)).to('kpc') # size of box
xmin = scaling*float(left_edges[0])
ymin = scaling*float(left_edges[1])
zmin = scaling*float(left_edges[2])
xmax = scaling*float(right_edges[0])
ymax = scaling*float(right_edges[1])
zmax = scaling*float(right_edges[2])

# store new halo info
new_halo_list = []

# main loop
for halo in halo_list:
    # find parameters of halo
    index = halo[0]
    x = halo[2]
    y = halo[3]
    z = halo[4]
    center = [(x/scaling).value, (y/scaling).value, (z/scaling).value]
    halo_mass = halo[5] # new mass
    radius = halo[7] # new radius
    isSatellite = halo[9]
    hstel_mass = halo[11]
    
    # check if halo is inside zoom-in box
    if xmin <= x < xmax and ymin <= y < ymax and zmin <= z < zmax:
        pass
    else:
        print('OutOfBounds')
        continue

    # check that radius is not 0
    if hstel_mass == 0:
        # print reason
        print('StelMassEqualZero')
        
        # ensure that running this code multiple times doesn't create a long list
        # check if halo already has more info
        if len(halo) > end_num_entries:
            # create copy and replace values
            new_halo = halo
            new_halo[start_num_entries] = 0
        else:
            # append 0's to halo_info if data not applicable
            new_halo = halo[:start_num_entries]
            new_halo.append(0)

        # append halo_info to new halo list
        new_halo_list.append(new_halo)
        # skip rest of tests
        continue
        
    # check that not a satellite 
    if isSatellite >= 0:
        # print reason
        print('IsSatellite')
        
        # ensure that running this code multiple times doesn't create a long list
        # check if halo already has more info
        if len(halo) > end_num_entries:
            # create copy and replace values
            new_halo = halo
            new_halo[start_num_entries] = 0
        else:
            # append 0's to halo_info if data not applicable
            new_halo = halo[:start_num_entries]
            new_halo.append(0)

        # append halo_info to new halo list
        new_halo_list.append(new_halo)
        # skip rest of tests
        continue
    
    # create a sphere data object with halo position and radius
    sp = ds.sphere(center, (radius.to('cm').value, 'cm'))
    
    # find stellar mass in second way 
    # find boolean mask for stellar particles
    stellar_mask = (sp[('all', 'particle_type')] == 2)
    mass_array = sp[('all', 'particle_mass')][stellar_mask] # array of masses
    x_array = sp[('all', 'particle_position_x')][stellar_mask] # arrays of position
    y_array = sp[('all', 'particle_position_y')][stellar_mask]
    z_array = sp[('all', 'particle_position_z')][stellar_mask]
    
    # find num values for halo center in cm to reduce calculations 
    xcenter = x.to('kpc').value
    ycenter = y.to('kpc').value
    zcenter = z.to('kpc').value
    
    # convert position arrays to kpc values
    x_array = x_array.in_units('kpc').value
    y_array = y_array.in_units('kpc').value
    z_array = z_array.in_units('kpc').value
    
    # find distances from center
    xdist = xcenter*np.ones(len(x_array)) - x_array
    ydist = ycenter*np.ones(len(y_array)) - y_array
    zdist = zcenter*np.ones(len(z_array)) - z_array
    netdistarray = (xdist**2 + ydist**2 + zdist**2)**0.5
    
    # create bins to check mass sum
    # in units kpc
    bins = np.geomspace(start=rad_min.to('kpc').value, stop=radius.to('kpc').value,
                        num=40)
    
    # check mass sum in each bin and find effective radius
    threshold = 0.5 * hstel_mass
    rad_eff = 0.0 * u.kpc
    stel_mass_eff = 0.0 * u.g
    for check_rad in bins:
        # no need to assign check_rad units, since all rad are in kpc
        
        # find stel mass at this radius
        mass_sum = 0
        
        # use bool mask to limit mass array
        check_limit_mask = netdistarray < check_rad
        mass_limit = mass_array[check_limit_mask]
        
        # find sum of masses
        mass_sum = mass_limit.sum()

        if mass_sum > threshold.to('g').value:
            print('break at: ', check_rad)
            rad_eff = check_rad * u.kpc
            stel_mass_eff = mass_sum.to('Msun')
            break
        else:
            rad_eff = check_rad * u.kpc
            stel_mass_eff = mass_sum.to('Msun')
    
    # the folowing block of code is to 
    # ensure that running this code multiple times doesn't create a long list
    # or eliminate further results
    
    # check if halos have further info
    if len(halo) > end_num_entries:
        # create copy with new info
        new_halo = halo
        new_halo[start_num_entries] = rad_eff.to('kpc')
    
    else:
        # append calculated ratios to truncated halo_info
        new_halo = halo[:start_num_entries]
        new_halo.append(rad_eff.to('kpc'))
            
    # print results
    print(index, radius, rad_eff, (stel_mass_eff/hstel_mass).value)
    
    # append halo_info to new halo list
    new_halo_list.append(new_halo)

# store new halo list to file
with open('calc_list_3000', 'wb') as outfile:
    pickle.dump(new_halo_list, outfile)
    
print('123AllDone123')
