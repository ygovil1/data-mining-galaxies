#specify input filenames
halo_catalog_filename = './halo_catalogs/catalog/catalog0070_thres160.0.h5'
dataset_filename = '~/../../tigress/cen/LG4_2048_40pc/RD0070/redshift0070'
redshift_filename = 'redshift0070'

# some constants
start_num_entries = 10 # number of entries halos have at start of this program
end_num_entries = 14 # min number of entries halos have at end of this program
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
# import halo catalogue func
from yt.analysis_modules.halo_analysis.api import *
# load halo dataset
halos_ds = yt.load(halo_catalog_filename)
# load raw dataset
ds = yt.load(dataset_filename)
# Instantiate a catalog using those two paramter files
hc = HaloCatalog(halos_ds=halos_ds, output_dir=os.path.join(tmpdir, 'halo_catalog'))
hc.load()

# load python halo list
with open('calc_list0070_2000', 'rb') as file1:
    halo_list = pickle.load(file1)

# load mass of smallest halo
with open('finest_particle0070', 'rb') as file2:
    min_mass = pickle.load(file2)
    min_mass = 100 * min_mass * u.g # 100*finest_particle

# load halo count from file
with open('count_halo0070_160', 'rb') as file3:
    count = pickle.load(file3)
    
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
# min = 1 kpc proper
# max = 0.5 Mpc comoving
# convert to centimeters value (without astropy units)
rad_min = 1 * u.kpc
rad_max = 0.5 * u.Mpc
rad_max = rad_max / (1 + redshift) # convert to physical

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

# --find ratio of stellar mass to halo mass
# first list uses first output of TotalMass funct
# second list uses second output
ratiolist1 = []
ratiolist2 = []

# update calc_halo list also
new_halo_list = []

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
    
    # check if halo is inside zoom-in box
    if xmin <= x < xmax and ymin <= y < ymax and zmin <= z < zmax:
        pass
    else:
        print('OutOfBounds')
        continue

    # check that radius is not 0
    if radius == 0:
        # print reason
        print('RadiusEqualZero')
        
        # check if halo already has more info
        if len(halo) > end_num_entries:
            # create copy and replace values
            new_halo = halo
            new_halo[start_num_entries] = 0
            new_halo[start_num_entries +1] = 0
            new_halo[start_num_entries +2] = 0
            new_halo[start_num_entries +3] = 0
            
            # append copy with new info
            new_halo_list.append(new_halo)
            # skip rest of tests
            continue
        
        # append 0's to halo_info if data not applicable
        # ensure that running this code multiple times doesn't ceate a long list
        new_halo = halo[:start_num_entries]
        new_halo.append(0)
        new_halo.append(0)
        new_halo.append(0)
        new_halo.append(0)

        # append halo_info to new halo list
        new_halo_list.append(new_halo)
        # skip rest of tests
        continue
        
    # check that not a satellite 
    if isSatellite >= 0:
        # print reason
        print('IsSatellite')
        
        # check if halo already has more info
        if len(halo) > end_num_entries:
            # create copy and replace values
            new_halo = halo
            new_halo[start_num_entries] = 0
            new_halo[start_num_entries +1] = 0
            new_halo[start_num_entries +2] = 0
            new_halo[start_num_entries +3] = 0
            
            # append copy with new info
            new_halo_list.append(new_halo)
            # skip rest of tests
            continue
        
        # append 0's to halo_info if data not applicable
        # ensure that running this code multiple times doesn't ceate a long list
        new_halo = halo[:start_num_entries]
        new_halo.append(0)
        new_halo.append(0)
        new_halo.append(0)
        new_halo.append(0)

        # append halo_info to new halo list
        new_halo_list.append(new_halo)
        # skip rest of tests
        continue
    
    # create a sphere data object with halo position and radius
    sp = ds.sphere(center, (radius.to('cm').value, 'cm'))

    # find the two output masses from TotalMass in Msun
    masses = sp.quantities.total_mass()
    gas_mass = masses[0] * u.g
    particle_mass = masses[1] * u.g
    gas_mass = gas_mass.to('Msun')
    particle_mass = particle_mass.to('Msun')

    # find stellar mass using total particle mass from TotalMass
    stellar_mass = particle_mass.to('Msun') - (halo_mass*omegas).to('Msun')
    
    # find stellar mass in second way 
    # find boolean mask for stellar particles
    stellar_mask = sp[('all', 'particle_type')] == 2
    mass_array = sp[('all', 'particle_mass')][stellar_mask] # array of masses
    # find stellar mass 2 by summing stellar particle masses
    stellar_mass2 = 0
    for imass in mass_array:
        stellar_mass2 = stellar_mass2 + imass
    # give units and convert stellar mass 2 to Msun
    stellar_mass2 = (stellar_mass2 * u.g).to('Msun')

    # find the two ratios
    ratio1 = gas_mass / halo_mass
    ratio2 = stellar_mass / halo_mass

    print(index, '\n', gas_mass.value, particle_mass.value, stellar_mass.value, stellar_mass2.value)
    print(ratio1, ratio2, '\n')
    
    # add ratios to list
    ratiolist1.append(ratio1.value)
    ratiolist2.append(ratio2.value)
    
    # the folowing block of code is to 
    # ensure that running this code multiple times doesn't ceate a long list
    
    # check if halos have further info
    if len(halo) > end_num_entries:
        # create copy with new info
        new_halo = halo
        new_halo[start_num_entries] = gas_mass.to('Msun')
        new_halo[start_num_entries +1] = particle_mass.to('Msun')
        new_halo[start_num_entries +2] = stellar_mass.to('Msun')
        new_halo[start_num_entries +3] = stellar_mass2.to('Msun')
    
    else:
        # append calculated ratios to truncated halo_info
        new_halo = halo[:start_num_entries]
        new_halo.append(gas_mass.to('Msun'))
        new_halo.append(particle_mass.to('Msun'))
        new_halo.append(stellar_mass.to('Msun'))
        new_halo.append(stellar_mass2.to('Msun'))
    
    # append halo_info to new halo list
    new_halo_list.append(new_halo)
        
# store ratio lists in files
with open('ratio_list0070_1.txt', 'wb') as ratiofile1:
    pickle.dump(ratiolist1, ratiofile1)
with open('ratio_list0070_2.txt', 'wb') as ratiofile2:
    pickle.dump(ratiolist2, ratiofile2)
    
# store new halo list to file
with open('calc_list0070_2000', 'wb') as outfile:
    pickle.dump(new_halo_list, outfile)