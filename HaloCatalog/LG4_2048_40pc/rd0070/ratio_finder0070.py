#specify input filenames
halo_catalog_filename = './halo_catalogs/catalog/catalog0070_thres160.0.h5'
dataset_filename = '~/../../tigress/cen/LG4_2048_40pc/RD0070/redshift0070'
redshift_filename = 'redshift0070'

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
        continue

    # check that radius is not 0
    if radius == 0:
        # append 0's to halo_info if data not applicable
        # ensure that running this code multiple times doesn't ceate a long list
        new_halo = halo[:10]
        new_halo.append(0)
        new_halo.append(0)

        # append halo_info to new halo list
        new_halo_list.append(new_halo)
        # skip rest of tests
        continue
        
    # check that not a satellite 
    if isSatellite >= 0:
        # append 0's to halo_info if data not applicable
        # ensure that running this code multiple times doesn't ceate a long list
        new_halo = halo[:10]
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

    # find the two ratios
    ratio1 = gas_mass / halo_mass
    ratio2 = stellar_mass / halo_mass

    print(index, '\n', gas_mass.value, particle_mass.value, stellar_mass.value)
    print(ratio1, ratio2, '\n')
    
    # add ratios to list
    ratiolist1.append(ratio1.value)
    ratiolist2.append(ratio2.value)
    
    # append calculated ratios to halo_info
    # ensure that running this code multiple times doesn't ceate a long list
    new_halo = halo[:10]
    new_halo.append(gas_mass.to('Msun'))
    new_halo.append(particle_mass.to('Msun'))
    new_halo.append(stellar_mass.to('Msun'))
    
    # append halo_info to new halo list
    new_halo_list.append(new_halo)
        
# store ratio lists in files
with open('ratio_list0070_1.txt', 'wb') as ratiofile1:
    pickle.dump(ratiolist1, ratiofile1)
with open('ratio_list0070_2.txt', 'wb') as ratiofile2:
    pickle.dump(ratiolist2, ratiofile2)
    
# store new halo list to file
with open('calc_list00070_2000', 'wb') as outfile:
    pickle.dump(new_halo_list, outfile)