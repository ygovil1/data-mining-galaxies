#specify input filenames
halo_catalog_filename = './halo_catalogs/catalog/catalog.0.h5'
dataset_filename = '~/../../tigress/cen/LG4_2048_40pc/RD0110/redshift0110'
redshift_filename = 'redshift0110'

# some constants
start_num_entries = 13 # number of entries halos have at start of this program
end_num_entries = 14 # min number of entries halos have at end of this program
# note- start_num_entires is index of first changed/appended value

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

# specify boundary on neighborhood
# chosen to be 1 comoving Mpc
nbhood = 1 * u.Mpc
nbhood = nbhood / (1 + redshift) # convert to physical

# store new halo info
new_halo_list = []

# loop through halos to find prox indicator
for halo1 in halo_list:
    # create prox value and sum of neighborhood mass
    M_other = 0 * u.Msun
    eta = 0 
    
    # find parameters of halo
    index1 = halo1[0]
    x1 = halo1[2]
    y1 = halo1[3]
    z1 = halo1[4]
    center1 = [(x1/scaling).value, (y1/scaling).value, (z1/scaling).value]
    halo_mass1 = halo1[5] # new mass
    radius1 = halo1[7] # new radius
    isSatellite1 = halo1[9]
    
    # skip halo if outside zoom box or satellite
    # check if halo is inside zoom-in box
    if xmin <= x1 < xmax and ymin <= y1 < ymax and zmin <= z1 < zmax:
        pass
    else:
        print('OutOfBounds')
        continue

    # check that mass is not 0
    if halo_mass1 == 0:
        # print reason
        print('HaloMassEqualZero')
        
        # ensure that running this code multiple times doesn't create a long list
        # check if halo already has more info
        if len(halo1) > end_num_entries:
            # create copy and replace values
            new_halo = halo1
            new_halo[start_num_entries] = 0
        else:
            # append 0's to halo_info if data not applicable
            new_halo = halo1[:start_num_entries]
            new_halo.append(0)

        # append halo_info to new halo list
        new_halo_list.append(new_halo)
        # skip rest of tests
        continue
        
    # check that not a satellite 
    if isSatellite1 >= 0:
        # print reason
        print('IsSatellite')
        
        # ensure that running this code multiple times doesn't create a long list
        # check if halo already has more info
        if len(halo1) > end_num_entries:
            # create copy and replace values
            new_halo = halo1
            new_halo[start_num_entries] = 0
        else:
            # append 0's to halo_info if data not applicable
            new_halo = halo1[:start_num_entries]
            new_halo.append(0)

        # append halo_info to new halo list
        new_halo_list.append(new_halo)
        # skip rest of tests
        continue
    
    # loop through other halos
    for halo2 in halo_list:
        # find parameters of halo
        index2 = halo2[0]
        x2 = halo2[2]
        y2 = halo2[3]
        z2 = halo2[4]
        center2 = [(x2/scaling).value, (y2/scaling).value, (z2/scaling).value]
        halo_mass2 = halo2[5] # new mass
        radius2 = halo2[7] # new radius
        isSatellite2 = halo2[9]
        
        # skip same halo again
        if index1 == index2:
            continue
        
        # find distance between halo centers
        xdist = x2 - x1
        ydist = y2 - y1
        zdist = z2 - z1
        netdist = (xdist**2 + ydist**2 + zdist**2)**0.5
        
        # check if halo2 is in neighborhood
        # if so, add mass2 to M_other
        if netdist < nbhood:
            M_other += halo_mass2
    
    # calculate eta = M_other / M_h
    # set to zero if halo_mass is zero
    if halo_mass1 > 0:
        eta = M_other.to('Msun') / halo_mass1
    else:
        eta = 0 * u.Msun / u.Msun
    
    # ensure eta is stored in proper location in list
    # check if halos have further info
    if len(halo1) > end_num_entries:
        # create copy with new info
        new_halo = halo1
        new_halo[start_num_entries] = eta.value
    else:
        # append calculated ratios to truncated halo_info
        new_halo = halo1[:start_num_entries]
        new_halo.append(eta.value)
    
    print(index1, eta)
    
    # append halo_info to new halo list
    new_halo_list.append(new_halo)

# store new halo list to file
with open('calc_list_3000', 'wb') as outfile:
    pickle.dump(new_halo_list, outfile)
    
print('123AllDone123')
