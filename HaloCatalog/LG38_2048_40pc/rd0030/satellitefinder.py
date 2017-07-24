#specify input filenames
halo_catalog_filename = './halo_catalogs/catalog/catalog.0.h5'
dataset_filename = '~/../../tigress/cen/LG38_2048_40pc/RD0030/redshift0030'
redshift_filename = 'redshift0030'

# some constants
start_num_entries = 9 # number of entries halos have at start of this program
end_num_entries = 10 # min number of entries halos have at end of this program
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
    hubb_now = hubb_now.to('s**-1')
    
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
hubb_z = (100* hubb_now) * ((omega_m * (1 + redshift)**3) + (1 - omega_m))**0.5

# specify boundaries of zoom-in box
# scaling factor multiplied by info from text file 
# units in cm
scaling =  ((box_size / hubb_now.value) / (1 + redshift)).to('cm') # size of box
xmin = scaling*float(left_edges[0])
ymin = scaling*float(left_edges[1])
zmin = scaling*float(left_edges[2])
xmax = scaling*float(right_edges[0])
ymax = scaling*float(right_edges[1])
zmax = scaling*float(right_edges[2])

# set min mass to be 8e6 Msun
min_mass = 8e6 * u.Msun

# fill up lists
new_halo_list = []
for halo1 in halo_list:
    # find first index
    index1 = halo1[0]

    # create int to check if this halo isn't satellite
    isSatellite = -1

    # find parameters of halo1
    x1 = halo1[2]
    y1 = halo1[3]
    z1 = halo1[4]
    mass1 = halo1[5]
    radius1 = halo1[7]

    # loop through second halo
    for halo2 in halo_list:
        # find second index
        index2 = halo2[0]
        
        # if same halo, skip
        if index1 == index2:
            continue

        # find parameters of halo2
        x2 = halo2[2]
        y2 = halo2[3]
        z2 = halo2[4]
        mass2 = halo2[5]
        radius2 = halo2[7]
        
        # if error with new rad (too large, too small, didn't cross threshold) 
        # skip this halo
        if radius2 > (40 * u.kpc) or radius2 == 0:
            continue

        # calculate dist between halo centers
        xdist = x2 - x1
        ydist = y2 - y1
        zdist = z2 - z1
        netdist = (xdist**2 + ydist**2 + zdist**2)**0.5
        
        # If already a satellite of a halo, skip all future halos
        # this is because halos are ordered by decreasing mass
        if isSatellite >= 0:
            continue

        # if halo1 inside halo2, halo1 satellite of halo2
        if netdist < radius2:
            isSatellite = index2

    # after loop, append isSatellite int to halo_info
    # ensure that running this code multiple times doesn't ceate a long list
    # check if list already has more entries
    if len(halo_list) > end_num_entries:
        new_halo = halo1
        new_halo[start_num_entries] = isSatellite
    else:
        new_halo = halo1[:start_num_entries]
        new_halo.append(isSatellite)
    
    # print result
    print(new_halo)
    
    # append halo_info to new halo list
    new_halo_list.append(new_halo)

# store new list to file
with open('calc_list_3000', 'wb') as outfile:
    pickle.dump(new_halo_list, outfile)

print('123AllDone123')
