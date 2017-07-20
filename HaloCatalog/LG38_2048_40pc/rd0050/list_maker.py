#specify input filenames
halo_catalog_filename = './halo_catalogs/catalog/catalog.0.h5'
dataset_filename = '~/../../tigress/cen/LG38_2048_40pc/RD0050/redshift0050'

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

# create py list with halo quantities
halo_list = []

index = -1 #keep count of index
for halo in hc.halo_list:
    # increment index
    index += 1
      
    # find parameters of halo
    x = halo.quantities.get('particle_position_x') * u.cm
    y = halo.quantities.get('particle_position_y') * u.cm
    z = halo.quantities.get('particle_position_z') * u.cm
    mass = halo.quantities.get('particle_mass').in_units('Msun') * u.Msun
    radius = halo.quantities.get('virial_radius') * u.cm
    
    # create halo info tuple
    # store masses in Msun and lengths in pkpc
    halo_info = (index, x.to('kpc'), y.to('kpc'), z.to('kpc'), mass, radius.to('kpc'))
    
    # append info tuple
    halo_list.append(halo_info)
    
# sort list by mass
halo_list.sort(key=itemgetter(4), reverse = True)

# store list to file
with open('halo_list', 'wb') as outfile:
     pickle.dump(halo_list, outfile)

print('123AllDone123')
