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

# --find ratio of stellar mass to halo mass
# first list uses first output of TotalMass funct
# second list uses second output
ratiolist1 = []
ratiolist2 = []

for halo in hc.halo_list:
    # find coord and mass of halo
    x = halo.quantities.get('particle_position_x') * u.cm
    y = halo.quantities.get('particle_position_y') * u.cm
    z = halo.quantities.get('particle_position_z') * u.cm
    halo_mass = halo.quantities.get('particle_mass').in_units('Msun') * u.Msun
    
    # check if halo is inside zoom-in box
    if xmin <= x < xmax and ymin <= y < ymax and zmin <= z < zmax:
        # create a sphere data object with halo position and radius
        radius = halo.quantities.get('virial_radius') * u.cm
        sp = ds.sphere([x.value/scaling, y.value/scaling, z.value/scaling], (radius.value, 'cm'))
        
        # find the two output masses from TotalMass in Msun
        masses = sp.quantities.total_mass()
        gas_mass = masses[0] * u.g
        particle_mass = masses[1] * u.g
        gas_mass = gas_mass.to('Msun')
        particle_mass = particle_mass.to('Msun')
        
        # find stellar mass using total particle mass from TotalMass
        stellar_mass = particle_mass - halo_mass
        
        # find the two ratios
        ratio1 = gas_mass / halo_mass
        ratio2 = stellar_mass / halo_mass
        
        # add ratios to list
        ratiolist1.append(ratio1.value)
        ratiolist2.append(ratio2.value)
        
        
# store lists in files
with open('ratio_list0070_1.txt', 'wb') as ratiofile1:
    pickle.dump(ratiolist1, ratiofile1)

with open('ratio_list0070_2.txt', 'wb') as ratiofile2:
    pickle.dump(ratiolist2, ratiofile2)