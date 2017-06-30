# import basic libraries
import yt
import pickle
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

# import halo catalogue func
from yt.analysis_modules.halo_analysis.api import *

# import functions to ensure hc works
import tempfile
import shutil
import os

# Create temporary directory for storing files
tmpdir = tempfile.mkdtemp()

# load data dataset
# ds = yt.load('~/../../tigress/cen/LG4_2048_40pc/RD0070/redshift0070')

# load halo dataset
halos_ds = yt.load('./halo_catalogs/catalog/catalog0070_thres160.0.h5')

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

# create lists for satellites and non-satellites
satellite = []
nonsatellite = []

# fill up lists
for halo1 in hc.halo_list:
    # find first index
    index1 = halo1.quantities.get('particle_identifier')

    # create bool to check if this halo isn't satellite
    isSatellite = False

    # find parameters of halo1
    x1 = halo1.quantities.get('particle_position_x').in_units('cm')
    y1 = halo1.quantities.get('particle_position_y').in_units('cm')
    z1 = halo1.quantities.get('particle_position_z').in_units('cm')
    radius1 = halo1.quantities.get('virial_radius').in_units('cm')
    mass1 = halo1.quantities.get('particle_mass').in_units('Msun')

    # loop through second halo
    # checking if halo1 is satellite of halo2
    for halo2 in hc.halo_list:
        
        # if mass1 > mass2, skip other checks
        mass2 = halo2.quantities.get('particle_mass').in_units('Msun')
        if mass1 > mass2:
            continue
        
        # find second index
        index2 = halo1.quantities.get('particle_identifier')

        # find parameters of halo2
        x2 = halo2.quantities.get('particle_position_x').in_units('cm')
        y2 = halo2.quantities.get('particle_position_y').in_units('cm')
        z2 = halo2.quantities.get('particle_position_z').in_units('cm')
        radius2 = halo2.quantities.get('virial_radius').in_units('cm')
        
        # calculate dist between halo centers
        xdist = np.absolute(x2 - x1)
        ydist = np.absolute(y2 - y1)
        zdist = np.absolute(z2 - z1)
        netdist = (xdist**2 + ydist**2 + zdist**2)**0.5
        
        # if netdist < radius2, halo1 must be satellite of halo2
        if netdist < radius2:
            satellite.append(index1)
            isSatellite = True

    # after loop, if satellite bool still false, and halo1 not in satellite
    # list, then add to non-satellite list
    if not isSatellite and index1 not in satellite:
        nonsatellite.append(index1)

# store lists in files
with open('satellite_list0070v2.txt', 'wb') as satellite_file:
    pickle.dump(satellite, satellite_file)

with open('satellite_non_list0070v2.txt', 'wb') as nonsatellite_file:
    pickle.dump(nonsatellite, nonsatellite_file)
