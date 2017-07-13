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
halos_ds = yt.load('./halo_catalogs/catalog/catalog0070_thres160.0.h5')

# load raw dataset
ds = yt.load('~/../../tigress/cen/LG4_2048_40pc/RD0070/redshift0070')

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
with open('redshift0070', 'rt') as param_file:
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

# --find ratio of stellar mass to halo mass
# first list uses first output of TotalMass funct
# second list uses second output
ratiolist1 = []
ratiolist2 = []

for halo in halo_list:
    # find parameters of halo
    x = halo[2]
    y = halo[3]
    z = halo[4]
    center = [x.value/scaling, y.value/scaling, z.value/scaling]
    halo_mass = halo[5] # new mass
    radius = halo[7] # new radius
    
    # check if halo is inside zoom-in box
    if xmin <= x < xmax and ymin <= y < ymax and zmin <= z < zmax:
        pass
    else:
        continue

    # check that radius is not 0
    if radius == 0:
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
    stellar_mass = particle_mass.to('Msun') - halo_mass.to('Msun')

    # find the two ratios
    ratio1 = gas_mass / halo_mass
    ratio2 = stellar_mass / halo_mass

    print(gas_mass.value, particle_mass.value, stellar_mass.value)
    print(x, y, z, '\n')
    
    # add ratios to list
    ratiolist1.append(ratio1.value)
    ratiolist2.append(ratio2.value)

        
# store lists in files
with open('ratio_list0070_1.txt', 'wb') as ratiofile1:
    pickle.dump(ratiolist1, ratiofile1)

with open('ratio_list0070_2.txt', 'wb') as ratiofile2:
    pickle.dump(ratiolist2, ratiofile2)