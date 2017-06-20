# import basic libraries
import yt
import numpy as np
import matplotlib.pyplot as plt

# import halo catalogue func
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from yt.analysis_modules.halo_finding.halo_objects import HaloFinder

# create dataset
ds = yt.load('~/../../tigress/cen/LG4_2048_40pc/redshift0110')

# specify parameters for halo finder
finder_pars = {'threshold': 160., 'dm_only':False, 'ptype':'io', 'padding': 0.02}

# create halo catalogue
hc = HaloCatalog(finder_kwargs=finder_pars, data_ds=ds, finder_method='hop', output_dir='/tigress/ygovil/Halos/io')
hc.create()
