# import basic libraries
import yt
import numpy as np
import matplotlib.pyplot as plt

# import halo catalogue func
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from yt.analysis_modules.halo_finding.halo_objects import HaloFinder

# create dataset
ds = yt.load('~/../../tigress/cen/LG38_2048_40pc/RD0050/redshift0050')

# specify parameters for halo finder
finder_pars = {'threshold': 160., 'dm_only':True, 'padding' : 0.02}

# create halo catalogue
hc = HaloCatalog(finder_kwargs=finder_pars, data_ds=ds, finder_method='hop')
hc.create()
hc.save_catalog()
