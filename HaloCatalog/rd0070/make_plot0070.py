import yt
from yt.analysis_modules.halo_analysis.api import *
import tempfile
import shutil
import os

# Create temporary directory for storing files
tmpdir = tempfile.mkdtemp()

# create dataset
ds = yt.load('~/../../tigress/cen/LG4_2048_40pc/RD0070/redshift0070')

# Load the rockstar data files
halos_ds = yt.load('./halo_catalogs/catalog/catalog0070.0.h5')

# Instantiate a catalog using those two paramter files
hc = HaloCatalog(data_ds=ds, halos_ds=halos_ds, output_dir=os.path.join(tmpdir, 'halo_catalog'))

# create projection plot
p = yt.ProjectionPlot(ds, "x", "Dark_Matter_Density", width=(2, 'Mpc'))
p.set_zlim(field="Dark_Matter_Density", zmin=1e23, zmax=1e31)
p.annotate_halos(hc, factor = "particle_mass")
p.save("./projplot0070.png")
