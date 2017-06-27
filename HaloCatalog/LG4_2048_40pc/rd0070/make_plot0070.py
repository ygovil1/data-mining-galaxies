import yt
from yt.analysis_modules.halo_analysis.api import *

# create dataset
ds = yt.load('~/../../tigress/cen/LG4_2048_40pc/RD0070/redshift0070')

# Load the rockstar data files
halos_ds = yt.load('./halo_catalogs/catalog/catalog0070_thres160.0.h5')

# Instantiate a catalog using those two paramter files
hc = HaloCatalog(data_ds=ds, halos_ds=halos_ds)

# Filter out less massive halos
hc.add_filter("quantity_value", "particle_mass", ">", 1e8, "Msun")
hc.load()

# create projection plot
p = yt.ProjectionPlot(ds, "x", "Dark_Matter_Density", width=(1.5, 'Mpc'), max_level=30)
p.set_zlim(field="Dark_Matter_Density", zmin=5e24, zmax=2e26)
#p.annotate_halos(hc, factor = "particle_mass")
p.annotate_halos(hc)
p.save("./projplot0070.png")