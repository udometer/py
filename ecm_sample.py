#!/usr/bin/env python
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import numpy as np
import pdb

plt.figure(figsize=(12,8))
 
grib = 'jul18/gfs.t06z.pgrb2.0p25.f000.ctl'
grbs = pygrib.open(grib)
 
grb = grbs.select()[0]
data = grb.values
 
# need to shift data grid longitudes from (0..360) to (-180..180)
lons = np.linspace(float(grb['longitudeOfFirstGridPointInDegrees']), \
float(grb['longitudeOfLastGridPointInDegrees']), int(grb['Ni']) )
lats = np.linspace(float(grb['latitudeOfFirstGridPointInDegrees']), \
float(grb['latitudeOfLastGridPointInDegrees']), int(grb['Nj']) )
#data, lons = shiftgrid(180., data, lons, start=False)
grid_lon, grid_lat = np.meshgrid(lons, lats) #regularly spaced 2D grid

clat = (lats[0] + lats[1])/2
clon = (lons[0] + lons[1])/2

print (" " + clat + ", " + clon)

#reason for choosing Lambert conformal projection:
#https://deeppradhan.heliohost.org/gis/indian-grid/
#https://www.researchgate.net/post/What_is_the_projection_system_for_Survey_of_India_maps_How_can_I_set_that_projection_system_in_ArcGIS
ax = plt.axes(projection=ccrs.PlateCarree(
		central_longitude = clon,
	 )
	)
ax.coastlines()
plt.title('ECMWF 2017 07 16') # Set the name of the variable to plot
plt.show()
#plt.savefig(grib+'.png') # Set the output file name
#vim ts=4
