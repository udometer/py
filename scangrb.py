#!/usr/bin/env python
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.ndimage.filters import maximum_filter, minimum_filter
import numpy as np
import pdb
import pprint

import sys

print (sys.argv)
if (len(sys.argv) > 1):
    grib = sys.argv[1]
else:
    print ("Usage: ", sys.argv[0], " <grib file>")
    sys.exit(1)

grbs = pygrib.open(grib)

grb = grbs.select()
for r in grbs.select(): print (r.level, r.typeOfLevel, r.name, r.shortName, r.units)

pdb.set_trace()
pp = pprint.PrettyPrinter(indent=4, compact=True)

sn = set()
for k in grb:
    try:
        sn.add((k.shortName, k.name, k.units))
    except:
        continue

for i in sn:
    print (i)

rh = grbs.select(shortName='r')
levs = set()
for k in rh:
    levs.add((k.typeOfLevel, k.level))
rev_levs =  sorted((list(levs)))
rev_levs.reverse()
pp.pprint(rev_levs)

