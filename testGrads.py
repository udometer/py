#!/usr/bin/env python2
import pdb
import numpy as np
import matplotlib.pyplot as plt
import sys
from grads.ganum import GaNum
import sys

print (sys.argv)
if (len(sys.argv) > 1):
    grib = sys.argv[1]
else:
    print ("Usage: ", sys.argv[0], " <grib file>")
    sys.exit(1)

gn = GaNum(Bin='grads', Window=False)
fh = gn.open(grib)

gn("""
""")

gn("""
set lev 1000
""")
#ts = gn.exp("HGTprs")  # export variable ts from GrADS

gn("""
set mpdraw off
set gxout stream
d UGRDprs; VGRDprs
""")

gn("""
draw shp Shapefiles/ne_10m_coastline/ne_10m_coastline.shp
draw shp Shapefiles/India.shp
printim /tmp/test.png x1280 y720 white """
)
#vim ai ts=4 sw=4 sts=4 ft=python
