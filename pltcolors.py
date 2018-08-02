#!/usr/bin/env python
import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import colorcet as cc

ph = np.linspace(0, 2*np.pi, 13)
x = np.cos(ph)
y = np.sin(ph)
u = np.cos(ph)
v = np.sin(ph)
direction = ((270 - np.arctan2(u, v) * 180 /np.pi) %360)/90
print (u, v)
print (direction)

norm = Normalize(vmin=0, vmax=3)
# we need to normalize our colors array to match it colormap domain
# which is [0, 1]

colormap = cc.cm['diverging_bkr_55_10_c35'] #cm.inferno
plt.figure(figsize=(6, 6))
c=colormap(norm(direction))
pdb.set_trace()
plt.quiver(x, y, u, v, color=c)
plt.savefig("/tmp/test.png")
