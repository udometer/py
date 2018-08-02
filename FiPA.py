#!/usr/bin/env python
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import Normalize
import colorcet as cc
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.ndimage.filters import maximum_filter, minimum_filter, gaussian_filter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pdb
import sys, os

print (sys.argv)

if (len(sys.argv) < 2):
    print ("Usage: ", sys.argv[0], " <grib file>")
    sys.exit(1)

grib = sys.argv[1]
myloc = os.path.dirname(sys.argv[0])

class FiPA():
    def __init__(self, gribFile, context_degrees=1, subplot_rows=1, subplot_cols=1, xsize=19.2, ysize=10.8):
        self.grib = gribFile
        self.subrows = subplot_rows
        self.subcols = subplot_cols
        self.subidx = 1
        self.grbs = None
        self.IndiaShapeFile = os.path.join(myloc, 'Shapefiles', 'India.shp')

        self.grbs = pygrib.open(grib)
        self.fig = plt.figure(figsize=(xsize, ysize), clear=False, tight_layout=False)
        self.context_degrees = context_degrees #degree

        prmsl = self.grbs.select(shortName='prmsl')[0]
        atime = prmsl.analDate.strftime("%Y-%m-%d %H%Mz")
        fh = prmsl.forecastTime

        self.fig.suptitle("(test FiPA small scale. 2° or 10° expanse) GFS for " + str(lat) + ", " + str(lon)
            + ". " +  atime +" analysis. Fcst hr " + str(fh))

    def newSubPlot(self):
        self.ax = self.fig.add_subplot(self.subrows, self.subcols, self.subidx, projection=ccrs.PlateCarree())
        self.subidx += 1
        self.fig.tight_layout(rect=[0, 0.03, 1, 0.97], h_pad=1.25)
        #self.ax = plt.axes(projection=ccrs.PlateCarree())
        return self.ax


    def trim_array_to_extent(self, gribs, values):
        #pdb.set_trace()
        lons = list(
                np.linspace(float(gribs['longitudeOfFirstGridPointInDegrees']), 
                    float(gribs['longitudeOfLastGridPointInDegrees']),
                    int(gribs['Ni'])))

        lats = list(np.linspace(float(gribs['latitudeOfFirstGridPointInDegrees']),
            float(gribs['latitudeOfLastGridPointInDegrees']),
            int(gribs['Nj'])))


        e = self.context_degrees
        (minlon, maxlon, minlat, maxlat) = (lon-e, lon+e, lat-e, lat+e)


        if 0:
            minlat = min(lats)
            maxlat = max(lats)
            minlon = min(lons)
            maxlon = max(lons)

        def index(l, v):
            idx = -1
            i = 0
            while i <= len(l):
                if (l[i] <= v):
                    idx = i
                else:
                    break
                i += 1

            return idx

        x1 = index(lons, minlon)
        x2 = index(lons, maxlon)

        y1 = index(lats, minlat)
        y2 = index(lats, maxlat)

        grid_lon, grid_lat = np.meshgrid(lons[x1:x2+1], lats[y1:y2+1]) #regularly spaced 2D grid
        pruned_values = []
        for row in values[y1:y2+1]:
            pruned_values.append(row[x1:x2+1])
           
        return (grid_lon, grid_lat, np.array(pruned_values))


    def setup_grid(self, minlon, maxlon, minlat, maxlat):
        e = self.context_degrees
        xrange = np.arange(minlon, maxlon, (maxlon - minlon + e)/10) 
        yrange = np.arange(minlat, maxlat, (maxlat - minlat + e)/10)

        if (self.ax):
            ax = self.ax
            ax.set_extent((minlon, maxlon, minlat, maxlat), crs=ccrs.PlateCarree())
            gl = ax.gridlines(draw_labels=False, crs=ccrs.PlateCarree(), linewidth=0.0,
                xlocs = xrange,
                ylocs = yrange
                )
            gl.xlabels_top = False
            gl.ylabels_right = False
            ax.set_xticks(xrange, crs=ccrs.PlateCarree())
            ax.set_yticks(yrange, crs=ccrs.PlateCarree())

    def draw_boundaries(self):
        self.ax.coastlines(resolution='10m')
        self.g = self.ax.add_geometries(shapereader.Reader(self.IndiaShapeFile).geometries(),
                          ccrs.PlateCarree(),
                          edgecolor='black',
                          facecolor='',
                          linewidth=0.25,
                          alpha=1
                          )
        self.g.linewidth=0.25
        plt.plot(lon, lat, marker='o', markerfacecolor=None, color='red', markersize=8, alpha=0.8, 
            transform=ccrs.PlateCarree())

    def plot_combined_rh_mslp(self, levA, levB, cint, context_degrees):
        self.context_degrees = context_degrees
        grbs = self.grbs

        known_isobaric_levels = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1]

        if (levA > levB):
            lev1 = levA
            lev2 = levB
        else:
            lev1 = levB
            lev2 = levA

        start_index = known_isobaric_levels.index(levA)
        end_index = known_isobaric_levels.index(levB)
        isobaric_levels_for_filter = known_isobaric_levels[start_index : end_index+1]

        rh_initialized = False
        rh_all = None
        rh_at_lev = None

        gribs_at_lev = grbs.select(shortName='r',typeOfLevel='isobaricInhPa', level=lambda l: l <= lev1 and l >= lev2, )
        for grb in gribs_at_lev:
            rh_at_lev = grb.values
            if (rh_initialized):
                rh_all += rh_at_lev
            else:
                rh_initialized = True
                rh_all = rh_at_lev

        if (not rh_initialized):
            return

        rh_all = rh_all/len(gribs_at_lev)

        ax = self.newSubPlot()

        grid_lon, grid_lat, rh_all = self.trim_array_to_extent(grb, rh_all)
        self.setup_grid(grid_lon.min(), grid_lon.max(), grid_lat.min(), grid_lat.max())

        shade_colors = LinearSegmentedColormap.from_list(
                name = "BrWGr",
                colors=[ #mimic tropicaltidbits.com
                "#0c352d",
                "#134a42", "#1a615b", "#297b75", "#41948e",
                "#62ada5", "#88c8c0", "#a9d8d2", "#e0eeed", "#f4f4f4",
                "#ffffff",
                "#f6ecde", "#f8e3c7", "#edcea5", "#e1b983", "#d19c5f", 
                "#c07c3b", "#a6642a", "#8b4c1c", "#6d3b14", "#512a0c"],
                N=21
                )
        c = ax.contourf(grid_lon, grid_lat, rh_all, 10, cmap=shade_colors.reversed(), 
                vmin=0, vmax=100, extend="both", alpha=1)
        plt.colorbar(c)

        mslp = grbs.select(shortName='prmsl')[0]
        grid_lon, grid_lat, mslpvalues = self.trim_array_to_extent(mslp, mslp.values)
        nlines = cint * 2 #int((mslp.max - mslp.min)/cint)
        c = ax.contour(grid_lon, grid_lat, mslpvalues, nlines, colors='k', linewidths=0.5, alpha=0.5)


        def extrema(mat, mode='nearest', size=3):
            mn = minimum_filter(mat, size=size, mode=mode)
            mx = maximum_filter(mat, size=size, mode=mode)
            return np.nonzero(mat == mn), np.nonzero(mat == mx)

        lmin, lmax = (extrema(mslpvalues/100, mode='nearest', size=20))
        xlows = grid_lon[lmin]
        ylows = grid_lat[lmin]

        xhighs = grid_lon[lmax]
        yhighs = grid_lat[lmax]

        lowvals = mslpvalues[lmin]
        highvals = mslpvalues[lmax]


        mxmax = grid_lon.max()
        mxmin = grid_lon.min()
        mymax = grid_lat.max()
        mymin = grid_lat.min()

        pdb.set_trace()

        xyplotted = []
        yoffset = 0.022 * (mymax - mymin)
        dmin = yoffset

        for x,y,p in zip(xlows, ylows, lowvals):
            if x < mxmax and x > mxmin and y < mymax and y > mymin:
                dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
                if not dist or min(dist) > dmin:
                    plt.text(x,y,'L',fontsize=14,fontweight='bold',
                            ha='center',va='center',color='b')
                    plt.text(x,y-yoffset,repr(int(p)),fontsize=9,
                            ha='center',va='top',color='b',
                            #bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5))
                            )
                    xyplotted.append((x,y))
        # plot highs as red H's, with max pressure value underneath.
        xyplotted = []
        for x,y,p in zip(xhighs, yhighs, highvals):
            if x < mxmax and x > mxmin and y < mymax and y > mymin:
                dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
                if not dist or min(dist) > dmin:
                    plt.text(x,y,'H',fontsize=14,fontweight='bold',
                            ha='center',va='center',color='r')
                    plt.text(x,y-yoffset,repr(int(p)),fontsize=9,
                            ha='center',va='top',color='r',
                            #bbox = dict(boxstyle="square", ec='None',fc=(1,1,1,0.5))
                            )
                    xyplotted.append((x,y))

        #plt.text(x, y, "L", transform=ccrs.PlateCarree())

        plt.clabel(c, fmt="%d")


        self.draw_boundaries()

        title = "MSLP (Pa) " + str(levA) + "-" + str(levB) + " hPa RH(%, shaded)"
        ax.set_title(title)

    def plot_vertical_shear(self, levA, levB, cint, context_degrees):
        #from https://www.wunderground.com/education/shear.asp
        #When one hears the phrase, "wind shear is 20 knots over the hurricane", this typically refers to the difference in wind speed between 200 mb (the top of the troposphere, 40,000 feet altitude) and a layer where a pressure of 850 mb is found--about 5,000 feet above the surface.

        self.context_degrees = context_degrees
        grbs = self.grbs

        u = grbs.select(shortName='u',typeOfLevel='isobaricInhPa',level=levA)[0]
        grid_lon, grid_lat, u1 = self.trim_array_to_extent(u, u.values)

        v = grbs.select(shortName='v',typeOfLevel='isobaricInhPa',level=levA)[0]
        grid_lon, grid_lat, v1 = self.trim_array_to_extent(v, v.values)

        u = grbs.select(shortName='u',typeOfLevel='isobaricInhPa',level=levB)[0]
        grid_lon, grid_lat, u2 = self.trim_array_to_extent(u, u.values)

        v = grbs.select(shortName='v',typeOfLevel='isobaricInhPa',level=levB)[0]
        grid_lon, grid_lat, v2 = self.trim_array_to_extent(v, v.values)


        #http://www.worldwindsinc.com/windshear.htm
        #http://www.drfitz.net/uploads/tips.pdf
        #vertical wind shear single point (VWSPT)
        #VWSPT = {[Du(r = 0)]^2 + [Dv(r = 0)]^2}1/2
        du = u1 - u2
        dv = v1 - v2
        w = np.sqrt(np.square(du) + np.square(dv)) * 1.94384

        nlines = cint #(w.max() - w.min()) / cint

        ax = self.newSubPlot()
        self.setup_grid(grid_lon.min(), grid_lon.max(), grid_lat.min(), grid_lat.max())

        clist = cc.rainbow_bgyr_35_85_c72
        clist.insert(0, '#ffffff');
        shade_colors = LinearSegmentedColormap.from_list(
            'whiteRainbow', clist)
        c = ax.contourf(grid_lon, grid_lat, w, nlines, cmap=shade_colors, alpha=1, vmin=5, extend="both")
        plt.colorbar(c)

        ax.streamplot(grid_lon, grid_lat, du, dv, density=1, transform=ccrs.PlateCarree(), linewidth=0.5, color=(0.0, 0.0, 0.0, 1))
        self.draw_boundaries()

        title = str(levA) + "hPa-" + str(levB) + "hpa Wind Shear (kt) (shaded & streamlined)"
        ax.set_title(title)


    def plot_gph_vort_wind(self, lev, cint, context_degrees):
        self.context_degrees = context_degrees
        grbs = self.grbs

        gh = grbs.select(shortName='gh',typeOfLevel='isobaricInhPa',level=lev)[0]
        grid_lon, grid_lat, ghvalues = self.trim_array_to_extent(gh, gh.values)

        r = grbs.select(shortName='r',typeOfLevel='isobaricInhPa',level=lev)[0]
        grid_lon, grid_lat, rvalues = self.trim_array_to_extent(r, r.values)

        u = grbs.select(shortName='u',typeOfLevel='isobaricInhPa',level=lev)[0]
        grid_lon, grid_lat, uvalues = self.trim_array_to_extent(u, u.values)

        v = grbs.select(shortName='v',typeOfLevel='isobaricInhPa',level=lev)[0]
        grid_lon, grid_lat, vvalues = self.trim_array_to_extent(v, v.values)

        absv = grbs.select(shortName='absv',typeOfLevel='isobaricInhPa',level=lev)[0]
        grid_lon, grid_lat, absvvalues = self.trim_array_to_extent(absv, absv.values)


        ax = self.newSubPlot()
        self.setup_grid(grid_lon.min(), grid_lon.max(), grid_lat.min(), grid_lat.max())
        
##############
        shade_colors = cc.cm['linear_kryw_0_100_c71_r']
        nlines = cint #int((absvvalues.max - absv.min)/cint)

        print (absvvalues.min(), absvvalues.max())

        c = ax.contourf(grid_lon, grid_lat, absvvalues * 100000, nlines, cmap=shade_colors, alpha=1, vmin=3, vmax=50, extend="both")
        plt.colorbar(c)
###############
        nlines = cint #int((gh.max - gh.min)/cint)
        c = ax.contour(grid_lon, grid_lat, ghvalues, nlines, cmap='hsv', linewidths=0.75, alpha=1)
        plt.clabel(c, fmt="%d")
###############
        #direction = ((((270 - np.arctan2(u.values, v.values) * 180 / np.pi) + 00) % 360)/90).astype(int)
        ukvalues = uvalues * 1.94384 #convert from meters per second to knots
        vkvalues = vvalues * 1.94384
        ax.streamplot(grid_lon, grid_lat, uvalues, vvalues, density=3, transform=ccrs.PlateCarree(), linewidth=0.5, arrowsize=0.00001, color=(0.5, 0.5, 0.5, 0.5))
        ax.barbs(grid_lon, grid_lat, ukvalues, vkvalues, transform=ccrs.PlateCarree(), linewidth=0.5, regrid_shape=10, length=5)
###############
        self.draw_boundaries()

        title = str(gh.level) + " hPa, gph(" + gh.units + " ), abs vort. (10^5 " + absv.units + ", shaded), wind (kt)"
        ax.set_title(title)
###############
        #plt.savefig(grib + str(gh.dataDate) + "_" + str(gh.dataTime) + "_" + str(gh.forecastTime) +  "_" + str(gh.shortName) + "_" + str(gh.level) + '.png') # Set the output file name

lat = 11.28
lon = 76.95
fipa = FiPA(grib, 1, 2, 3, 19.2, 10.8 )
fipa.plot_gph_vort_wind(850, 10, 1)
fipa.plot_gph_vort_wind(500, 10, 1)
fipa.plot_combined_rh_mslp(700, 400, 10, 5)
#fipa.plot_355k_pv()
fipa.plot_vertical_shear(200, 800, 10, 5)
#fipa.plot_skew_t()

plt.savefig(grib + "_FiPA"+ '.') # Set the output file name
#pdb.set_trace()
#vim ts=4
