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

known_isobaric_levels = (1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1)

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

        self.fig.suptitle("(test FiPA small scale. 2° or 10° expanse) GFS for " + str(latitude) + ", " + str(longitude)
            + ". " +  atime +" analysis. Fcst hr " + str(fh))

    def newSubPlot(self):
        self.ax = self.fig.add_subplot(self.subrows, self.subcols, self.subidx, projection=ccrs.PlateCarree())
        xh_pad = 0
        if (self.subidx > (self.subrows * self.subcols)/2):
            xh_pad = 0.25
        self.fig.tight_layout(rect=[0, 0.03, 1, 0.97], h_pad=1.25 + xh_pad)
        self.subidx += 1
        #self.ax = plt.axes(projection=ccrs.PlateCarree())
        return self.ax


    def trim_array_to_extent(self, gribs, values):
        lons = list(
                np.linspace(float(gribs['longitudeOfFirstGridPointInDegrees']), 
                    float(gribs['longitudeOfLastGridPointInDegrees']),
                    int(gribs['Ni'])))

        lats = list(np.linspace(float(gribs['latitudeOfFirstGridPointInDegrees']),
            float(gribs['latitudeOfLastGridPointInDegrees']),
            int(gribs['Nj'])))


        e = self.context_degrees
        (minlon, maxlon, minlat, maxlat) = (longitude-e, longitude+e, latitude-e, latitude+e)


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
        plt.plot(longitude, latitude, marker='o', markerfacecolor=None, color='red', markersize=8, alpha=0.8, 
            transform=ccrs.PlateCarree())

    def plot_pv(self, tempk, cint, context_degrees):
        self.context_degrees = context_degrees
        grbs = self.grbs

        def trim_3d_data_to_extent(td): 
            d_all = []
            for d in td:
                grid_lon, grid_lat, d = self.trim_array_to_extent(d, d.values)
                d_all.append(d)
            return grid_lon, grid_lat, np.array(d_all)
        
        t_in = grbs.select(shortName='t', level=known_isobaric_levels, typeOfLevel='isobaricInhPa')
        grid_lon, grid_lat, t_in = trim_3d_data_to_extent(t_in)

        u_in = grbs.select(shortName='u', level=known_isobaric_levels, typeOfLevel='isobaricInhPa')
        grid_lon, grid_lat, u_in = trim_3d_data_to_extent(u_in)

        v_in = grbs.select(shortName='v', level=known_isobaric_levels, typeOfLevel='isobaricInhPa')
        grid_lon, grid_lat, v_in = trim_3d_data_to_extent(v_in)

        hgt_in = grbs.select(shortName='gh', level=known_isobaric_levels, typeOfLevel='isobaricInhPa')
        grid_lon, grid_lat, hgt_in = trim_3d_data_to_extent(hgt_in)


        #############
        """
        Borrowed from 
        https://github.com/udometer/potential-vorticity/blob/master/gfs_pv_1.2.py
        """
        
        #constants
        re=6.37e6
        g=9.81
        cp=1004.5
        r=2*cp/7
        kap=r/cp
        omega=7.292e-5
        pi=3.14159265

        tpdef=2

        lon1, lon2, lat1, lat2 = (grid_lon.min(), grid_lon.max(), grid_lat.min(), grid_lat.max())

        lat_in  = grid_lat[0:grid_lat.size-1,0:1].reshape(grid_lat[0].size)
        lon_in  = grid_lon[0][:]
        lev = list(known_isobaric_levels)
        #lev.reverse()
        lev = np.array(lev)

        # get array indices for latitude-longitude range
        # specified above
        iy1 = np.argmin( np.abs( lat_in - lat1 ) )
        iy2 = np.argmin( np.abs( lat_in - lat2 ) )
        ix1 = np.argmin( np.abs( lon_in - lon1 ) )
        ix2 = np.argmin( np.abs( lon_in - lon2 ) )


        lon=lon_in[ix1:ix2+1]
        lat=lat_in[iy1:iy2+1]


        t=t_in[:,iy1:iy2+1,ix1:ix2+1]
        u=u_in[:,iy1:iy2+1,ix1:ix2+1]
        v=v_in[:,iy1:iy2+1,ix1:ix2+1]
        hgt=hgt_in[:,iy1:iy2+1,ix1:ix2+1]

        # some prep work for derivatives
        xlon,ylat=np.meshgrid(lon,lat)
        dlony,dlonx=np.gradient(xlon)
        dlaty,dlatx=np.gradient(ylat)
        dx=re*np.cos(ylat*pi/180)*dlonx*pi/180
        dy=re*dlaty*pi/180

        # define potential temperature and Coriolis parameter
        theta=t*(1.E5/(lev[:,np.newaxis,np.newaxis]*100))**kap
        f=2*omega*np.sin(ylat*pi/180)

        # calculate derivatives
        # (np.gradient can handle 1D uneven spacing,
        # so build that in for p, but do dx and dy 
        # external to the function since they are 2D)
        ddp_theta=np.gradient(theta,lev*100,axis=0)
        ddx_theta=np.gradient(theta,axis=2)/dx
        ddy_theta=np.gradient(theta,axis=1)/dy
        ddp_u=np.gradient(u,lev*100,axis=0)
        ddp_v=np.gradient(v,lev*100,axis=0)
        ddx_v=np.gradient(v,axis=2)/dx
        ddy_ucos=np.gradient(u*np.cos(ylat*pi/180),axis=1)/dy

        # calculate contributions to PV and PV
        absvort=ddx_v-(1/np.cos(ylat*pi/180))*ddy_ucos+f
        pv_one=g*absvort*(-ddp_theta)
        pv_two=g*(ddp_v*ddx_theta-ddp_u*ddy_theta)
        pv=pv_one+pv_two

        # calculate pressure of tropopause, Fortran-style (alas!)
        # as well as potential temperature (theta) and height
        #
        # starting from 10hPa and working down, to avoid
        # more complicated vertical structure higher up
        #
        nx=ix2-ix1+1
        ny=iy2-iy1+1
        nz=lev.size
        nzs=np.argwhere(lev==10.0)[0,0]
        tp=np.empty((ny,nx))*np.nan   # initialize as undef
        tp_theta=np.empty((ny,nx))*np.nan   # initialize as undef
        tp_hgt=np.empty((ny,nx))*np.nan   # initialize as undef

        for ix in range(0,nx):
            for iy in range(0,ny):
                for iz in range(nzs,0,-1):
                    if pv[iz,iy,ix]/1e-6<=tpdef:
                        if np.isnan(tp[iy,ix]):
                            tp[iy,ix]=(
                            (lev[iz]*(pv[iz+1,iy,ix]-tpdef*1e-6)
                            -lev[iz+1]*(pv[iz,iy,ix]-tpdef*1e-6))/
                            (pv[iz+1,iy,ix]-pv[iz,iy,ix])
                            )
            
                            tp_theta[iy,ix]=(
                            ((lev[iz]-tp[iy,ix])*theta[iz+1,iy,ix]+
                            (tp[iy,ix]-lev[iz+1])*theta[iz,iy,ix])/
                            (lev[iz]-lev[iz+1])
                            )
                            
                            tp_hgt[iy,ix]=(
                            ((lev[iz]-tp[iy,ix])*hgt[iz+1,iy,ix]+
                            (tp[iy,ix]-lev[iz+1])*hgt[iz,iy,ix])/
                            (lev[iz]-lev[iz+1])
                            )

        # calculate PV on the tempkK isentropic surface
        # (also not in a pythonic way)
        nx=ix2-ix1+1
        ny=iy2-iy1+1
        nz=lev.size
        pvtempk=np.zeros((ny,nx))  # initialize as undef
        for ix in range(0,nx):
            for iy in range(0,ny):
                for iz in range(nz-2,0,-1):
                    #print (theta[iz, iy, ix], theta[iz-1, iy, ix])
                    if theta[iz,iy,ix]>=tempk:
                        if theta[iz-1,iy,ix]<=tempk:
                            if pvtempk[iy,ix] == 0:
                                pvtempk[iy,ix]=(
                                ((tempk-theta[iz-1,iy,ix])*pv[iz,iy,ix]+
                                (theta[iz,iy,ix]-tempk)*pv[iz-1,iy,ix])/
                                (theta[iz,iy,ix]-theta[iz-1,iy,ix])
                                )

        ax = self.newSubPlot()
        self.setup_grid(grid_lon.min(), grid_lon.max(), grid_lat.min(), grid_lat.max())
        clevs=np.arange(-10,11,1)
        shade_colors = 'RdBu_r'
        c = ax.contourf(lon, lat, pvtempk/1e-6, clevs, cmap=shade_colors, alpha=1, extend="both")
        plt.colorbar(c, ticks=clevs)

        #############
        self.draw_boundaries()

        title = str(tempk) + "K PV "
        ax.set_title(title)


    def plot_combined_rh_mslp(self, levA, levB, cint, context_degrees):
        self.context_degrees = context_degrees
        grbs = self.grbs

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

        grid_lon, grid_lat, rh_all = self.trim_array_to_extent(grb, rh_all)

        ax = self.newSubPlot()
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
        clevs=np.arange(0,105,5)
        c = ax.contourf(grid_lon, grid_lat, rh_all, clevs, cmap=shade_colors.reversed(), 
                vmin=0, vmax=100, extend="both", alpha=1)
        plt.colorbar(c, ticks=clevs)

        mslp = grbs.select(shortName='prmsl')[0]
        grid_lon, grid_lat, mslpvalues = self.trim_array_to_extent(mslp, mslp.values)
        nlines = cint * 2 #int((mslp.max - mslp.min)/cint)
        c = ax.contour(grid_lon, grid_lat, mslpvalues, nlines, colors='k', linewidths=0.5, alpha=0.5)

        ######################
        """
        This section was adapted from:
        https://matplotlib.org/basemap/users/figures/plothighsandlows.py
        via
        https://matplotlib.org/basemap/users/examples.html
        """

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
        ######################

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
        clevs=np.arange(5,105,5)
        c = ax.contourf(grid_lon, grid_lat, w, clevs, cmap=shade_colors, alpha=1, vmin=5, extend="both")
        plt.colorbar(c, ticks=clevs)

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

        #print (absvvalues.min(), absvvalues.max())

        clevs=np.arange(5, 80,5)
        clevs=np.insert(clevs, 0, 3)
        c = ax.contourf(grid_lon, grid_lat, absvvalues * 100000, clevs, cmap=shade_colors, alpha=1, vmin=3, extend="both")
        plt.colorbar(c, ticks=clevs)
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

latitude = 11.28
longitude = 76.95
fipa = FiPA(grib, 1, 2, 3, 19.2, 10.8 )
#"""
fipa.plot_gph_vort_wind(850, 10, 1)
fipa.plot_gph_vort_wind(500, 10, 1)
fipa.plot_combined_rh_mslp(700, 400, 10, 5)
#"""
fipa.plot_pv(355, 10, 5)
fipa.plot_pv(330, 10, 5)
fipa.plot_vertical_shear(200, 800, 10, 5)
#fipa.plot_skew_t()

plt.savefig(grib + "_FiPA"+ '.') # Set the output file name
#vim ts=4
