#!/usr/bin/env python
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import Normalize
import colorcet as cc
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage.filters import maximum_filter, minimum_filter, gaussian_filter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pdb
import sys, os, math

print (sys.argv)

if (len(sys.argv) < 2):
    print ("Usage: ", sys.argv[0], " <grib file>")
    sys.exit(1)


knots_per_m = 1.94384
grib = sys.argv[1]
latitude = 11.28
longitude = 76.95
if (len(sys.argv) == 4):
    longitude = float(sys.argv[2])
    latitude = float(sys.argv[3])

myloc = os.path.dirname(sys.argv[0])

known_isobaric_levels = (1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1)
known_isobaric_levels_to_100 = (1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100)

class FiPA():
    def __init__(self, gribFile, context_degrees=1, subplot_rows=1, subplot_cols=1, xsize=19.2, ysize=10.8):
        self.grib = gribFile
        self.subrows = subplot_rows
        self.subcols = subplot_cols
        self.subidx = 1
        self.grbs = None
        self.IndiaShapeFile = os.path.join(myloc, 'Shapefiles', 'India.shp')

        self.grbs = pygrib.open(grib)
        self.context_degrees = context_degrees #degree

        prmsl = self.grbs.select(shortName='prmsl')[0]
        self.atime = prmsl.analDate.strftime("%Y-%m-%d %H%Mz")
        self.fh = prmsl.forecastTime

        self.suptitle = "(test FiPA small scale. 2° or 10° expanse) GFS for " + str(latitude) + ", " + str(longitude) + ". " +  self.atime +" analysis. Fcst hr " + str(self.fh)

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
        first_longitude = float(gribs['longitudeOfFirstGridPointInDegrees'])
        last_longitude = float(gribs['longitudeOfLastGridPointInDegrees'])
        Ni = gribs['Ni']

        first_latitude = float(gribs['latitudeOfFirstGridPointInDegrees'])
        last_latitude = float(gribs['latitudeOfLastGridPointInDegrees'])
        Nj = gribs['Nj']

        lons = list(np.linspace(first_longitude, last_longitude, Ni))
        lats = list(np.linspace(first_latitude, last_latitude, Nj))


        e = self.context_degrees
        def index(l, v):
            idx = -1
            i = 0
            while i < len(l):
                if (l[i] <= v):
                    idx = i
                else:
                    break
                i += 1

            return idx

        step_lat = step_lon = e
        (minlon, maxlon, minlat, maxlat) = (longitude-step_lon, longitude+step_lon, latitude-step_lat, latitude+step_lat)

        def get_extant_indices():
            x1 = index(lons, minlon)
            x2 = index(lons, maxlon)
            
            y1 = index(lats, minlat)
            y2 = index(lats, maxlat)
            return (x1, x2, y1, y2)

        (x1, x2, y1, y2) = get_extant_indices()

        if (e == 0 and (lons[x1] != longitude or lats[y1] != latitude)):
            import scipy
            #the location does not match exactly.
            #interpolate
            step_lon = (last_longitude - first_longitude)/(Ni-1)
            step_lat = (last_latitude - first_latitude)/(Nj-1)
            k = step_lon * 0
            minlon = (longitude // step_lon) * step_lon - k
            maxlon = minlon + step_lon + k
            
            minlat = (latitude//step_lat) * step_lat - k
            maxlat = minlat + step_lat + k

            (x1, x2, y1, y2) = get_extant_indices()

            #method recommended by matplotlib for 
            #regular sized rectangular grids.
            tlons = np.array(lons[x1:x2+1])
            tlats = np.array(lats[y1:y2+1]),
            tvals = values[y1:y2+1, x1:x2+1]

            f = scipy.interpolate.interp2d(tlons, tlats, tvals)
            nval = f([longitude], [latitude])

            lons = np.linspace(longitude, longitude, 0)
            lats = np.linspace(latitude, latitude, 0)
            x1=x2=y1=y2=0

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

    def trim_3d_data_to_extent(self, td): 
        d_all = []
        for d in td:
            grid_lon, grid_lat, dt = self.trim_array_to_extent(d, d.values)
            d_all.append((d.level, dt))
        d_all = sorted(d_all)
        d_all.reverse()
        ds = []
        for d in d_all:
            ds.append(d[1])
        return grid_lon, grid_lat, np.array(ds)
    
    def plot_skew_t(self):
        import FiPAskewt as skewt
        grbs = self.grbs

        self.context_degrees = 0;

        p = list(known_isobaric_levels_to_100)
        p = np.array(p)

        def get_1d_field(sn):
            f = grbs.select(shortName=sn, level=known_isobaric_levels_to_100, typeOfLevel='isobaricInhPa')
            grid_lon, grid_lat, fvalues = self.trim_3d_data_to_extent(f)
            f = fvalues.reshape(fvalues.size)
            return grid_lon, grid_lat, f

        grid_lon, grid_lat, RH = get_1d_field('r')
        grid_lon, grid_lat, Tc = get_1d_field('t')
        Tc -= 273.15
        grid_lon, grid_lat, z = get_1d_field('gh')
        grid_lon, grid_lat, u = get_1d_field('u')
        grid_lon, grid_lat, v = get_1d_field('v')

        #th =  (273.15 + Tc) * ( 1000 / p ) ** 0.286 + (3 * (RH * (3.884266 * 10 ** (( 7.5 * Tc ) / ( 237.7 + Tc )) ) /100 ))

        kappa_d = 0.2854131450730834
        th = (Tc+273.15)/((1000/p)**-kappa_d)

        def rh2mixr(RHp, PhPa, Tc):
            Es = 6.11*10**(7.5*Tc/(237.7+Tc))
            Ms = 0.622 * Es/(PhPa - Es)
            M= RHp*Ms/100
            return M #g/kg

        qv = rh2mixr(RH, p, Tc)

        #skewt.ptop = 150
        #skewt.pbottom = 1010
        self.suptitle = "(test FiPA) GFS " +  self.atime +" analysis. Fcst hr " + str(self.fh)
        skewt.plot(str(latitude) + ", " + str(longitude), z, th, p * 100, qv, u, v, grib + '_skewt_' + str(latitude) + "_" + str(longitude) + '.png', title=self.suptitle);

    def plot_pv(self, tempk, cint, context_degrees):
        self.context_degrees = context_degrees
        grbs = self.grbs

        t_in = grbs.select(shortName='t', level=known_isobaric_levels, typeOfLevel='isobaricInhPa')
        grid_lon, grid_lat, t_in = self.trim_3d_data_to_extent(t_in)

        u_in = grbs.select(shortName='u', level=known_isobaric_levels, typeOfLevel='isobaricInhPa')
        grid_lon, grid_lat, u_in = self.trim_3d_data_to_extent(u_in)

        v_in = grbs.select(shortName='v', level=known_isobaric_levels, typeOfLevel='isobaricInhPa')
        grid_lon, grid_lat, v_in = self.trim_3d_data_to_extent(v_in)

        hgt_in = grbs.select(shortName='gh', level=known_isobaric_levels, typeOfLevel='isobaricInhPa')
        grid_lon, grid_lat, hgt_in = self.trim_3d_data_to_extent(hgt_in)


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
        pi=math.pi

        tpdef=2

        lon1, lon2, lat1, lat2 = (grid_lon.min(), grid_lon.max(), grid_lat.min(), grid_lat.max())

        lat_in  = grid_lat[0:grid_lat.size-1,0:1].reshape(grid_lat[0].size)
        lon_in  = grid_lon[0][:]
        lev = list(known_isobaric_levels)
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

        # calculate PV on the tempkK isentropic surface
        # (also not in a pythonic way)
        nx=ix2-ix1+1
        ny=iy2-iy1+1
        nz=lev.size
        pvtempk=np.zeros((ny,nx))  # initialize as undef
        #utempk=np.zeros((ny,nx))  # initialize as undef
        #vtempk=np.zeros((ny,nx))  # initialize as undef
        #ptempk=np.zeros((ny,nx))  # initialize as undef
        for ix in range(0,nx):
            for iy in range(0,ny):
                for iz in range(nz-2,0,-1):#iz - increasing order of pressure
                    #print (theta[iz, iy, ix], theta[iz-1, iy, ix])
                    if (theta[iz,iy,ix]>=tempk and
                       theta[iz-1,iy,ix]<=tempk and
                       pvtempk[iy,ix] == 0):
                            ti = (iz, iy, ix)
                            ni = (iz-1, iy, ix)
                            ti2 = (iy, ix)
                            
                            pvtempk[ti2]=(
                            ((tempk-theta[ni])*pv[ti]+
                            (theta[ti]-tempk)*pv[ni])/
                            (theta[ti]-theta[ni])
                            )

                            """
                            #an average is probably not the most accurate calculation
                            utempk[ti2] = (u[ni] + u[ti])/2 * knots_per_m
                            vtempk[ti2] = (v[ni] + v[ti])/2 * knots_per_m
                            ptempk[ti2] = (lev[iz] + lev[iz-1])/2
                            """

        ax = self.newSubPlot()
        self.setup_grid(grid_lon.min(), grid_lon.max(), grid_lat.min(), grid_lat.max())
        clevs=np.arange(-10,11,1)
        shade_colors = 'RdBu_r'
        c = ax.contourf(lon, lat, pvtempk/1e-6, clevs, cmap=shade_colors, alpha=1, extend="both")
        plt.colorbar(c, ticks=clevs)

        """
        ax.barbs(lon, lat, utempk, vtempk, transform=ccrs.PlateCarree(), linewidth=0.5, regrid_shape=10, length=5)
        nlines = cint * 2
        c = ax.contour(lon, lat, ptempk, nlines, colors='k', alpha=0.75)
        """

        #############
        self.draw_boundaries()

        title = str(tempk) + "K PV, Wind (kt) "
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
        w = np.sqrt(np.square(du) + np.square(dv)) * knots_per_m

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
        ukvalues = uvalues * knots_per_m #convert from meters per second to knots
        vkvalues = vvalues * knots_per_m
        ax.streamplot(grid_lon, grid_lat, uvalues, vvalues, density=3, transform=ccrs.PlateCarree(), linewidth=0.5, arrowsize=0.00001, color=(0.5, 0.5, 0.5, 0.5))
        ax.barbs(grid_lon, grid_lat, ukvalues, vkvalues, transform=ccrs.PlateCarree(), linewidth=0.5, regrid_shape=10, length=5)
###############
        self.draw_boundaries()

        title = str(gh.level) + " hPa, gph(" + gh.units + " ), abs vort. (10^5 " + absv.units + ", shaded), wind (kt)"
        ax.set_title(title)
###############
        #plt.savefig(grib + str(gh.dataDate) + "_" + str(gh.dataTime) + "_" + str(gh.forecastTime) +  "_" + str(gh.shortName) + "_" + str(gh.level) + '.png') # Set the output file name


if __name__ == '__main__':
    xsize = 19.2
    ysize = 10.8
    #"""
    fipa = FiPA(grib, 1, 2, 3, xsize, ysize)
    fipa.fig = plt.figure(figsize=(xsize, ysize), clear=False, tight_layout=False)
    fipa.fig.suptitle(fipa.suptitle)
    fipa.plot_gph_vort_wind(850, 10, 1)
    fipa.plot_gph_vort_wind(500, 10, 1)
    fipa.plot_combined_rh_mslp(700, 400, 10, 5)
    fipa.plot_pv(355, 10, 5)
    fipa.plot_pv(330, 10, 5)
    fipa.plot_vertical_shear(200, 800, 10, 5)
    plt.savefig(grib + '_FiPA_' + str(latitude) + "_" + str(longitude) + '.png')
    plt.close()
    #"""
    fipa2 = FiPA(grib, 1, 1, 1, xsize, ysize)
    fipa2.plot_skew_t()

#vim ts=4
