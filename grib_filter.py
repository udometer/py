#!/usr/bin/env python
import requests, time, sys, os
from grib_filter_params import get_params as params
import pdb
def get_grib(dirname):
    (y, m, d, a, fs, l, r, t, b) = params()
    for f in fs:
        time.sleep(1)
        filename = os.path.join(dirname, "gfs_" + y + m + d + a + f + "_" + l + "_" + r + "_" + t + "_" + b)
        if (os.path.isfile(filename)):
            print ("File already exists: " + filename)
            return True
        #CMD="wget -t 0 --read-timeout=15 --timeout=15 -nc -O " + filename
        PRE="http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t" + a + "z.pgrb2.0p25.f" + f
        #OPT1="&lev_1000_mb=on&lev_100_m_above_ground=on&lev_100_mb=on&lev_10_m_above_ground=on&lev_200_mb=on&lev_300_mb=on&lev_400_mb=on&lev_500_mb=on&lev_600_mb=on&lev_700_mb=on&lev_750_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_950_mb=on&lev_mean_sea_level=on&lev_planetary_boundary_layer=on&lev_surface=on&var_DPT=on&var_HGT=on&var_POT=on&var_PRES=on&var_PRMSL=on&var_PWAT=on&var_RH=on&var_TMP=on&var_UGRD=on&var_VGRD=on"
        OPT="&all_lev=on&all_var=on"
        #OPT="&lev_1000_mb=on&lev_100_mb=on&lev_200_mb=on&lev_300_mb=on&lev_400_mb=on&lev_500_mb=on&lev_600_mb=on&lev_700_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_950_mb=on&lev_mean_sea_level=on&lev_planetary_boundary_layer=on&lev_PV%3D-2e-06_%5C%28Km%5C%5E2%2Fkg%2Fs%5C%29_surface=on&lev_PV%3D2e-06_%5C%28Km%5C%5E2%2Fkg%2Fs%5C%29_surface=on&lev_surface=on&var_HGT=on&var_POT=on&var_PRES=on&var_PRMSL=on&var_RH=on&var_TMP=on&var_UGRD=on&var_VGRD=on"
        REG="&subregion=&leftlon=" + l + "&rightlon=" + r + "&toplat=" + t + "&bottomlat=" + b
        DIR="&dir=%2Fgfs." + y + m + d + a
        URL = PRE + OPT + REG + DIR
        #cmd= CMD  + URL
        print("Getting file from: ", URL)
        
        response = requests.get(URL)
        if (response.status_code != 200):
            print ("Failed to get the file. HTTP request status code: ", response.status_code)
            return False

        with open(filename + ".part", 'wb') as fd:
            for chunk in response.iter_content(chunk_size=128):
                fd.write(chunk)

        os.rename(filename + ".part", filename)

        return True

if __name__ == "__main__":
    #print ("Please use as a module.")
    print (sys.argv)
    if (len(sys.argv) < 2):
        print ("Usage: ", sys.argv[0], " <grib file>")
        sys.exit(1)

    dirname = sys.argv[1]
    get_grib(dirname)
