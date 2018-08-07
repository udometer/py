#!/usr/bin/env bash
y=2018
m=08
d=06
a=12
fs="006"
lat=11.28
lon=76.95
ctxh=5
l=$(expr $lon - $ctxh)
r=$(expr $lon + $ctxh)
t=$(expr $lat + $ctxh)
b=$(expr $lat - $ctxh)

for f in $fs
do
sleep 1
CMD="wget -t 0 --read-timeout=15 --timeout=15 -nc -O gfs_"$y$m$d$a$f$l$r$t$b
PRE="http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t"$a"z.pgrb2.0p25.f"$f
#OPT1="&lev_1000_mb=on&lev_100_m_above_ground=on&lev_100_mb=on&lev_10_m_above_ground=on&lev_200_mb=on&lev_300_mb=on&lev_400_mb=on&lev_500_mb=on&lev_600_mb=on&lev_700_mb=on&lev_750_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_950_mb=on&lev_mean_sea_level=on&lev_planetary_boundary_layer=on&lev_surface=on&var_DPT=on&var_HGT=on&var_POT=on&var_PRES=on&var_PRMSL=on&var_PWAT=on&var_RH=on&var_TMP=on&var_UGRD=on&var_VGRD=on"
OPT="&all_lev=on&all_var=on"
#OPT="&lev_1000_mb=on&lev_100_mb=on&lev_200_mb=on&lev_300_mb=on&lev_400_mb=on&lev_500_mb=on&lev_600_mb=on&lev_700_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_950_mb=on&lev_mean_sea_level=on&lev_planetary_boundary_layer=on&lev_PV%3D-2e-06_%5C%28Km%5C%5E2%2Fkg%2Fs%5C%29_surface=on&lev_PV%3D2e-06_%5C%28Km%5C%5E2%2Fkg%2Fs%5C%29_surface=on&lev_surface=on&var_HGT=on&var_POT=on&var_PRES=on&var_PRMSL=on&var_RH=on&var_TMP=on&var_UGRD=on&var_VGRD=on"
REG="&subregion=&leftlon="$l"&rightlon="$r"&toplat="$t"&bottomlat="$b
DIR="&dir=%2Fgfs."$y$m$d$a
cmd="$CMD $PRE$OPT$REG$DIR"
echo $cmd
$cmd
sleep 10
done
