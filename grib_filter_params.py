#!/usr/bin/env python
def get_params():
    y="2018" #analysis year
    m="08"   #analysis month
    d="06"   #analysis date
    a="12"   #analysis hour
    fs=["006", ] #forecast hours

    lat=11.28   #central point latitude
    lon=76.95   #central point longitude
    ctxh=5      #context degrees in one direction. So context will be ctxh * 2
     
    l=lon - ctxh
    r=lon + ctxh
    t=lat + ctxh
    b=lat - ctxh

    lat = str(lat)
    lon = str(lon)
    l = str(l)
    r = str(r)
    t = str(t)
    b = str(b)
  
    params  = (y, m, d, a, fs, l, r, t, b)
    return params
