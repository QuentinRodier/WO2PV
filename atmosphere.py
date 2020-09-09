# -*- coding: utf-8 -*-

import numpy as np
import sys
sys.path.append("Tools")
from tools import get_sun

class Atmosphere():   
    def __init__(self, lat, lon, dates, swdir, swdiff, ta, wind, albdir = None, albdiff = None):          
        self.lat = lat
        self.lon = lon
        self.dates = dates
        self.nt = len(np.atleast_1d(dates)) # to ensure it works for single time steps
        self.swdir = swdir
        self.swdiff = swdiff
        self.ta = ta
        self.wind = wind   
        
        if not albdir is None:              # albedo initialized only if provided (optional since can be external)
            self.albdir = albdir
            self.albdiff = albdiff
        else:
            self.albdir = np.zeros_like(swdir)
            self.albdiff = np.zeros_like(swdir)
    
    def set_sun_geometry(self, sza, saa, sw_toa):
        """Fill the sun variables directly without computing again"""
        self.sza = sza
        self.saa = saa
        self.sw_toa = sw_toa
        
    def compute_sun(self): 
        """Compute solar zenith and azimuth angles,
        as well as local solar constant from lat, lon, date"""
        sza = np.zeros_like(self.ta)
        saa = np.zeros_like(self.ta)
        sw_toa = np.zeros_like(self.ta)
        ncol, nt = np.shape(self.ta)
        for i in range(ncol):
            for j in range(nt):
                sza[i,j], saa[i,j], sw_toa[i,j] = get_sun(self.lat[i], self.lon[i], self.dates[j])

        self.sza = sza
        self.saa = saa
        self.sw_toa = sw_toa
        
    def get_column(self, ic):
        """Extract one column (index ic) out of a spatio-temporal series"""
        return self.lat[ic], self.lon[ic], self.dates, self.swdir[ic,:,:], self.swdiff[ic,:,:], self.ta[ic,:], self.wind[ic,:], self.albdir[ic,:,:], self.albdiff[ic,:,:]
    
    def get_column_time_series(self, n1, n2):
        """Extract a time series out of a single column full time series"""
        return self.lat, self.lon, self.dates[n1:n2], self.swdir[n1:n2,:], self.swdiff[n1:n2,:], self.ta[n1:n2], self.wind[n1:n2], self.albdir[n1:n2,:], self.albdiff[n1:n2,:]
        
    def get_map(self, it):
        """Extract one time step out of a spatio-temporal series"""    
        return self.lat, self.lon, self.dates[it], self.swdir[:,it,:], self.swdiff[:,it,:], self.ta[:,it], self.wind[:,it], self.albdir[:,it,:], self.albdiff[:,it,:]
    
    def get_map_subdomain(self, n1, n2):
        """Extract a subdomain out of a full domain"""    
        return self.lat[n1:n2], self.lon[n1:n2], self.dates, self.swdir[n1:n2,:], self.swdiff[n1:n2,:], self.ta[n1:n2], self.wind[n1:n2], self.albdir[n1:n2,:], self.albdiff[n1:n2,:]
    
    

            
            
    