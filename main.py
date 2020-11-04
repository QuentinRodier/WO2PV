# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 14:32:38 2018

@author: lindsayn
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import tools
from atmosphere import Atmosphere
from namelist import Namelist # just create a file called namelist_yourname to import another Namelist
from panel import Panel
from pv_interface import pv_interface_timeseries, pv_interface_atlas # handles all inputs and calls pv_production
import time
import cartopy.crs as ccrs
from matplotlib.mlab import griddata

def get_gridded_data(lon,lat,var):
    """Interpolation of spatial data on a regular grid"""
    numcols, numrows = 400, 400
    xi = np.linspace(-1,4, numcols)
    yi = np.linspace(42,45, numrows)
    xi, yi = np.meshgrid(xi, yi)
    
    zi = griddata(lon,lat,var, xi, yi,interp='linear')
    
    return xi,yi,zi

# 1 Load options from namelist
def get_pv(forcing, output, n1 = 0, n2 = 10000, option = "timeseries"):
    """ forcing : nc file assumed to be located in Forcings/ e.g. "SIRTA_all" 
        output : name of output file containing PV production, stored in PV_Outputs directory       
        n1 and n2 : left and right indicess for partial computations along timeseries
        option : specify the type of forcing to optimize vectorial computations
    """
    print("----------Namelist loading----------------")
    nam = Namelist()                      # all options should be hardcoded in the Namelist() attributes
    
    print("----------Atmosphere loading----------------")   
    # build an Atmosphere object from an nc file. Such object can be built directly from np.arrays
    atm = Atmosphere(*tools.get_atm_from_ncfile(forcing = forcing, albedo = nam.albedo, n1 = n1, n2 = n2)) 
    atm.compute_sun()                     # to compute the full series of sza, saa, sw_toa
    
    # define PV panel from information contained in Namelist
    print("----------Module initialization----------------")
    panel = Panel(name = nam.panel_name, technology = nam.panel_technology, beta = nam.panel_beta, gamma = nam.panel_gamma, method_FF = nam.method_FF, SR_method = nam.SR_method, namelist = nam)       # get physical characteristics of the panel from a description file - here SIRTA panel used for consistency
                                
    # call pv_interface to handle large arrays and successive (and possibly multi-thread) calls to pv_production
    t1 = time.time()
    if option == "timeseries":
        power_simul = pv_interface_timeseries(nam, atm, panel)
    elif option == "atlas":
        power_simul = pv_interface_atlas(nam, atm, panel)
        
    # save PV outputs    
    np.savetxt("PV_Outputs/%s.dat"%output,power_simul)
    t2 = time.time()
    print("Time to run the PV code = %s seconds"%(t2-t1))
    
    return atm.dates, power_simul, atm.lat, atm.lon


if __name__ == "__main__":
    # The code below is used if directly calling this main.py file
    forcing = "atlas_AROME_subset" # atlas_AROME_subset, SIRTA_June, MNH etc.
    output = "atlas_AROME_subset"
    option = "atlas" # "timeseries", "atlas", controls how the data are splitted up (along time or space) for optimized computations
    n1 = 0           # only apply on subset along time (from n1 to n2 time steps)
    n2 = 1000
    data = sys.argv  # to read arguments if called via command-line
    
    # used if command-line call
    if len(data) > 3:       
        forcing = data[1]
        output = data[2]
        option = data[3]
        
        if len(data) > 4:
            n1 = int(data[4])
            n2 = int(data[5])
            dates, power_simul, lat, lon = get_pv(forcing, output, n1 , n2, option = option)
            
        else: 
            dates, power_simul, lat, lon = get_pv(forcing, output, option = option)
            
    else:
        dates, power_simul, lat, lon = get_pv(forcing, output, n1 , n2, option = option)
        
    
    # Showing the atlas or time series of PV power  
    if option == "atlas":
        "below works for forcing atlas_AROME_subset.nc" 
        lat2d,lon2d,power_simul2d = get_gridded_data(lon,lat,power_simul[:,0]) # interpolating on a regular lat/lon grid
          
        plot_lat_min = np.min(lat)
        plot_lat_max = np.max(lat)
        plot_lon_min = np.min(lon)
        plot_lon_max = np.max(lon)
  
        plt.subplots(1,1,figsize=(16,10))
        proj = ccrs.PlateCarree()
        ax = plt.axes(projection=proj)
        c = plt.pcolormesh(lat2d, lon2d, power_simul2d)
        ax.coastlines(resolution='10m', color='black')
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
        ax.set_extent((plot_lon_min-0.5, plot_lon_max+0.5, plot_lat_min-0.5, plot_lat_max+0.5))
        cb = plt.colorbar(c,pad=0.05)
        cb.set_label(label=r"W m$^{-2}$",size=22)
        plt.title("PV power production for SW France - 2020-06-06:12:00",size = 24,pad=25)
        plt.xlabel(r"Longitude ($^\circ$)",size=24)
        plt.ylabel(r"Latitude ($^\circ$)",size=24)
        
    else:    
        plt.figure(2,figsize=(18,10))
        plt.plot(dates, power_simul[0,:],"bo", label = "Simulation")
        plt.xlabel(r"Time",size=24)
        plt.ylabel(r"PV production (W m$^{-2}$)",size=24)
    
        if "SIRTA" in forcing:
            power_meas = np.loadtxt("Tests/%s_power_meas.dat"%forcing)
            plt.plot(dates, power_meas[n1:n2],"r", label = "Measurement")
        
    plt.legend(loc=0)
    plt.savefig("Figures/%s.jpg"%forcing,format = "jpg", dpi=200)
    plt.show()    


  
