# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 14:32:38 2018

@author: lindsayn
"""

import numpy as np
from pv_power import pv_power
from atmosphere import Atmosphere

def pv_prod(namelist, atm, panel, nt, output):
    """process called for multiprocessing to put result in Queue()"""
    P = pv_power(namelist, atm, panel, nt)    
    output.put(P)

def pv_interface_timeseries(namelist, atm, panel):
    """Return the spatio-temporal PV power from atmospheric inputs, panel characteristics and PV code options"""
    
    ncol, nt = np.shape(atm.ta)  # to get the dimensions of the input
    nthread = namelist.nthread   # number of parallel processes  
    power = []                   # initialization of the power to be returned 
    
    for n in range(ncol):
        # The code is called for each point of the domain, and a time series is returned
        atm_column = Atmosphere(*atm.get_column(n))                                 # extract one column from full atmosphere
        atm_column.set_sun_geometry(atm.sza[n,:], atm.saa[n,:], atm.sw_toa[n,:])    # sun geometry for the point is taken from full sun geometry
        
        if panel.technology == "tracker":
            panel.set_geometry(atm_column.sza, atm_column.saa, atm_column)          # 
            
        # Multi-porcessing    
        if nthread > 1:
            #----------------------------
            import multiprocessing as mp
            #----------------------------
            print("Multi-processes with %s processes"%nthread)
                                                     
            outputs = [mp.Queue() for x in range(nthread)]          # independent Queue objects to avoid out-of-order return at the end
            m = int(nt // nthread)                                  # how to equally slice the input atmosphere
            processes = []
            # Preparing all processes
            for k in range(nthread-1):
                n1 = k*m                                            # start of the slice
                n2 = (k+1)*m                                        # end of the slice
                atm_temp = Atmosphere(*atm_column.get_column_time_series(n1, n2))   # new atmosphere built as temporal slice of full atmosphere
                atm_temp.set_sun_geometry(atm_column.sza[n1:n2], atm_column.saa[n1:n2], atm_column.sw_toa[n1:n2])  # sun geometry set for this Atmosphere object
                processes+= [mp.Process(target = pv_prod, args = (namelist, atm_temp, panel, n2 - n1, outputs[k]))]    # processes call pv_production vi pv_prod interface and output the result to the dedicated Queue
                     
            # Residual for last thread
            n1 = (nthread-1)*m
            n2 = nt
            atm_temp = Atmosphere(*atm_column.get_column_time_series(n1, n2))
            atm_temp.set_sun_geometry(atm_column.sza[n1:n2], atm_column.saa[n1:n2], atm_column.sw_toa[n1:n2])
            processes+= [mp.Process(target=pv_prod, args=(namelist, atm_temp, panel, n2 - n1, outputs[-1]))]
            
            # Start processes
            for p in processes:
                p.start()

            # Exit the completed processes
            for p in processes:
                p.join()
                
            # Concatenate all Queues
            power0 = np.array([])
            for k,p in enumerate(processes):
                power0 = np.concatenate((power0, outputs[k].get()))                       
        
        # No multi-processing but splitting input according to maxcol if provided                                
        else:             
            try:
                mt = namelist.maxlen                        # split full timeseries into slices of length maxlen
                m = int(nt // mt )                          # this corresponds to m + 1 independent computations (nt = m * mt + res)
                power0 = np.array([])
                print("Time series split in", m)
                for k in range(m):
                    n1 = k*mt                                            # start of the slice
                    n2 = (k+1)*mt                                       # end of the slice
                    atm_temp = Atmosphere(*atm_column.get_column_time_series(n1, n2))   # new atmosphere built as temporal slice of full atmosphere
                    atm_temp.set_sun_geometry(atm_column.sza[n1:n2], atm_column.saa[n1:n2], atm_column.sw_toa[n1:n2])  # sun geometry set for this Atmosphere object
                    power0 = np.concatenate((power0, pv_power(namelist, atm_temp, panel, n2-n1)))
                    
                # Residual for last thread
                n1 = m*mt
                n2 = nt
                atm_temp = Atmosphere(*atm_column.get_column_time_series(n1, n2))
                atm_temp.set_sun_geometry(atm_column.sza[n1:n2], atm_column.saa[n1:n2], atm_column.sw_toa[n1:n2])
                power0 = np.concatenate((power0, pv_power(namelist, atm_temp, panel, n2-n1)))
                
                
            except:               
                power0 = pv_power(namelist, atm_column, panel, nt)
        
        power+= [power0] 
        
    return np.array(power)    
    
def pv_interface_atlas(namelist, atm, panel):
    """Return the spatio-temporal PV power from atmospheric inputs, panel characteristics and PV code options"""
    
    ncol, nt = np.shape(atm.ta)  # to get the dimensions of the input
    nthread = namelist.nthread   # number of parallel processes  
    power = []                   # initialization of the power to be returned 
    
    for n in range(nt): # loop over time, one atlas per timestep
        # The code is called for timestep n, and an atlas is returned
        atm_map = Atmosphere(*atm.get_map(n))                                    # extract one column from full atmosphere
        atm_map.set_sun_geometry(atm.sza[:,n], atm.saa[:,n], atm.sw_toa[:,n])    # sun geometry for all points in domain
        
        if panel.technology == "tracker":
            panel.set_geometry(atm_column.sza, atm_column.saa, atm_column)
            
        # Multi-porcessing    
        if nthread > 1:
            #----------------------------
            import multiprocessing as mp
            #----------------------------
            print("Multi-processes with %s processes"%nthread)
                                                     
            outputs = [mp.Queue() for x in range(nthread)]          # independent Queue objects to avoid out-of-order return at the end
            m = int(ncol // nthread)                                # how to equally slice the input atmosphere
            processes = []
            # Preparing all processes
            for k in range(nthread-1):
                n1 = k*m                                            # start of the slice
                n2 = (k+1)*m                                        # end of the slice
                atm_temp = Atmosphere(*atm_map.get_map_subdomain(n1, n2))   # new atmosphere built as temporal slice of full atmosphere
                atm_temp.set_sun_geometry(atm_map.sza[n1:n2], atm_map.saa[n1:n2], atm_map.sw_toa[n1:n2])  # sun geometry set for this Atmosphere object
                processes+= [mp.Process(target = pv_prod, args = (namelist, atm_temp, panel, n2 - n1, outputs[k]))]    # processes call pv_production vi pv_prod interface and output the result to the dedicated Queue
                     
            # Residual for last thread
            n1 = (nthread-1)*m
            n2 = ncol
            atm_temp = Atmosphere(*atm_map.get_map_subdomain(n1, n2))
            atm_temp.set_sun_geometry(atm_map.sza[n1:n2], atm_map.saa[n1:n2], atm_map.sw_toa[n1:n2])
            processes+= [mp.Process(target=pv_prod, args=(namelist, atm_temp, panel, n2 - n1, outputs[-1]))]
            
            # Start processes
            for p in processes:
                p.start()

            # Exit the completed processes
            for p in processes:
                p.join()
                
            # Concatenate all Queues
            power0 = np.array([])
            for k,p in enumerate(processes):
                power0 = np.concatenate((power0, outputs[k].get()))                       
        
        # No multi-processing but splitting input according to maxcol if provided                                
        else:             
            try:
                mt = namelist.maxlen                        # split full domain into slices of length maxlen
                m = int(ncol // mt )                        # this corresponds to m + 1 independent computations (nt = m * mt + res)
                power0 = np.array([])
                print("Time series split in", m)
                for k in range(m):
                    n1 = k*mt                                            # start of the slice
                    n2 = (k+1)*mt                                       # end of the slice
                    atm_temp = Atmosphere(*atm_map.get_map_subdomain(n1, n2))   # new atmosphere built as temporal slice of full atmosphere
                    atm_temp.set_sun_geometry(atm_map.sza[n1:n2], atm_map.saa[n1:n2], atm_map.sw_toa[n1:n2])  # sun geometry set for this Atmosphere object
                    power0 = np.concatenate((power0, pv_power(namelist, atm_temp, panel, n2-n1)))
                    
                # Residual for last thread
                n1 = m*mt
                n2 = ncol
                atm_temp = Atmosphere(*atm_map.get_map_subdomain(n1, n2))
                atm_temp.set_sun_geometry(atm_mapn.sza[n1:n2], atm_map.saa[n1:n2], atm_map.sw_toa[n1:n2])
                power0 = np.concatenate((power0, pv_power(namelist, atm_temp, panel, n2-n1)))
                
                
            except:               
                power0 = pv_power(namelist, atm_map, panel, ncol)
        
        power+= [power0] 
        
    return np.array(power).transpose()   # to keep ncol, nt dimensions