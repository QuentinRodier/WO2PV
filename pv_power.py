# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:05:58 2018

@author: lindsayn
"""

import numpy as np
import os
import scipy.interpolate
import scipy.integrate
import tools
import matplotlib.pyplot as plt

def pv_power(namelist, atm, panel, nt, debug = 0):   
    """Return power time series from atmospheric time series"""
    
    # Checking spectral resolution of irradiance input 
    nswband = np.shape(atm.swdir)[1] 
    print("Number of spectral bands of inputs: %s"%nswband)
    
    if nswband < 1:
        print('Inappropriate number of bands')
    
    if len(namelist.swbands) != nswband + 1:     # namelist.swbands gives the limits of the bands
        print("WARNING: Spectral bands inconsistent with irradaince inputs")
                
    #======================================   SPECTRALISATION    =================================================
    # From broadband or narrowband inputs compute high spectral resolution spectra using ARTDECO reference spectra
    #=============================================================================================================
    
    swdiff_hres = []       # Diffuse irradiance high spectral resolution
    swdir_hres = []        # Direct irradiance high spectral resolution
    swref_hres = []        # Reflected irradiance high spectral resolution
          
    wvl_artdeco = tools.wvl_artdeco()  # ARTDECO wavelengths of reference spectra
    
    if namelist.albedo == "weather":                         # albedo read from Atmosphere so nominal resolution
        albdir_hres = np.zeros([nt, len(wvl_artdeco)])
        albdiff_hres = np.zeros([nt, len(wvl_artdeco)])
        for k, wvl_inf in enumerate(namelist.swbands[:-1]):
            wvl_sup = namelist.swbands[k+1]
            i, j = np.searchsorted(wvl_artdeco, [wvl_inf, wvl_sup])
            albdir_hres[:,i:j] = np.transpose(np.tile(atm.albdir[:,k],(j-i,1)))
            albdiff_hres[:,i:j] = np.transpose(np.tile(atm.albdiff[:,k],(j-i,1)))
                    
    else:
        # Interpolation already performed on ARTDECO resolution when initializing Atmosphere
        albdir_hres = atm.albdir
        albdiff_hres = atm.albdiff
    
    if debug > 0:         # extra runtime outputs to be used for bedugging
        plt.step(namelist.swbands[:-1], atm.albdir[0,:])    
        plt.plot(wvl_artdeco, albdir_hres[0,:]) 
        plt.show()
    
    for i in range(nt): 
        # spectralisationis performed iteratively for each time step, this takes time   
        if nswband == 1:
            swdir_hres0 = tools.hres(atm.swdir[i], bands = namelist.swbands,sza = atm.sza[i],method = 'single-band', source = "direct")
            swdir_hres+= [swdir_hres0]
            swdiff_hres0 = tools.hres(atm.swdiff[i], bands = namelist.swbands,sza = atm.sza[i],method = 'single-band', source = "diffuse")
            swdiff_hres+= [swdiff_hres0]
            swref_hres+= [albdiff_hres[i,:]*swdiff_hres0 + albdir_hres[i,:]*swdir_hres0]
                    
        elif nswband > 1:
            swdir_hres0 = tools.hres(atm.swdir[i,:], bands = namelist.swbands, sza = atm.sza[i],method = 'multi-band', source = "direct")
            swdir_hres+= [swdir_hres0]
            swdiff_hres0 = tools.hres(atm.swdiff[i,:], bands = namelist.swbands, sza = atm.sza[i],method = 'multi-band', source = "diffuse")
            swdiff_hres+= [swdiff_hres0]
            swref_hres+= [albdiff_hres[i,:]*swdiff_hres0 + albdir_hres[i,:]*swdir_hres0]
  
    swdir_hres = np.array(swdir_hres)
    swdiff_hres = np.array(swdiff_hres)
    swref_hres = np.array(swref_hres)
    
    # Compute broadband incident radiation
    swdir_bb = np.sum(atm.swdir,axis=1)
    swdiff_bb = np.sum(atm.swdiff,axis=1)
    
    #=============================  TRANSPOSITION  ================================
    # How to get irradiance in the panel plane ?
    #==============================================================================
    
    # Direct irradiance
    # Angle of incidence of direct irradiance on the panel
    theta = tools.incident_angle_panel(panel.beta, panel.gamma, atm.sza, atm.saa)
    
    # Transposition coefficients
    if namelist.method_POA == 'default':
        # Transposition method different for different components, choice made based on analysis of SIRTA data
        transpose_dir, transpose_diff, transpose_ref = tools.transposition(theta, atm.sza, panel.beta, swdir_bb = swdir_bb, swdiff_bb = swdiff_bb, sw_bb = swdir_bb + swdiff_bb, sw_toa = atm.sw_toa)
    
    else:        
        transpose_dir = tools.transposition_direct(theta, atm.sza)
        transpose_diff = tools.transposition_diffuse(panel.beta, method = namelist.method_POA, theta = theta, sza = atm.sza, swdir_bb = swdir_bb, swdiff_bb = swdiff_bb, sw_bb = swdir_bb + swdiff_bb, sw_toa = atm.sw_toa)
        # Reflected forced to be isotropic if not isotropic3D
        transpose_ref = tools.transposition_reflected(panel.beta, method = 'isotropic'*(namelist.method_POA != 'isotropic') + 'isotropic'*(namelist.method_POA == 'isotropic'))
    
    # Plane-of-array (POA) irradiance
    swdir_poa_hres, swdiff_poa_hres, swref_poa_hres = tools.sw_poa(swdir_hres, swdiff_hres, swref_hres,transpose_dir, transpose_diff, transpose_ref)
    sw_poa_hres = swdir_poa_hres + swdiff_poa_hres + swref_poa_hres  # High resolution POA
    sw_poa_bb = scipy.integrate.simps(sw_poa_hres, wvl_artdeco)
    
#    plt.figure(1)
#    plt.plot(sw_poa_hres[200,:])
#    plt.show()
   
#    swdir_poa_bb = scipy.integrate.simps(swdir_poa_hres, wvl_artdeco)
#    swdiff_poa_bb = scipy.integrate.simps(swdiff_poa_hres, wvl_artdeco)
#    swref_poa_bb = scipy.integrate.simps(swref_poa_hres, wvl_artdeco)    
    
    # Optical losses - not spectral yet
    ol_dir, ol_diff, ol_ref = tools.optical_losses(theta, panel.beta, method = namelist.method_OL, ar = namelist.ar, c1 = namelist.c1, c2 = namelist.c2, n0 = namelist.n0, n1 = namelist.n1, K = namelist.K, L = namelist.L)
    
    #Cell-incident irradiance
    swdir_cell_hres, swdiff_cell_hres, swref_cell_hres = tools.sw_cell(swdir_poa_hres, swdiff_poa_hres, swref_poa_hres, ol_dir, ol_diff, ol_ref)  # multiplying by OL to account for optical losses
    sw_cell_hres = swdir_cell_hres + swdiff_cell_hres + swref_cell_hres
    
    
    #========================  CELL TEMPERATURE  ==================================
    # How to get cell temperature from air temperature ? 
    #==============================================================================   

    if namelist.method_tc == 'measured':# cell temperature directly measured
        tc = panel.tcell
        
    elif namelist.method_tc == 'measuredKing2004': # module temperature measured, converted to cell temperature
        # Correction following King  2004       
        tc = tools.mod2cell_temperature(panel.tcell, sw_poa_bb, mount_type = namelist.mount_type_King)
       
    else: 
        # Computed from weather parameters
        # Methos available: Simple, King2004, Skoplaki2008_1, Skoplaki2008_2, Skoplaki2008_3, Skoplaki2008_4
        tc = tools.cell_temperature(atm.ta, sw_poa_bb, method = namelist.method_tc, NOCT = panel.NOCT, ta_NOCT = panel.ta_NOCT, sw_NOCT = panel.sw_NOCT, CtPmax = panel.Ct_Pmax, eta_STC = panel.Pmax_STC/panel.sw_STC, \
                                   wind = atm.wind, mount_type = namelist.mount_type_King, tau_alpha = namelist.tau_alpha, ta_STC = panel.ta_STC)

    #==============================================================================
    #=============   OPERATING ONE-DIODE MODEL CIRCUIT VARIABLES  =================
    #==============================================================================

    # ----------- Short-circuit Intensity  -----------    
    # Eq. (3) of Lindsay et al. (2019)    
    Isc_theo = tools.Isc(sw_cell_hres, panel.SR, panel.alpha_SR, tc, panel.tcell_STC, panel.Ct_Isc, resolution = namelist.resolution)   

    # ------------ Open circuit voltage  -----------------  
    # Eq. (6) of Lindsay et al. (2019)
    Vt = tools.thermal_v(tc)
    Voc_theo = tools.Voc(Vt, panel.Voc_STC, tc, panel.tcell_STC, method = namelist.method_Voc, Isc = Isc_theo, Isc_STC = panel.Isc_STC, Ct_Voc = panel.Ct_Voc, lambda_gap = panel.lambda_gap)
    
    # ------------ Fill factore -------------      
    # FF = Pmpp/(Voc*Isc)    
    rs = panel.Rs/(Voc_theo/Isc_theo)
    Voc_norm = Voc_theo/Vt
    FF0 = tools.FF0(Voc_norm)
    FF = tools.FF(FF0, rs, method = namelist.method_FF)
       
    # ----------------- Power  ----------------   
    P = tools.power_module(FF, Isc_theo, Voc_theo, panel.A, panel.Nmodules)
       
    return P
