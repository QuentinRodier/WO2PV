# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 10:00:36 2017

@author: lindsayn
"""

import numpy as np
# import matplotlib.pyplot as plt
import netCDF4 as nc
import ephem
import scipy
import datetime
import sys
sys.path.append("..")
from namelist import Namelist

namelist = Namelist()
datadir = namelist.datadir

#==============================================================================
#=========================   CONSTANTS     ====================================
#==============================================================================

#Planck constant
h_eV = 4.135667662e-15    # eV⋅s
h = 6.626070040e-34       # J⋅s

#Speed of light
c = 299792458             #m s-1

#elemental charge
q = 1.60217662e-19        #C

#Boltzmann constant
kB=1.38064852e-23         #J⋅K−1
kB_eV=8.6173303e-5        #eV⋅K−1

# Solar constant
S0 = 1366.                #  W m-2
#==============================================================================

#==============================================================================
#========================    SUN GEOMETRY    ==================================
#==============================================================================
def get_sun(lat,lon,date):
    obs = ephem.Observer()
    obs.lat = str(lat)
    obs.long = str(lon)
    obs.date = date
    sun = ephem.Sun(obs)
    sun.compute(obs)
    sza = np.pi/2 - sun.alt
    saa = float(sun.az)
    distance = sun.earth_distance
    sw_toa = S0*(1./distance)**2
       
    return sza, saa, sw_toa

#==============================================================================
#=================  ATMOSPHERE EXTRACTION FROM NETCDF  ========================
#==============================================================================
def get_atm_from_ncfile(forcing, albedo, n1 = 0, n2 = -1):
    data =  nc.Dataset("Forcings/%s.nc"%forcing, "r") 
    lat = data.variables["lat"][:]
    lon = data.variables["lon"][:]
    time = data.variables["time"][:]
    n2 = min(n2,max(n2, len(time)))
    dates = [datetime.datetime(2000,1,1,0,0,0) + datetime.timedelta(seconds=int(t)) for t in time[n1:n2]]
    ta = data.variables["ta"][:,n1:n2]
    wind = data.variables["wind"][:,n1:n2]
    swdir = data.variables["swdir"][:,n1:n2,:]
    swdiff = data.variables["swdiff"][:,n1:n2,:]
    

    if albedo == "weather":
        albdir = data.variables["albdir"][:,n1:n2,:]
        albdiff = data.variables["albdiff"][:,n1:n2,:]
             
    else:
        ncol, nt = np.shape(ta)
        data_alb = np.loadtxt("%s/albedo/%s.dat"%(datadir,albedo))
        wvl_albedo = data_alb[:,0]
        albedo = data_alb[:,1]
        albedo = np.interp(artdeco_wls, wvl_albedo, albedo)
        nwvl = len(albedo)
        albdir = np.zeros([ncol, nt, nwvl])
        albdiff = np.zeros([ncol, nt, nwvl])
        albdir[:,:,:] = albedo
        albdiff[:,:,:] = albedo
        
    return lat, lon, dates, swdir, swdiff, ta, wind, albdir, albdiff
    

#==============================================================================
#========================    REF IRRADIANCE  ==================================
#==============================================================================

# --------------- ASTM spectrum --------------------
def ASTM_spectrum():
    data_ASTM = np.loadtxt('%s/ASTM_spectrum.txt'%datadir)
    wvl_ASTM = data_ASTM[:,0]  # nm
    sw_ASTM = data_ASTM[:,1] # W m-2 nm-1
   
    return wvl_ASTM, sw_ASTM
    
wvl_ASTM, sw_ASTM = ASTM_spectrum()

# -----Create reference irradiance file from files computed with ARTDECO-------
def compute_reference_spectra():
    if namelist.artdeco_resolution == "high":
        artdeco_file = "%s/ARTDECO_Spectra_surf.nc"%datadir
    else: 
        artdeco_file = "%s/ARTDECO_Spectra_surf_lowres.nc"%datadir
        
    data = nc.Dataset(artdeco_file, "r")
    wls = 1e3*data.variables["Wv"][:]  # microns -> nanometers
    dwls = np.diff(wls)                # bandwidth
    artdeco_wls = np.asarray(0.5*(wls[:-1]+wls[1:]))
    flux_down = data.variables["flux_dw"][:,:] # wl, sza
    flux_dir = data.variables["flux_dir"][:,:] # W m-2 per spectral interval
    flux_diff = flux_down - flux_dir
    flux_diff = flux_diff[1:,:]/dwls[:, None] # W m-2 nm-1
    flux_dir = flux_dir[1:,:]/dwls[:, None]
    
    artdeco_sza = data.variables["Sza"][:]
    # [ 0.  5. 10. 15. 20. 25. 30. 35. 40. 45. 50. 55. 60. 65. 70. 75. 80. 85.]
    nsza = len(artdeco_sza)
    nwl = len(artdeco_wls)                            
    artdeco_spec = np.zeros([nsza,nwl,2])
    for i in range(nsza):
        artdeco_spec[i,:,0] = flux_dir[:,i]
        artdeco_spec[i,:,1] = flux_diff[:,i]    
                
    return artdeco_wls, artdeco_spec, artdeco_sza

artdeco_wls, artdeco_spec, artdeco_sza = compute_reference_spectra()

# read the ASTM  wavelengths    
def wvl_astm():
    return wvl_ASTM

# read the MODTRAN wavelengths   
def wvl_artdeco():
    return artdeco_wls

def compute_ASTM_spectrum():
    """Simulated ASTM direct and diffuse spectra with ARTDECO"""
    if namelist.artdeco_resolution == "high":
        artdeco_file = "%s/ARTDECO_ASTM_Spectrum_surf.nc"%datadir
    else:
        artdeco_file = "%s/ARTDECO_ASTM_Spectrum_surf_lowres.nc"%datadir   
        
    data = nc.Dataset(artdeco_file, "r")
    wls = 1e3*data.variables["Wv"][:]  # microns -> nanometers
    dwls = np.diff(wls) # bandwidth
    artdeco_wls = np.asarray(0.5*(wls[:-1]+wls[1:]))
    
    flux_down = data.variables["flux_dw"][:,0] # wl, sza, surface
    flux_dir = data.variables["flux_dir"][:,0] # W m-2 per spectral interval
    flux_diff = flux_down - flux_dir
    artdeco_diff = flux_diff[1:]/dwls # W m-2 nm-1
    artdeco_dir = flux_dir[1:]/dwls
                            
    return artdeco_wls, artdeco_dir, artdeco_diff
    

# ------ Increasing spectral resolution ------------------  
def hres(sw_in, bands = [0,4000], sza=0., method='single-band', source = "direct"):   
    """Take a low resolution spectrum over a few bands and
    return a spectrum with resolution of ASTM or ARTDECO"""
    
    # artdeco_sza [ 0.  5. 10. 15. 20. 25. 30. 35. 40. 45. 50. 55. 60. 65. 70. 75. 80. 85.]
    # Selection of high resolution reference spectrum
    if source == "direct":
        index = 0
        
    elif source == "diffuse":
        index = 1
  
    sw = np.zeros((np.shape(sw_in)[0], len(artdeco_wls))) 
     
    inf_sza = np.where(sza<=artdeco_sza[0])[0]
    sw[inf_sza,:] = artdeco_spec[0,:,index]
    
    sup_sza = np.where(sza>=artdeco_sza[-1])[0]
    sw[sup_sza,:] = artdeco_spec[-1,:,index]
           
    # interpolation on spectra
    in_sza = np.where((sza<artdeco_sza[-1]) & (sza>artdeco_sza[0]))[0]
    indx_temp = np.searchsorted(artdeco_sza,sza)
    f1 = artdeco_spec[indx_temp,:,index]  # interpolation between successive spectra
    f2 = artdeco_spec[indx_temp+1,:,index]
    sw[in_sza,:] = (sza[in_sza,None] - artdeco_sza[indx_temp][in_sza,None])/(artdeco_sza[indx_temp+1][in_sza,None] - artdeco_sza[indx_temp][in_sza,None])*(f2[in_sza,:]-f1[in_sza,:])+f1[in_sza,:]      

    # Converting from low to high resolution
    sw_hres = np.zeros((np.shape(sw_in)[0], len(artdeco_wls))) 
    
    for i in range(len(bands)-1):
        i1, i2 = np.searchsorted(artdeco_wls, (bands[i],bands[i+1]))
        sw_band = np.zeros(np.shape(sw_in)[0])

        if i1 != i2:         
            sw_band[:] = scipy.integrate.simps(sw[:,i1:i2], artdeco_wls[i1:i2])  # energy in the band
    
        ze = np.where(sw_band != 0)[0]
        if method == 'multi-band':
            sw_hres[ze,i1:i2] = sw_in[ze,i,None]/sw_band[ze,None]*sw[ze,i1:i2]
            
        if method == 'single-band': # one single band
            sw_hres[ze,i1:i2] = sw_in[ze]/sw_band*sw[ze,i1:i2]
            
    return sw_hres   


# -----Get ordered spectral fluxes from ecRad output -----
def get_spectral_flux_ecrad(filename, nswband):
    spec_liste = ["spectral_flux_dn_direct_sw_surf", "spectral_flux_dn_sw_surf"]
    dim_spec=len(spec_liste)
    spec_flux = np.zeros([nswband,2])
    data = nc.Dataset(filename,'r')
    
    for k,spec in enumerate(spec_liste):
        data0 = data.variables[spec][0,:]
        spec_flux[:,k] = np.append(output0[-1],output0[0:-1]) # switch last band
        spec_flux[:,k] = spec_flux[:,k][::-1]                 # increasing wavelength
        
    return spec_flux # direct and total fluxes

#==============================================================================
#==================     TRANSPOSITION AND OPTICAL LOSSES     ==================
#==============================================================================

# ------ Incident angle on an inclined surface  ---------
def incident_angle_panel(beta,gamma,sza,saa):
    # Duffie and Beckman [2013], eqn (1.6.3) p.14. (Solar Engineering of Thermal Processes)
    # Should all be radians
    cos_theta = np.cos(sza)*np.cos(beta)+np.sin(sza)*np.sin(beta)*np.cos(saa-gamma) # saa and gamma should have same reference (North = 0)
    theta = np.arccos(cos_theta)
    
    return theta

# Effective incident angles to compute optical losses
# Eq .(6) of Brandemuehle and Beckman (1980)  
def effective_ground_reflected_incident_angle(beta):
    return (90 - 0.5788*beta/np.pi*180. + 0.002693*(beta/np.pi*180.)**2)/180*np.pi  

#effective beam radiation incident angle for ground_reflected radiation
# Eq .(7) of Brandemuehle and Beckman (1980)
def effective_diffuse_incident_angle(beta):
    return (59.68-  0.1388*beta/np.pi*180. + 0.001497*(beta/np.pi*180.)**2)/180*np.pi        
    
# Transposition from horizontal to tilted plane - return coefficient to multiply input irradiance by
#==============================================================================

def Erbs_model(Kt): # Eq. (1) of Erbs et al. (1982)
	if Kt<=0.22:
		fd = 1-0.09*Kt
	elif (0.22<Kt)*(Kt<=0.80):
		fd = 0.9511-0.1604*Kt+4.388*Kt**2-16.638*Kt**3+12.338*Kt**4
	else:
		fd= 0.165 # clear sky
	return fd

def anisotropy_index(DNI,sw_toa):
    return DNI/sw_toa

# ---- Air mass ------------
# Kasten & Young 1989
def air_mass(sza):
    airmass = 1/(np.cos(sza)+0.50572*(96.07995-sza/np.pi*180)**(-1.6364))
    if np.isscalar(sza) :
        if (sza<0 and sza>(np.pi/2*1.01)):
            airmass = np.nan
    else :
        airmass[sza<0] = np.nan
        airmass[sza>(np.pi/2*1.01)] = np.nan
        
    return airmass

def transposition_direct(theta,sza):
    # Fdir/cos(sza)*cos(theta) ->  1/cos(sza)*cos(theta) if theta in [-pi/2; pi/2]; 0 otherwise
    # ((cos(theta))>0) is to ensure that the beam radiation is not incident behind the panel
    return (np.cos(theta)/np.cos(sza))*((np.cos(theta))>0)

def transposition_diffuse(beta, method='ISOTROPIC3D', theta = None, sza = None, swdir_bb = None, swdiff_bb = None, sw_bb = None, sw_toa = 1361):
    # POA depends on beta and sza
    if method == 'isotropic': # Eq. (4) of Badescu (2002)
        # Solid angle seen from panel
        return (1. + np.cos(beta))/2
    
    elif method == 'isotropic3D': # Eq. (11) of Badescu 2002
        return (3. + np.cos(2*beta))/4
        
    elif method == 'Klucher': # Eq. (3) of Klucher 1979
        F = 1 - (swdir_bb/sw_bb)**2
        return (1. + np.cos(beta))/2*(1.+F*np.sin(beta/2)**3)*(1.+F*np.cos(theta)**2*sin(sza)**3)  
        
    elif method=='Hay-Davies': # Eq. (7) of Loutzenhiser et al. (2007) (Hay and Davies 1980)
        A = anisotropy_index(swdir_bb/np.cos(sza),sw_toa)   # DNI computation      
        return transposition_direct(theta, sza)*A+(1-A)*(1.+np.cos(beta))/2
        
    elif method=='Reindl': # # Eq. (7) of Loutzenhiser et al. (2007) (Reindl 1990)
        A = anisotropy_index(swdir_bb/np.cos(sza),sw_toa)
        return transposition_direct(theta, sza)*A+(1-A)*(1.+np.cos(beta))/2*(1.+np.sqrt(swdiff_bb/sw_bb)*np.sin(beta/2)**3)
        
    elif method=='Perez': # Reindl 1990
        sza = np.minimum(sza, np.pi/2-1e-12)
        
        a = np.maximum(0., np.cos(theta)) # a=max(0.,cos(theta))
        b = np.maximum(0.087, np.cos(sza))

        epsilon = ((swdiff_bb + swdir_bb/np.cos(sza))/swdiff_bb +1.041*sza**3)/(1+1.041*sza**3) # Eq.(1) of Perez (1990)
                
        delta = air_mass(sza)*swdiff_bb/(sw_toa) # Eq.(2) of Perez (1990) warning: simplification - relative approx absolute air mass (i.e. no altitude correction) - not clear if I0 is S0 or S0cos(SZA) but the latter does not work
                
        # Clearness based on epsilon value
        categories_epsilon = np.loadtxt('%s/Perez_clearness_index.txt'%datadir)   # Table 1 of Perez (1990)
        categories_epsilon[-1,-1] = np.inf # to be sure that any large value falls in the last bin
        
        if type(epsilon) is float:
            index_cat = 0  # to have the same dimensions as epsilon
            
        else:
            index_cat = np.zeros_like(epsilon) # to have the same dimensions as epsilon
            
        for i in range(len(categories_epsilon)):
            # To append a vector with index_cat for each epsilon of the input (selected for epsilon in range)
            index_cat+= i*((categories_epsilon[i,1] <= epsilon)*(epsilon < categories_epsilon[i,2]))
       
        index_cat = index_cat.astype(int)       
        coeff_file = np.transpose(np.loadtxt('%s/Perez_LUT.txt'%datadir))   # Table 6 of Perez (1990)
        f = coeff_file[1:,index_cat] # each colum of f is the right combination of parameters to be used
        
        F1 = np.maximum(0, f[0]+f[1]*delta + sza*f[2]) # Circumsolar Brightening Coefficient (f[0] is the first line of f)
        F2 = f[3] + f[4]*delta + sza*f[5] #Horizon Brightening Coefficient  
                
        return (1-F1)*(1+np.cos(beta))/2 + F1*a/b + F2*np.sin(beta) # Eq. (9) of Perez (1990)
            
      
def transposition_reflected(beta, method = 'isotropic'):
    if method == 'isotropic': 
        return (1. - np.cos(beta))/2    # Eq. (7) of Badescu 2002
        
    elif method == 'isotropic3D': 
        return (1. - np.cos(2*beta))/4  # Eq. (14) of Badescu 2002

def transposition(theta, sza, beta, swdir_bb = None, swdiff_bb = None, sw_bb = None, sw_toa = None):
    # default transposition
    # return Tb, Td, Tr with distinction according to SZA, based on observations at SIRTA
    # different choice (isotropic3D or Perez) for diffuse method (SZA>70 or < 70 deg)
    method_temp = (sza>70./180*np.pi)

    # 0 if SZA<=70, 1 if SZA>70 -> 0 = ISTROPIC 3D, 1= PEREZ (diffuse) & ISOTROPIC (reflected)
    TD_ISO3D = transposition_diffuse(beta, method = 'isotropic3D')
    TD_PEREZ = transposition_diffuse(beta, method='Perez', theta = theta, sza = sza, swdir_bb = swdir_bb, swdiff_bb = swdiff_bb, sw_bb = sw_bb, sw_toa = sw_toa)
#    input()
    TR_3D = transposition_reflected(beta, method='isotropic3D')
    TR_2D = transposition_reflected(beta, method='isotropic')
    
    return transposition_direct(theta,sza), method_temp*TD_PEREZ + (1-method_temp)*TD_ISO3D, method_temp*TR_2D + (1-method_temp)*TR_3D   


#------------ POA IRRADIANCE  ------------
def sw_poa(swdir_hres, swdiff_hres, swref_hres, transpose_dir, transpose_diff, transpose_ref):
    # coefficients have dimension of times, need to add spectral dimension for multiplication
    if not np.isscalar(transpose_dir):
        transpose_dir = transpose_dir[:,None]
   
    if not np.isscalar(transpose_diff):
        transpose_diff = transpose_diff[:,None]
           
    if not np.isscalar(transpose_ref):
        transpose_ref = transpose_ref[:,None]

    return swdir_hres*transpose_dir, swdiff_hres*transpose_diff, swref_hres*transpose_ref


# -----  Optical losses ----------
    
# Direct irradiance, Eq. (6a) of Martin (2001)
def ft_b(theta, ar):
    return (np.exp(-np.cos(theta)/ar)-np.exp(-1/ar))/(1-np.exp(-1/ar)) 

# Diffuse irradiance, Eq. (6c) of Martin (2001)
def ft_d(beta, ar, c1, c2):
    return np.exp(-1/ar*(c1*(np.sin(beta)+(np.pi-beta-np.sin(beta))/(1 + np.cos(beta))) \
    +c2*(np.sin(beta) + (np.pi - beta - np.sin(beta))/(1 + np.cos(beta)))**2))

# Reflected irradiance, Eq. (6b) of Martin (2001)
def ft_r(beta, ar, c1, c2):
    output= np.exp(-1/ar*(c1*(np.sin(beta) + (beta-np.sin(beta))/(1 - np.cos(beta))) \
    +c2*(np.sin(beta) + (beta - np.sin(beta))/(1 - np.cos(beta)))**2))
    
    indx = np.where(beta==0)[0]
    if indx!= []:
        output[indx] = np.ones(len(indx)) # limit for beta->0 of eq. 6b
    return output


#------ Based on Fresnel  ----------

# Simple air-glass model,  reflection and absorption. De Soto 2004
# Model for 1 layer
# Note : can also be applied for other thickness, etc by changing default parameters.
# Transmittance = (1-Reflection)*(1-Absorption)
def reflection_absorption(theta, K=4., L=2.*1e-3, n0=1, n1=1.526, absorption = False): # De Soto 2006,IN TERMS OF TRANSMITTED
    
    theta = np.minimum(theta, np.pi/2)
    
    theta_r = np.arcsin(n0/n1*np.sin(theta)) # angle of refraction, Snell's law
    
    if isinstance(theta,(float)):
        if theta == 0.:
            temp = 1-(1.-n0/n1)**2/(1.+n0/n1)**2 # Taylor series when theta->0
        else:
            temp = (1-0.5*((np.tan(theta-theta_r)**2)/(np.tan(theta+theta_r)**2))-0.5*((np.sin(theta - theta_r)**2)/(np.sin(theta + theta_r)**2)))    # Eq. (2.8) of DeSoto (2004)          
    else:
        temp = ( 1-0.5*((np.tan(theta - theta_r)**2)/(np.tan(theta+theta_r)**2))-0.5*((np.sin(theta - theta_r)**2)/(np.sin(theta + theta_r)**2)))
        temp[theta == 0.] = (1-(1.-n0/n1)**2/(1.+n0/n1)**2)*np.ones(len(temp[theta==0.]))
        
    return temp * (1-absorption) + absorption*np.exp(-K*L/np.cos(theta_r))*temp
    


# --------  Wraper for transmission coefficient --------------
     
def optical_losses(theta, beta, method = 'reflection', ar = 0.157, c1 = 4./(3*np.pi), c2 = -0.074, n0 = 1., n1 = 1.526, K = 4., L = 2.*1e-3):  
    if method == 'Martin': #-N Martin and JM Ruiz [2001]
        return 1 - ft_b(theta,ar), 1-ft_d(beta,ar,c1,c2), 1-ft_r(beta,ar,c1,c2) # Eq. (5) of Martin (2001)
        
    elif method == 'reflection_absorption': # De Soto [2004] # ol_diff and ol_ref are scalars
        return reflection_absorption(theta, K, L, n0, n1, absorption = True), reflection_absorption(effective_diffuse_incident_angle(beta),K,L,n0,n1, absorption = True), reflection_absorption(effective_ground_reflected_incident_angle(beta),K,L,n0,n1, absorption = True)
        
    elif method == 'reflection': # Sjerps [1995]
        return reflection_absorption(theta, K, L, n0, n1, absorption = False), reflection_absorption(effective_diffuse_incident_angle(beta),K,L,n0,n1, absorption = False), reflection_absorption(effective_ground_reflected_incident_angle(beta),K,L,n0,n1, absorption = False)

        
#-------- Incoming flux on cell ----------- (POA * transmission)        

def sw_cell(swdir_poa_hres, swdiff_poa_hres, swref_poa_hres, ol_dir, ol_diff, ol_ref):
    if not np.isscalar(ol_dir):
        ol_dir = ol_dir[:,None]
    
    if not np.isscalar(ol_diff):
        ol_diff = ol_diff[:,None]      
    
    if not np.isscalar(ol_ref):
        ol_ref = ol_ref[:,None]  
    
    return swdir_poa_hres*ol_dir, swdiff_poa_hres*ol_diff, swref_poa_hres*ol_ref    

# -------- Spectral flux impacting the cell in STC ----------------
def sw_cell_stc(method_OL = 'reflection',ar = 0.157, c1 = 4./(3*np.pi), c2 = -0.074, n0 = 1., n1 = 1.526, K = 4., L = 2.*1e-3):
    # Goal is to get the 3 components of the ASTM POA spectrum from ARTEDCO simulations to compute propertly optical losses
    # STC illumination conditions
    sza_stc = 48.19/180*np.pi
    beta_stc = 37./180*np.pi
    sw_albedo_stc = 0.2
#    sw_bb_stc = 1000.
    
    theta_stc = incident_angle_panel(beta_stc, 0, sza_stc, 0)
    
    wvl_ASTM, sw_ASTM = ASTM_spectrum()
    wvl_artdeco, swdir_artdeco, swdiff_artdeco = compute_ASTM_spectrum()
    
    swref_artdeco = sw_albedo_stc*(swdir_artdeco + swdiff_artdeco)
    swdir_artdeco_bb = scipy.integrate.simps(swdir_artdeco, wvl_artdeco)
    swdiff_artdeco_bb = scipy.integrate.simps(swdiff_artdeco, wvl_artdeco)
    
    Tdir_stc, Tdiff_stc, Tref_stc = transposition(theta_stc, sza_stc, beta_stc, swdir_bb = swdir_artdeco_bb, swdiff_bb = swdiff_artdeco_bb, sw_bb = swdir_artdeco_bb + swdiff_artdeco_bb, sw_toa = 1361.)
    
    swdir_poa, swdiff_poa, swref_poa = sw_poa(swdir_artdeco, swdiff_artdeco, swref_artdeco, Tdir_stc, Tdiff_stc, Tref_stc)
    sw_poa_hres = swdir_poa + swdiff_poa + swref_poa
    i1, i2 = np.searchsorted(wvl_artdeco, (wvl_ASTM[0], wvl_ASTM[-1]))

    sw_poa_bb = scipy.integrate.simps(sw_poa_hres[i1:i2], wvl_artdeco[i1:i2])
        
#    print("SW POA ARTDECO", sw_poa_bb, "(should be 1000 W m-2)")
    
    if method_OL == 'None':
        ol_stc = 1.   
        
    else: 
        ol_dir, ol_diff, ol_ref = optical_losses(theta_stc, beta_stc, method = method_OL, ar = ar, c1 = c1, c2 = c2, n0 = n0, n1 = n1, K = K, L = L)
        swdir_cell, swdiff_cell, swref_cell = sw_cell(swdir_poa, swdiff_poa, swref_poa, ol_dir, ol_diff, ol_ref)
        sw_cell_hres = swdir_cell + swdiff_cell + swref_cell
        swcell_bb = scipy.integrate.simps(sw_cell_hres[i1:i2], wvl_artdeco[i1:i2])
        
        ol_stc = swcell_bb/sw_poa_bb # average losses

    return ol_stc*sw_ASTM  # the spectrum is provided at ASTM resolution


#==============================================================================
#========================    PV CALCULATIONS    ===============================
#==============================================================================

# -------------- Energy gap -------------
# Lambda Gap in nm from the band gap energy Eg in eV
def energy_gap_j(lambda_g):
    # Eq. (1) of Lindsay (2019)
    return (h*c)/(lambda_g*1e-9)
 
def lambda_gap_eV2nm(Egap):
    return (h_eV*c)/Egap*1e9
 
def lambda_gap(celltype='default', Egap=0):
    if celltype=='default':
        lgap = lambda_gap_eV2nm(Egap)
        
    else:
        data = np.loadtxt('%s/lambda_gap.txt'%datadir, dtype=np.str)
        lgap = 'not in list'
        for i in range(len(data)):
            if data[i][0] == celltype:
                lgap = float(data[i][2])
                
    return lgap


# ------------ Spectral response ------------

#l_mid = midwidth spectral length for each band (array)
#l_lim = limit under which not absorbed (cf glass : approx 400nm) (scalar)
#l_gap = gap defined by the type of cell (scalar)
    
def qe_theo(wvl, wvl_inf, wvl_sup):
    temp = (wvl >= wvl_inf)*(wvl<wvl_sup) # ideal efficiency (0 or 1)
    return temp.astype(np.int)    

def spectral_response_theo(wvl,QE):
    # Eq. (2) of Lindsay (2019)
    return q/(h*c)*10**-9*QE*wvl

# Reading empirical spectral response from typical file
def spectral_response_empirical(celltype):
    filedir = '%s/SR_empirical'%datadir
    list_panels = ['CdTe', 'CIGS', 'c-Si', 'DSC', 'GaAs', 'mc-Si', 'OSC']
    
    if celltype in list_panels:
        indx = list_panels.index(celltype)
        fname='%s/SR_%s.csv'%(filedir, list_panels[indx])
        data = np.loadtxt(fname)
        wls = data[0]
        SR = data[1]
        
    else:
        SR = 'not in files'
    return wls, SR


def spectral_response(celltype = 'theoretical', wvl_inf = 400, wvl_sup = 1100, resolution = 'astm', wvl = wvl_ASTM):
    """Return the spectral response at the required spectral resolution (ASTM or ARTDECO)
    Either ideal spectral response or empirical"""
    
    if celltype == 'theoretical':
        if resolution == 'artdeco':
            wvl = wvl_artdeco()
            
        SR = spectral_response_theo(wvl, qe_theo(wvl, wvl_inf, wvl_sup))
        
    else:
        wls, SR = spectral_response_empirical(celltype)
        SR = np.interp(wvl, wls, SR)
        
    return SR
    
def correction_sr(SR, Isc_ref, spectrum_ref = sw_ASTM, wvl_ref = wvl_ASTM):  
    # Eq. (10) of Lindsay (2019)
    return Isc_ref/scipy.integrate.simps(SR*spectrum_ref, wvl_ref) 


# ------------ Cell temperature -------------

def mod2cell_temperature(Tmod,POA,mount_type=5):
    King_table = np.loadtxt('%s/King2004_cell_temperature.txt'%datadir,skiprows=1,usecols=[2,3,4])
    return Tmod+POA/1000*King_table[mount_type,2]      # Eq. (12) of King (2004)

def cell_temperature(ta, sw_poa_bb, method='SIMPLE', NOCT = 45+273.15, ta_NOCT = 20+273.15, sw_NOCT = 800., CtPmax = None, eta_STC = None, wind = None, mount_type = 5, tau_alpha = 0.9, ta_STC = 25+273.15):
    if method == 'Simple': # "NOCT-Standard-formula"
        return ta + sw_poa_bb/sw_NOCT*(NOCT - ta_NOCT)
    
    elif method == 'King2004': 
        data_king = np.loadtxt('%s/King2004_cell_temperature.txt'%datadir, skiprows=1, usecols = [2,3,4])
        Tmod = sw_poa_bb*np.exp(data_king[mount_type,0] + data_king[mount_type,1]*wind) + ta # Eq. (11) of King (2004)
        return Tmod + sw_poa_bb/1000*data_king[mount_type, 2] # Eq. (12) of King (2004)
    
    elif method == 'Skoplaki2008_1': # Skoplaki [2008]
        return ta + sw_poa_bb/sw_NOCT*(NOCT-ta_NOCT)*(8.91+2.0*1)/(8.91+2.0*wind)*(1-(eta_STC/tau_alpha)*(1-CtPmax/100*(ta_STC-273.15))) # Eq. (16) and  (20)
    
    elif method == 'Skoplaki2008_2': # Skoplaki [2008]
        return ta + sw_poa_bb/sw_NOCT*(NOCT - ta_NOCT)*(5.7+2.8*1)/(5.7+2.8*wind)*(1-(eta_STC/tau_alpha)*(1-CtPmax/100*(ta_STC-273.15))) # Eq. (16) and  (22)
    
    elif method == 'Skoplaki2008_3': # Skoplaki [2008]
        return ta + 0.25/(5.7+3.8*wind)*sw_poa_bb # Eq. (23) of Skoplaki (2008)
    
    elif method == 'Skoplaki2008_4': # Skoplaki [2008]
        return ta + sw_poa_bb/sw_NOCT*(NOCT - ta_NOCT)*(5.7+2.8*1)/(5.7+2.8*wind)*(1-(eta_STC/tau_alpha)*(1-CtPmax/100*(ta_STC-273.15)+0.12*np.log(sw_poa_bb))) # Eq. (16) and  (22) and (4)
    
    else:
        raise("No such method for cell temperature")

     
#==============================================================================
# Electrical variables / output calculation
#==============================================================================

def series_resistance(FF, FF0, Voc, Isc, method='complex'):
    if method == 'simple':    
        Rs_temp=(1-FF/FF0)*Voc/Isc   # Eq. (6) of Green (1982)
        
    elif method == 'complex':
        # print("Rs possible values, 2nd chosen",np.roots([1/5.4,-1.1*FF0,FF0-FF])[1])
        Rs_temp = np.roots([1/5.4,-1.1*FF0,FF0-FF])[1]*(Voc/Isc)    # Eq. (7) of Green (1982)
        
    return  Rs_temp

def thermal_v(T, n=1):
    # n = cell ideality factor, set to 1 here
    return n*kB*T/q

def I0(Isc_STC, Voc_STC, lambda_g, tc, tc_STC):
    Vt_STC_temp = thermal_v(tc_STC)
    Vt_temp = thermal_v(tc)
    
    Egap = energy_gap_j(lambda_g)
    return Isc_STC*np.exp(-Voc_STC/Vt_STC_temp)*(tc/tc_STC)**3*np.exp(Egap/q*(1./Vt_STC_temp-1./Vt_temp)) # Eq. (6) of Lindsay (2019)

def FF0(Voc_norm):
    # Eq. (5) of Green (1982) 
    return (Voc_norm - np.log(Voc_norm+0.72))/(Voc_norm+1)

# Green [1982]      
def FF(FF0,rs,method='simple'):
    if method == 'simple':
        temp=FF0*(1-rs)
    elif method== 'complex':
        temp=FF0*(1-1.1*rs)+rs**2/5.4
    return temp

# Green [1981]     
def rs(Rs, Voc, Isc):
    return Rs/(Voc/Isc)

# Lorenzo [2003]
def Isc(sw_cell_hres, SR, alpha_SR, tc, tc_STC, Ct_Isc, resolution = 'artdeco'):
    if resolution == 'artdeco':
        wvl = wvl_artdeco()
        
    elif resolution =='astm':
        wvl = wvl_astm()
    
    else:
        print("Unkonwn resolution option")
        
    if sw_cell_hres.ndim == 1:
        Iref = scipy.integrate.simps(alpha_SR*SR*sw_cell_hres, wvl)
        
    else:      
        Iref = scipy.integrate.simps(alpha_SR*SR*sw_cell_hres, wvl[None,:])
               
    return Iref*(1 + (tc - tc_STC)*Ct_Isc/100) # Eq. (3) of Lindsay (2019)


# Lorenzo [2003]
def Voc(Vt, Voc_STC, tc, tc_STC, method='STC', Ct_Voc = None, Isc = None, Isc_STC = None, lambda_gap = None):
    if method == 'STC':
        # Go to Voc_STC and add temperature and irradiance contributions
        Voc_temp = (Voc_STC + thermal_v(tc_STC)*np.log(Isc/Isc_STC))*(1+(tc-tc_STC)*Ct_Voc/100) 
        
    if method == 'I0':
        # Eqs. (5) and (6) of Lindsay (2020)
        I = I0(Isc_STC, Voc_STC, lambda_gap, tc, tc_STC)
        Voc_temp = Vt*np.log(Isc/I)
        
    return Voc_temp       

def power(FF, Isc, Voc):
    return FF*Isc*Voc
    
def power_module(FF, Isc, Voc, A, M):
    # power is per cell and per m2 of cell so need to multiply by Ncell * Acell (ie A) to get Power
    return FF*Isc*A*M*Voc

def pv_all(Jeff,Tc,Jeff_STC,Jsc_STC,Voc_STC,Tc_STC,Ct_Jsc,Ct_Voc,Rs):
    # Return all PV outputs
    # to be used for detailed analysis
    out_Vt = thermal_v(Tc)
    out_Jsc = Jsc(Jsc_STC,Jeff,Jeff_STC,Tc,Tc_STC,Ct_Jsc)
    out_Voc = Voc(Voc_STC,Jeff, Jeff_STC,Tc,Tc_STC,Ct_Voc,out_Vt)
    out_rs = rs(Rs,out_Voc,out_Jsc)
    out_voc = voc(out_Voc,out_Vt)
    out_FF = FF(out_FF0,out_rs)
    out_P = POWER(out_FF, out_Isc, out_Voc)
    
    return np.asarray([out_P,out_Vt,out_Jsc,out_Voc,out_rs,out_voc,out_FF0,out_FF])



