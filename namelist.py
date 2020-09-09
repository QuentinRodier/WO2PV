# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:01:17 2018

@author: lindsayn
"""

import numpy as np

# read the RRTM-G SW bands from the corresponding text file. Note that this is not monotonically icnreasing
def lambda_sw(datadir):
    data = np.loadtxt("%s/ecrad_bands_SW.txt"%datadir)
    lambda_max = data[0,::-1]
    lambda_min = data[1,::-1]
    lambda_mid = data[2,::-1]
    lambda_width = data[3,::-1]
    
    return lambda_max, lambda_min, lambda_mid, lambda_width

# Return ecRad bands for increasing wavelengths
def bands_ecrad(datadir):
    temp = lambda_sw(datadir)
    return np.append(temp[1][0],temp[0])      

class Namelist():  
    """Namelist object only has attributes that can be accessed from other modules"""
    def __init__(self):
        self.datadir = '/home/liboisq/Documents/Stagiaires/Nicole_Lindsay/PVCode_clean/Tools/data'
        self.swbands = [200., 4000.] # wavelength limits of the atmospheric swdir and swdiff inputs (nm) ; self.swbands = bands_ecrad(self.datadir) if usng ecrad bands
#        self.swbands = bands_ecrad(self.datadir)
        
        self.method_gap = 'default' # 'default", otherwise should be provided by the user
              
        # transposition method for diffuse irradiance: 'default', 'isotropic', 'isotropic', 'Klucher','HHay-Davies', 'Reindl', 'Perez', 'ISOTROPIC3D'
        # default uses different models depending on the SZA value
        # for reflected radiation, isotropic3D is used if prescribed, otherwise isotropic
        self.method_POA = 'default'
                
        # optical losses at the panel surface
        self.method_OL = 'reflection' # 'reflection_absorption', 'reflection', 'Martin'
        
        # Martin and Ruiz (2001) parameters
        self.ar = 0.157                 # see Table 1 of Martin (2001)
        self.c1 = 4./(3*np.pi)             # see Table 3 of Martin (2001)
        self.c2 = -0.074   
        
        # fresnel parameters for 'Reflection_Absorption', 'Reflection'        
        self.n0 = 1.              # air refractive index       
        self.n1 = 1.526           # module cover refractive index (generally glass)       
        self.K = 4.               # module cover absorption coefficient (m-1)        
        self.L = 2.*1e-3          # module cover thickness (m)
        
        # Cell temperature model can be set by the user, computed from module temperature set by the user, or computed from atmospheric parameters
        self.method_tc = 'King2004' # 'measured', 'measuredKing2004','King2004', 'Skoplaki2008_1', 'Skoplaki2008_2', 'Skoplaki2008_3', 'Skoplaki2008_4', 'simple'
        self.mount_type_King = 5    # only if King2004 is used, integer between 0 and 5, 5 is set by default

        # Skoplaki parameters
        self.eta_STC = 0.12
        self.tau_alpha = 0.9
        
        # Electrical model
        self.method_Voc = 'I0'                # 'I0' or 'STC'
        self.method_FF = 'complex'            # 'simple', 'complex'
        self.resolution = "artdeco"           # 'artdeco', 'astm'
        self.artdeco_resolution = "low"       # 'low', 'high', spectral resolution of reference ARTDECO spectra
        
        self.lambda_inf = 100.                # lower limit of spec tral response; change if SR to be cut before reaching zero
        self.SR_method = "user"               # ideal, real, user
    
        
        self.nthread = 10                      # time seris can be split into multi processes to distribute on various processors
        self.maxlen = 10000                   # only used if self.nthread = 1 to split large files (not necessarily useful, depends on memory resources), can be removed
        self.albedo = "weather"               # 'weather' -> albedo read from Atmosphere, otherwise should be the name of a surface defined in Tools/data/albedo (e.g. 'black', 'sand' etc.)
        
        # PANEL
        self.panel_name = "FranceWatts"       # 'FranceWatts', 'Sharp', 'Panasonic', 'FirstSolar', 'SolarFrontier';
                                              # can be appended in Tools/data/modules.dat; "measured test files correspond to FranceWatts"
                                              # alternatively, cell technology can be given : 'c-Si', 'a-Si/c-Si', 'HIT', 'CdTe', 'CIS'
        self.panel_technology = "fixed"       # fixed, tracker
        self.panel_beta = 27./180*np.pi       # inclination in degrees
        self.panel_gamma = np.pi              # orientation, [0-pi] is east, [pi-2pi] is west, pi is South