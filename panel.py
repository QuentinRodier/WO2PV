# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:09:48 2018

@author: lindsayn
"""

import numpy as np
import tools

#==============================================================================
# INITIALISATION
#==============================================================================
class Panel():    
    def __init__(self, name = "FranceWatts", technology = "fixed", beta = 0, gamma = 0, lambda_gap = 1100., tcell = 0, method_FF = "simple", SR_method = "ideal", atm = None, Nmodules = 1., namelist = None):          
        """ Define all parameters required for the panel
        Technical characteristics are from datasheet modules.txt
        Geometry is provided by the user"""
        

        
        # Panel hardware
        datadir = namelist.datadir
        self.set_module(name, datadir)
        self.technology = technology                                                 # technology of the panel (fixed or tracker)       
        
        if self.technology == "fixed":  
            self.set_geometry(beta, gamma)                                           # done only once, otherwise done for each column in pv_interface
            
        
        self.Nmodules = Nmodules                                                     # number of modules
        print("The %s panel with %s technology was selected"%(self.name,self.type))
        
        # User defined characteristics
        self.tcell = tcell                                                           # cell temperature if not computed from air temperature                            
        
        # Standard test conditions
        self.ta_NOCT = 20 + 273.15                                                   # Air temperature at Nominal operating cell temperature  
        self.sw_NOCT = 800.                                                          # GHI at NOCT (W m-2)
        self.sw_STC = 1000.                                                          # Standard Test conditions POA (W m-2)
        self.ta_STC = 25. + 273.15                                                    # Standard Test conditions temperature (°C)
        
        # Computed PV quantities
        self.Vt_STC = tools.thermal_v(self.tcell_STC)                                # thermal voltage, Eq. (5) of Lindsay et al. (2019)
        self.set_rs(method_FF = method_FF)            
        
        # Spectral response
        if namelist.method_gap == 'default':                                         # cutoff wavelength of panel
            self.lambda_gap = tools.lambda_gap(celltype = self.type)                 # standard cut off wavelength for this technology
            
        else:                                                          
            # User defined or default value taken           
            print("Check that lambda_gap was properly set by the user")
            self.lambda_gap = lambda_gap                                             # provided by the user
                      
        self.set_SR(namelist) 
        
    def set_module(self, name, datadir):
        filename_modules = '%s/modules.txt'%datadir # contains a list of standard panels
        data_modules = np.genfromtxt(filename_modules,skip_header = 1,usecols = range(2,14))
        type_modules = list(np.genfromtxt(filename_modules,dtype=np.str,skip_header=1,usecols=[1]))
        name_modules = list(np.genfromtxt(filename_modules,dtype=np.str,skip_header=1,usecols=[0]))
         # Select panel according to technology or name of the panel
        try:
            index = name_modules.index(name)
            pass
        except:
            try:
                index = type_modules.index(name)    
                pass
            except:
                print("%s: No such module found"%name)
                raise                       
                 
        # Panel technical characteristics
        self.name = name_modules[index]                                              # name of the module
        self.type = type_modules[index]                                              # technology of the cells
        self.A = data_modules[index,0]                                               # total surface of the module (m2)
        self.Ncells = data_modules[index,1]                                          # number of cells per panel   
        self.Acell = self.A/self.Ncells                                              # area of one cell
        
        # Normalized quantities for the performances
        self.Pmax_STC = data_modules[index,2]/self.A                                 # max power per m2   
        self.Vmpp_STC = data_modules[index,3]/self.Ncells                            # voltage per cell at Pmpp (cells in series) is physically meaningful 
        self.Impp_STC = data_modules[index,4]/self.Acell                             # intensity at Pmpp for one cell, normalized by cell area, this ensures that P = U * I
        self.Voc_STC = data_modules[index,5]/self.Ncells                             # voltage open circuit per cell
        self.Isc_STC = data_modules[index,6]/self.Acell                              # intensity short circuit
        self.FF_STC = self.Vmpp_STC*self.Impp_STC/(self.Voc_STC*self.Isc_STC)        # fill factor
        
        # Variations with temperature
        self.Ct_Pmax = data_modules[index,7]                                         # Pmax = Pmax(Tref)+Ct(T-Tref) 
        self.Ct_Voc = data_modules[index,8]                                          # Voc = Voc(Tref)+Ct(T-Tref) 
        self.Ct_Isc = data_modules[index,9]                                          # Isc = Isc(Tref)+Ct(T-Tref)  
        
        # Standard test conditions and nominal operationg conditions
        self.tcell_STC = data_modules[index,10] + 273.15                             # Tref
        self.NOCT = data_modules[index,11] + 273.15                                  # Nominal operating cell temperature                                      
       
    def set_geometry(self, beta, gamma, atm = None):
        if self.technology == "fixed":
            self.beta = beta
            self.gamma = gamma
            
        elif self.technology == "tracker":
            """ The panel follows the Sun"""
            self.beta = np.minimum(np.pi/2,atm.sza)
            self.gamma = atm.saa
               
        else:
            print("Unknown panel technology")
        
    def set_SR(self, namelist):
        """Set the spectral response depending on the options chosen"""    
        wvl_astm = tools.wvl_astm()   
        wvl_artdeco = tools.wvl_artdeco()
        
        if namelist.SR_method == "ideal":
            self.SR = tools.spectral_response(wvl_inf = namelist.lambda_inf , wvl_sup = self.lambda_gap, resolution = 'artdeco')             # ARTDECO resolution
            self.SR_astm = tools.spectral_response(wvl_inf = namelist.lambda_inf, wvl_sup = self.lambda_gap, resolution = 'astm')                    # ASTM resolution
            
        elif namelist.SR_method == "real":
            self.SR_astm = tools.spectral_response(celltype = self.type, wvl_inf = namelist.lambda_inf, wvl_sup = 0, resolution='astm', wvl = wvl_astm)
            self.SR = np.interp(wvl_artdeco, wvl_astm, self.SR_astm) 
            
        elif namelist.SR_method == "user":
            data_sr = np.loadtxt("%s/SR_user.dat"%namelist.datadir)
            wvl_user = data_sr[:,0]
            sr_user = data_sr[:,1]
            self.SR_astm = np.interp(wvl_astm, wvl_user, sr_user)
            self.SR = np.interp(wvl_artdeco, wvl_user, sr_user)
            
        # Computing the scaling factor needed to match Isc_STC from provided spectral response
        sw_cell_stc_hres = tools.sw_cell_stc(method_OL = namelist.method_OL, ar = namelist.ar, c1 = namelist.c1, c2 = namelist.c2, n0 = namelist.n0, n1 = namelist.n1, K = namelist.K, L = namelist.L)
        self.alpha_SR = tools.correction_sr(self.SR_astm, self.Isc_STC, spectrum_ref = sw_cell_stc_hres, wvl_ref = wvl_astm)    
        print("For this panel the spectral response is scaled by %s"%self.alpha_SR)       

    def set_rs(self, method_FF = "simple"):
        # Computing expected behaviour in STC to deduce Rs that is supposed constant        
        Voc_STC_norm = self.Voc_STC/self.Vt_STC                                      # Eq. (1) of Green (1982) for one solar cell
        FF0_STC = tools.FF0(Voc_STC_norm)                                            # Eq. (5) of Green (1982), for Rs = 0
      
        # Under the hypothesis that Rs is constant FF is related to FF0 (Green 1982)
        self.Rs = tools.series_resistance(self.FF_STC, FF0_STC, self.Voc_STC, self.Isc_STC, method = method_FF)        

