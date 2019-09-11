#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:37:12 2019

@author: pallavipatil
"""


import numpy as np
import pandas as pd
from Radio_Models_func import radio_models_fits
from file_prep_radio_fitting import load_tables, file_prep_utils

#from Radio_Models_func import radio_models as rm
#from matplotlib.widgets import Slider, Button, RadioButtons, Cursor
#from astropy.io import fits,ascii
#import matplotlib.pyplot as plt
#from astropy.table import Table, unique, join
#import os
#from astropy.coordinates import SkyCoord
#import astropy.units as u
#import aplpy as apl
#from astropy.visualization import ZScaleInterval
#import scipy


##Accessing all different data files and creating astropy tables.
    
'''Edit the path here and make sure all the following files are in your path.
WISE_Xmatch.fits:: This is table Carol provided and has information from WISE to low
radio frequencies, mainly archival information.
JMFIT_CASA_A_results.dat:: Results from JVLA AX observations.
JMFIT_CASA_B_results.dat:: Results from JVLA BX observations. 
'''



loadt = load_tables()
fprep = file_prep_utils()
rmfit = radio_models_fits()
atscat, atgcat, vla_ax_grp, vla_bx_grp = loadt.get_tabs()

'''
plots_dir = './Plots_all/'
fitsdir_AX = '/Users/pallavipatil/Desktop/VLA/VLA-NewWork/AX/fitsfiles/'
fitsdir_BX = '/Users/pallavipatil/Desktop/VLA/VLA-NewWork/BX/fitsfiles/'

scan_info = pd.ExcelFile('JVLA_Scans_info.xlsx')
dfsinfo = pd.read_excel(scan_info, 'Summary')
scan_info.close()
ax_list = pd.Series(dfsinfo['12B-217'])
bx_list = pd.Series(dfsinfo['12A-064'])
'''

sp_class = loadt.get_sp_info()

#fitTab = ascii.read('GuessPar_Radiofits.csv', format='basic', delimiter=',')
#fitTab.add_index('Name')




for wisen in atscat['WISEname'][:10]:
    best_fits = []
    chi_sqs = []
    par_err = []

    scat = atscat.loc[wisen]
    glmcat = atgcat.loc[wisen]
    jvla_AX = vla_ax_grp.groups[vla_ax_grp.groups.keys['WISEname'] == wisen]
    jvla_BX = vla_bx_grp.groups[vla_bx_grp.groups.keys['WISEname'] == wisen]
    freq_arr, flux_arr, eflux_arr, alpha_AX, alpha_BX, alpha_GL, ALMA, \
            OIR, sp_flux, labels= fprep.data_prep(wisen, glmcat, scat,  jvla_AX, jvla_BX)
    for alpha in [alpha_AX, alpha_BX, alpha_GL]:
        freq_arr, flux_arr, eflux_arr = rmfit.prep_fit_arr(freq_arr, flux_arr, 
                                                           eflux_arr, alpha)
    #Guess parameters:
    guess_pars = fitTab.loc[wisen]    
    models = ['PL', 'CPL', 'EFFA', 'IFFA', 'SSA']
    s0 = np.max(flux_arr)
    nu_t = np.min(freq_arr[freq_arr>0])
    alpha = -0.7
    guess_cpl = [s0, alpha,nu_t]
    guess_pl = [s0, alpha]
    # Perform radio fitting:
    for model in models:
        if model =='PL':
            guess_pars = guess_pl
        else:
            guess_pars = guess_cpl
        fit_res, perr_res, chi_res = rmfit.chi_sq_fit(freq_arr, flux_arr, eflux_arr,
                                                      guess_pars, model)
        best_fits.append(fit_res)
        par_err.append(perr_res)
        chi_sqs.append(chi_res)
        
    print(best_fits)
        
    
    
                
    
    
