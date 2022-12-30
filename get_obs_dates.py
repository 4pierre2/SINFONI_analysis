#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 12:20:20 2021

@author: pierre
"""


import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from sinfobj import sinfobj


rep = '/media/pierre/Disque_2/SNR/archive/'

dirs = glob.glob(rep+'/*')

hdus = []
dates =  []
progs_ids = []
for di in dirs:
    prog_ids = []
    prog_dates = []
    print(di)
    fits_of_interest = glob.glob(di+'/*.fits')
    for fit in fits_of_interest:
        hdu = fits.open(fit)
        try:
            dates.append(hdu[0].header['DATE-OBS'][:10])
            prog_dates.append(hdu[0].header['DATE-OBS'][:10])
            progs_ids.append(hdu[0].header['HIERARCH ESO OBS PROG ID'])
            prog_ids.append(hdu[0].header['HIERARCH ESO OBS PROG ID'])
        except:
            b = 0
    print(list(set(prog_ids)))
    print(sorted(list(set(prog_dates))))
            
        
# dates =  []
# for hdu in hdus:
#     dates.append(hdu[0].header['DATE-OBS'][:10])
    
print(list(set(dates)))
print(list(set(progs_ids)))