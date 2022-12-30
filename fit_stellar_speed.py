#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 11:54:03 2021

@author: pierre
"""


import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from scipy.signal import medfilt
from sinfobj import sinfobj
from copy import deepcopy
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from scipy.constants import c, pi, G, parsec
from scipy.optimize import curve_fit
from astropy.constants import M_sun
from scipy.constants import G, parsec, c
from scipy.special import iv, kn
from astropy.cosmology import FlatLambdaCDM
# import pysynphot as S

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

hdu_1808_k = fits.open('/media/pierre/Disque_2/SNR/CLEANED_1/NGC1808_K_125/cleaned.fits')
obj = sinfobj(hdu_1808_k)

list_of_spectra = glob.glob('/home/pierre/Documents/STELLAR_SPECTRA/SAMPLE/*.fits')

spectra = []
wls = []
names = []
for file in list_of_spectra:
    hdu = fits.open(file)
    data = hdu[0].data
    header = hdu[0].header
    wl = np.arange(len(data))+header['CRVAL1']
    spectra.append(data)
    wls.append(wl)
    
spectra_interps = []
for spectrum, wl in zip(spectra, wls):
    print(len(wl), len(spectrum))
    spectrum_interp = np.interp(obj.lam*10000, wl, spectrum)
    spectra_interps.append(spectrum_interp)
    
    

