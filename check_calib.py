#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 16:32:27 2021

@author: pierre
"""

import glob
from astropy.io import fits
import numpy as np
import os
import scipy.signal
from sinfobj import sinfobj
from get_mag import get_coord_sinfoni, get_2MASS_mags
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import numpy.ma as ma
rep_save = '/media/pierre/Disque_2/SNR/CALIBRATED_1/'


def get_obj_header(filename):
    hdu = fits.open(filename)
    obj = hdu[0].data
    header = hdu[0].header
    hdu.close()
    return np.nan_to_num(obj), header

def get_fov(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    pscale = header['CDELT2']*3600
    fov = pscale*np.shape(hdu[0].data)[1]
    hdu.close()
    return fov

def get_coord(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    RA = header['RA']
    dec = header['DEC']
    hdu.close()
    return RA, dec

def get_band(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    wl = header['CRVAL3']
    if wl < 1.3:
        band ='j'
    elif wl < 1.8:
        band='h'
    else:
        band='k'
    return band

def get_zmin_zmax(filter_interp):
    zmin = 0
    zmax = len(filter_interp)
    first_half = True
    done = False
    for x in range(len(filter_interp)):
        # print(x)
        y = filter_interp[x]
        if y == 0 and first_half:
            zmin = x
        elif first_half:
            first_half = False
        elif  y == 0 and (not done):
            zmax = x
            done = True
    return zmin, zmax
            

# for rep1 in glob.glob(rep_save+'/*'):
#     for rep2 in glob.glob(rep1+'/*'):
#         for rep3 in glob.glob(rep2+'/*'):
#             for file in glob.glob(rep3+'/calibrated.fits'):
#                 if 'NGC' in file:
#                     print(file)
#                     obj, header = get_obj_header(file)
#                     raobj, decobj = get_coord(file)
#                     fov = get_fov(file)
#                     band = get_band(file)
#                     mags_obj = np.array(get_2MASS_mags(raobj, decobj, fov=fov, band=band, verbose=False))
#                     mag_obj = np.mean(mags_obj[~np.isnan(mags_obj)])
                    
#                     xmax, ymax = np.unravel_index(np.argmax(scipy.signal.medfilt(np.median(np.nan_to_num(obj),0),(5,5))), np.shape(obj[0]))
    
#                     limits = [np.max([0,xmax-28]), np.min([xmax+28, np.shape(obj[0])[0]]),np.max([0,ymax-28]), np.min([ymax+28, np.shape(obj[0])[1]])]
            
#                     if band == 'j':
#                         F_0 = 3.129e-13*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
#                         bwidth = 0.162
#                         filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/J.txt'
#                         filt = np.loadtxt(filter_name).T
                    
#                     if band == 'h':
#                         F_0 = 3.129e-13*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
#                         bwidth = 0.162
#                         filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/J.txt'
#                         filt = np.loadtxt(filter_name).T
                        
#                     if band =='k':
#                         F_0 = 4.283e-14*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
#                         bwidth = 0.261
#                         filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/K.txt'
#                         filt = np.loadtxt(filter_name).T
                
#                     wl = (np.arange(len(obj))-header['CRPIX3'])*header['CDELT3']+header['CRVAL3']
                
#                     filter_interp = np.interp(wl, filt[0], filt[1])
                    
#                     zmin, zmax = get_zmin_zmax(filter_interp)
                    
#                     obj_2MASSed = obj*filter_interp[:, np.newaxis, np.newaxis]
#                     tot_adu = np.sum(np.nan_to_num(obj_2MASSed[:,limits[0]:limits[1],limits[2]:limits[3]]))
#                     med_adu = np.median(np.nan_to_num(obj_2MASSed[zmin:zmax,limits[0]:limits[1],limits[2]:limits[3]]))*len(obj_2MASSed[zmin:zmax, limits[0]:limits[1], limits[2]:limits[3]].flatten())
#                     if abs(2.5*np.log10(tot_adu/bwidth/F_0)-2.5*np.log10(med_adu/bwidth/F_0))>1.2:
#                         print(band, fov, mag_obj, -2.5*np.log10(tot_adu/bwidth/F_0), -2.5*np.log10(med_adu/bwidth/F_0), '(', med_adu, ' and not ', tot_adu, ')')
#                     else:
#                         print('OK')
                        

filenames = glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/*Circi*/*50/N_*/*.fits')

filenames = glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/*1808*/H*/N_*/*.fits')
filename = '/media/pierre/Disque_2/SNR/archive/076B0098A/reflex_end_products/2021-04-23T11:22:26/SINFO.2005-03-24T00:07:57.876_tpl/calibrated/calibrated.fits'
specs = []
for filename in filenames:
    hdu = fits.open(filename)
    obj = sinfobj(hdu)
    obj.plot_spec(-0.5,0.5,-0.5,0.5, newfig=False)
    spec = obj.get_spec(-0.5,0.5,-0.5,0.5)[1]
    im = obj.get_im(1.6,1.7,-1,1,-1,1)
    xmax, ymax = np.unravel_index(np.argmax(medfilt(im,(5,5))), np.shape(im))
    if xmax > 5:        
        specs.append(ma.masked_invalid(spec))
        print(xmax, ymax)
        print(obj.header['OBJECT'])
    
plt.figure()
plt.plot(np.mean(specs, 0))
