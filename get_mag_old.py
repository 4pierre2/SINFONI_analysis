#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 13:39:27 2021

@author: pierre
"""

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from glob import glob

im_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/1808/aK_asky_981231s0540080.fits'
#xmax, ymax = 429, 854
im_name = glob('/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/NGC5253/*K_*.fits')[0]
xmax, ymax = 321, 397
im_name = glob('/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/NGC7552/aK_asky_990916s0970245.fits')[0]
xmax, ymax = 56, 115
im_name = glob('/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/M83/aK_*.fits')[0]
xmax, ymax = 269, 652
im_name = glob('/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/NGC6240/aJ*s071*.fits')[0]
xmax, ymax = 420, 500
# im_name = glob('/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/NGC1068/aJ_*4s068*.fits')[0]
# xmax, ymax = 71, 433
# im_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/HIP_026080/aJ_askyw_990208s0140198.fits'
# im_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/HIP_026080/aH_askyw_990208s0140198.fits'
# im_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_ims/BD-00413/aJ_asky_981004s0750267.fits'
hdu = fits.open(im_name)
im = hdu[0].data
w = WCS(hdu[0].header)
zp = hdu[0].header['MAGZP']

pixscale = 3600*abs(hdu[0].header['CDELT1'])

def phot_2mass(im,x0,y0,r=14,rin=14,rout=20):
    y = np.arange(len(im[:,0]))
    x = np.arange(len(im[0,:]))
    xx, yy = np.meshgrid(x,y)
    radius = ((xx-x0)**2+(yy-y0)**2)**0.5
    mask = radius<(r/pixscale)
    mask_a = radius<(rin/pixscale)
    mask_b = np.array((radius<(rout/pixscale))*1.-1.*mask_a, dtype='bool')
    im_no_bg = im-np.mean(np.nan_to_num(im[mask_b]))
    phot = np.sum(im_no_bg[mask])
    mag = -2.5*np.log10(phot)+zp
    return mag

def phot_sinf(im,x0,y0,fov=7):
    y = np.arange(len(im[:,0]))
    x = np.arange(len(im[0,:]))
    xx, yy = np.meshgrid(x,y)
    mask = ((((xx-x0)**2)**0.5)<=fov/2) *((((yy-y0)**2)**0.5)<=fov/2)
    im_no_bg = im-np.median(np.nan_to_num(im))
    phot = np.sum(im_no_bg[mask])
    
    mag = -2.5*np.log10(phot)+zp
    return mag

print(phot_sinf(im, xmax, ymax, fov=3))





