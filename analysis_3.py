#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 09:23:41 2021

@author: pierre
"""


import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from sinfobj import sinfobj

hdu_1808_j= fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/J.fits')
hdu_1808_h = fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/H.fits')
hdu_1808_k = fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/K.fits')
ngc1808_j = sinfobj(hdu_1808_j)
ngc1808_h = sinfobj(hdu_1808_h)
ngc1808_k = sinfobj(hdu_1808_k)

hdu_7552_j= fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_7552/J.fits')
hdu_7552_h = fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_7552/H.fits')
hdu_7552_k = fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_7552/K.fits')
ngc7552_j = sinfobj(hdu_7552_j)
ngc7552_h = sinfobj(hdu_7552_h)
ngc7552_k = sinfobj(hdu_7552_k)

# lj, sj = ngc7552_j.get_spec(-2,2,-2,2)
# lh, sh = ngc7552_h.get_spec(-2,2,-2,2)
# lk, sk = ngc7552_k.get_spec(-2,2,-2,2)

# lj, sj = ngc1808_j.get_spec(-2,2,-2,2)
# lh, sh = ngc1808_h.get_spec(-2,2,-2,2)
# lk, sk = ngc1808_k.get_spec(-2,2,-2,2)

# plt.plot(np.concatenate([lj, lh, lk]), np.concatenate([sj, sh, sk]))

a18, b18, c18, d18, mask18 = ngc1808_j.fit_all_gauss(1.286,1.285,1.29)
ab18, bb18, cb18, db18, maskb18= ngc1808_k.fit_all_gauss(2.173,2.170,2.178)
ext1808 =  11.52*np.log10(5.89/(a18/ab18))

mask1808 = np.zeros(np.shape(ext1808))

mask1808[:,:6]=1
mask1808[:,61:]=1
inds = np.array([[59, 17],[58, 17],[51, 43],[50, 43],[43, 55],[42, 55],[41, 56],[40, 56],[17, 59],[16, 59]])
for ind in inds:
    print(ind[0], ind[1])
    mask1808[ind[0], ind[1]]=1
# mask1808[ext1808>7] = 1
# mask1808[ext1808<0] = 1

map1808 = np.ma.masked_array(ext1808, mask1808)
img = plt.imshow(map1808, origin='lower', vmin=0, vmax=6, cmap='jet')
plt.imshow(ext1808, origin='lower', vmin=0, vmax=6, cmap='jet')
plt.colorbar()



a, b, c, d, mask = ngc7552_j.fit_all_gauss(1.286,1.285,1.29)
ab, bb, cb, db, maskb = ngc7552_k.fit_all_gauss(2.173,2.171,2.182)
ext7552 =  11.52*np.log10(5.89/(a/ab))

plt.figure()
plt.imshow(ext7552, origin='lower', vmin=0, vmax=7, cmap='jet')
plt.colorbar()


a_ab = 5.89/10**(7/11.52)

def A_l(wl, Av, Rv=4.05):
    return Av/Rv*(2.659*(-1.857+1.040/wl)+Rv)
    
ah218, bh218, ch218, dh218, maskh218 = ngc1808_k.fit_all_gauss(2.130,2.120,2.140)
ah275, bh275, ch275, dh275, maskh275 = ngc7552_k.fit_all_gauss(2.130,2.120,2.140)

plt.figure()
plt.imshow(10**(0.4*A_l(2.12,ext1808))*ah218*(np.pi/ch218)**0.5, cmap='jet', vmin = 0, vmax=1.5e-19, origin='lower')
plt.colorbar()
plt.figure()
plt.imshow(bh218, cmap='jet', vmin = 2.1282, vmax=2.1302, origin='lower')
plt.colorbar()
plt.figure()
plt.imshow(10**(0.4*A_l(2.12,ext7552))*ah275*(np.pi/ch275)**0.5, cmap='jet', vmin = 0, vmax=3.2e-20, origin='lower')
plt.colorbar()
plt.figure()
plt.imshow(bh275, cmap='jet', vmin = 2.1327, vmax=2.1341, origin='lower')
plt.colorbar()

# ngc75.data /= 
