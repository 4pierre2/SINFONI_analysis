#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 31 11:12:36 2021

@author: pierre
"""

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
from scipy.signal import medfilt
from sinfobj import sinfobj
from copy import deepcopy
from scipy.constants import c, pi, G, parsec
from scipy.optimize import curve_fit
from astropy.constants import M_sun
from matplotlib.backends.backend_pdf import PdfPages
import scipy.signal 
G2 = G*parsec**3/M_sun.value

z = 1.003343
Dl = 12.8
pix_scale = 0.125
def A_l(wl, Av, Rv=4.05):
    return Av/Rv*(2.659*(-1.857+1.040/wl)+Rv)

def multipage(filename, figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()
    

def manypdf(filename, figs=None, dpi=200, forma='pdf'):
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    k=0
    for fig in figs:
        pp = PdfPages(filename+'_'+str(k)+forma)
        fig.savefig(pp, format=forma)
        pp.close()
        k+=1
    
def manypng(filename, figs=None, dpi=200):
    forma = 'png'
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    k=0
    for fig in figs:
        fig.savefig(filename+'_'+str(k)+'.'+forma, format=forma)
        k+=1
        
hdu_1808_j= fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/J.fits')
hdu_1808_h = fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/H.fits')
hdu_1808_k = fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/K.fits')
ngc1808_j = sinfobj(hdu_1808_j)
ngc1808_h = sinfobj(hdu_1808_h)
ngc1808_k = sinfobj(hdu_1808_k)


plt.figure(figsize=(12,6))
ngc1808_j.plot_spec(-1,1,-1,1,newfig=False)
ngc1808_h.plot_spec(-1,1,-1,1, newfig=False)
ngc1808_k.plot_spec(-1,1,-1,1, newfig=False)
plt.axis([1.11,2.45,0,9e-14])



a18, b18, c18, d18, mask18 = ngc1808_j.fit_all_gauss(1.286,1.285,1.29)
ab18, bb18, cb18, db18, maskb18= ngc1808_k.fit_all_gauss(2.173,2.170,2.178)
ext1808 =  11.52*np.log10(5.89/(a18/ab18))


mask1808 = np.zeros(np.shape(ext1808))

mask1808[:,:6]=1
mask1808[:,61:]=1
# inds = np.array([[59, 17],[58, 17],[51, 43],[50, 43],[43, 55],[42, 55],[41, 56],[40, 56],[17, 59],[16, 59]])
# for ind in inds:
#     print(ind[0], ind[1])
#     mask1808[ind[0], ind[1]]=1
# inds = np.array([[23, 29], [22, 29], [63,  9],[62,  9],[53, 12],[52, 12],[55, 17],[54, 17],[17, 42],[16, 42],[17, 51],[16, 51],[11, 52],[10, 52],[ 5, 48],[ 4, 48],[ 3, 26],[ 2, 26],[ 3, 27],[ 2, 27],[ 7, 38],[ 6, 38],[ 5, 48],[ 4, 48]])
# for ind in inds:
#     print(ind[0], ind[1])
#     mask1808[ind[0], ind[1]]=1
    
    
    
map1808 = np.ma.masked_array(ext1808, mask1808)

map1808 = scipy.signal.medfilt(map1808,(3,3))
map1808 = np.ma.masked_array(map1808, mask1808)
# plt.imshow(map1808, vmin=0, vmax=5, cmap='jet')
plt.figure()
plt.imshow(map1808, origin='lower', vmin=0, vmax=6, cmap='jet')
plt.xlabel('Relative Right Ascension (")')
plt.ylabel('Relative Declination (")')
cbar = plt.colorbar()
cbar.set_label('$A_v$')

# plt.colorbar()



a_ab = 5.89/10**(7/11.52)

 
ah218, bh218, ch218, dh218, maskh218 = ngc1808_k.fit_all_gauss(2.130,2.120,2.132)

plt.figure()
map_h2flux = np.ma.masked_array(10**(0.4*A_l(2.12,ext1808))*ah218*(np.pi/ch218)**0.5, mask1808)
plt.imshow(map_h2flux, cmap='jet', vmin = 8e-21, vmax=1e-19, origin='lower', norm=LogNorm())
plt.xlabel('Relative Right Ascension (")')
plt.ylabel('Relative Declination (")')
cbar = plt.colorbar()
cbar.set_label('Flux ($W.m^{-2}.pixel^{-1}$)')


plt.figure()
map_h2hotmass = 5.0776e16*Dl**2*map_h2flux/pix_scale**2
plt.imshow(map_h2hotmass, cmap='jet', origin='lower')
plt.xlabel('Relative Right Ascension (")')
plt.ylabel('Relative Declination (")')
cbar = plt.colorbar()
cbar.set_label('Mass ($M_\odot.arcsec^{-1}$)')

# plt.figure()
map_h2lam = np.ma.masked_array(bh218, mask1808)
# img = plt.imshow(map_h2lam, cmap='jet', vmin = 2.1282, vmax=2.1302, origin='lower')
map_h2speed = c*(map_h2lam/z-2.122)/2.122
plt.figure()
plt.imshow(map_h2speed*1e-3, cmap='jet', origin='lower', vmin=-200, vmax=200, extent=[np.min(ngc1808_k.RA), np.max(ngc1808_k.RA), np.min(ngc1808_k.dec), np.max(ngc1808_k.dec)])
plt.xlabel('Relative Right Ascension (")')
plt.ylabel('Relative Declination (")')
cbar = plt.colorbar()
cbar.set_label('Velocity ($km.s^{-1}$)')
# as18, bs18, cs18, ds18, masks18 = ngc1808_k.fit_all_gauss_abs(2.330,2.325,2.345)

# map_slam = np.ma.masked_array(bs18, mask1808)
# map_sspeed = c*(map_slam/z-2.3227)/2.3227
# plt.figure()
# plt.imshow(map_sspeed*1e-3, cmap='jet', origin='lower', vmin=-200, vmax=200)
# plt.colorbar()


def plummer_full(mesh, Vs, M, A, i0, psi0, x0, y0):
    if i0>2*np.pi:
        return 1e25*i0*mesh[0]
    if i0<0*np.pi:
        return 1e25*i0*mesh[0]
    if psi0>2*np.pi:
        return 1e25*i0*mesh[0]
    if psi0<0*np.pi:
        return 1e25*i0*mesh[0]
    xs = mesh[0]
    ys = mesh[1]
    R = ((xs-x0)**2+(ys-y0)**2)**0.5
    psi = np.arctan2(ys, xs)
    psi[(xs==0)*(ys==0)] = 0
    T1 = Vs
    T2 = ((R**2*G2*M)/(R**2+A**2)**1.5)**0.5
    T3 = (np.sin(i0)*np.cos(psi-psi0))/((np.cos(psi-psi0)**2+np.sin(psi-psi0)**2/np.cos(i0)**2)**0.75)
    resulting_map = np.array(T1+T2*T3, dtype='float')
    return resulting_map.flatten()#*0+psi.flatten()


def plummer(mesh, Vs, M, A, i0, psi0):
    if i0>2*np.pi:
        return 1e25*i0*mesh[0]
    if i0<0*np.pi:
        return 1e25*i0*mesh[0]
    if psi0>2*np.pi:
        return 1e25*i0*mesh[0]
    if psi0<0*np.pi:
        return 1e25*i0*mesh[0]
    xs = mesh[0]
    ys = mesh[1]
    R = ((xs)**2+(ys)**2)**0.5
    psi = np.arctan2(ys, xs)
    psi[(xs==0)*(ys==0)] = 0
    T1 = Vs
    T2 = ((R**2*G2*M)/(R**2+A**2)**1.5)**0.5
    T3 = (np.sin(i0)*np.cos(psi-psi0))/((np.cos(psi-psi0)**2+np.sin(psi-psi0)**2/np.cos(i0)**2)**0.75)
    resulting_map = np.array(T1+T2*T3, dtype='float')
    return resulting_map.flatten()#*0+psi.flatten()


mesh = np.meshgrid(ngc1808_k.RA*62, ngc1808_k.dec*62 )
m = mask1808
pf, cf = curve_fit(plummer_full, [mesh[0][m==0],mesh[1][m==0]], map_h2speed[m==0], p0 = [0, 1e5, 10, 50*np.pi/180, 130*np.pi/180, 0, 0], maxfev=10000)

p, c = curve_fit(plummer, [mesh[0][m==0],mesh[1][m==0]], map_h2speed[m==0], p0 = [0, 1e5, 100, 50*np.pi/180, 130*np.pi/180], maxfev=10000)

plt.figure()
plt.imshow(np.reshape(plummer_full(mesh, *pf)*1e-3, np.shape(mesh[0])), cmap='jet', origin='lower', vmin=-200, vmax=200, extent=[np.min(ngc1808_k.RA), np.max(ngc1808_k.RA), np.min(ngc1808_k.dec), np.max(ngc1808_k.dec)])
plt.xlabel('Relative Right Ascension (")')
plt.ylabel('Relative Declination (")')
cbar = plt.colorbar()
cbar.set_label('Velocity ($km.s^{-1}$)')
# res_min=1e30
# for M in 10**np.arange(10):
#     for A in 5**np.arange(5):
#         for i0 in np.arange(0,2*np.pi,np.pi/8):
#             for psi0 in np.arange(0,2*np.pi,np.pi/8):
#                 rm = plummer([mesh[0][m==0],mesh[1][m==0]], 0, M, A, i0, psi0)
#                 res = np.sum((rm-map_h2speed[m==0])**2)
#                 if res < res_min:
#                     res_min = res
#                     print(M, A, i0, psi0)
    


mask1808_h = np.zeros(np.shape(ext1808))

mask1808_h[:,:6]=1
mask1808_h[:,61:]=1
# inds = np.array([[59, 17],[58, 17],[51, 43],[50, 43],[43, 55],[42, 55],[41, 56],[40, 56],[17, 59],[16, 59]])
# for ind in inds:
#     print(ind[0], ind[1])
#     mask1808[ind[0], ind[1]]=1

afe218, bfe218, cfe218, dfe218, maskfe18 = ngc1808_h.fit_all_gauss(1.65, 1.645, 1.655)

pir2 = 4*np.pi*(Dl*1e6*parsec)**2

snr = 8.08*10**(0.4*A_l(1.644,ext1808))*afe218*(np.pi/cfe218)**0.5*pir2*1e-35

snr = np.ma.masked_array(snr, mask1808_h)
plt.figure()
plt.imshow(1e5*snr, origin='lower', vmin=0, vmax=10, extent=[np.min(ngc1808_h.RA), np.max(ngc1808_h.RA), np.min(ngc1808_h.dec), np.max(ngc1808_h.dec)])
plt.xlabel('Relative Right Ascension (")')
plt.ylabel('Relative Declination (")')
cbar = plt.colorbar()
cbar.set_label('Supernova rate ($10^{-5}\ yr^{-1}.pixel^{-1}$)')


def multipage(filename, figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()
    

def manypdf(filename, figs=None, dpi=200, forma='pdf'):
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    k=0
    for fig in figs:
        pp = PdfPages(filename+'_'+str(k)+'.'+forma)
        fig.savefig(pp, format=forma)
        pp.close()
        k+=1
    
def manypng(filename, figs=None, dpi=200):
    forma = 'png'
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    k=0
    for fig in figs:
        fig.savefig(filename+'_'+str(k)+'.'+forma, format=forma)
        k+=1
multipage('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/ANALYSIS/plots.pdf')
manypdf('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/ANALYSIS/plot')
manypng('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/ANALYSIS/plot')
plt.close('all')
