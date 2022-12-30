#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 13:40:38 2021

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
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

from sinfobj import fit_line_cube

def get_coord(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    RA = header['RA']
    dec = header['DEC']
    hdu.close()
    return RA, dec

def get_simbad_info(filename, radius='0d0m10s'):
    name = ' '
    ra, dec = get_coord(filename)
    z = 0
    J = 0
    H = 0
    K = 0
    try:
        customSimbad = Simbad()
        customSimbad.add_votable_fields('z_value')
        customSimbad.add_votable_fields('flux(J)')
        customSimbad.add_votable_fields('flux(H)')
        customSimbad.add_votable_fields('flux(K)')
        result_table = customSimbad.query_region(SkyCoord(ra, dec, unit='deg'), radius=radius)
        x = np.argmin(result_table["FLUX_K"])
        name = result_table["MAIN_ID"][x].decode()
        z = result_table["Z_VALUE"][x]
        J = result_table["FLUX_J"][x]
        H = result_table["FLUX_H"][x]
        K = result_table["FLUX_K"][x]
    except Exception as e:
        print(e)
    return name, z, J, H, K

def radialspeed_exp_disk(r, rd, sigma_0):
    y = r/rd
    A = 4*np.pi*G*sigma_0*rd*y**2
    B1 = iv(0,y)*kn(0,y)
    B2 = iv(1,y)*kn(1,y)
    return (A*(B1-B2))**0.5

def los_velocity_exp_disk(xy, v0, x0, y0, PA, i, rd, sigma_0):
    PA = -PA
    x = xy[0]-x0
    y = xy[1]-y0
    xp = x*np.cos(PA)+y*np.sin(PA)
    yp = -x*np.sin(PA)+y*np.cos(PA)
    # xp0 = x0*np.cos(PA)+y0*np.sin(PA)
    # yp0 = -x0*np.sin(PA)+y0*np.cos(PA)
    r = ((xp)**2+(yp)**2)**0.5
    los_v = radialspeed_exp_disk(r, rd, sigma_0)*xp*np.sin(i)/r
    los_v[r==0]=0
    los_v += v0
    return los_v

def get_random_params(ranges):
    liste = []
    for ran in ranges:
        liste.append(np.random.random()*(ran[1]-ran[0])+ran[0])
    return liste
    

# filenames = glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/*1808*/H*/*/*.fits')
# filename = filenames[0]
# name, z, J, H, K = get_simbad_info(filename)
# pc_per_sec = cosmo.kpc_proper_per_arcmin(z).value*1e3/60
# hdu = fits.open(filename)
# obj = sinfobj(hdu, z=z)
# obj.data[~np.isfinite(obj.data)]=0
# hdu.close()
# data = obj.data[1070:1130,7:61,7:61]
# RA = obj.RA[7:61]
# dec = obj.dec[7:61]
# l = obj.lam[1070:1130]
# cube_p, cube_err, im_snr = fit_line_cube(l, data, mode='gauss')

# kernel = (3, 3)
# pos_med = medfilt(cube_p[:,:,1], kernel)
# v_map = c*(pos_med-np.median(pos_med))/np.median(pos_med)
# v_map = v_map[1:-1,1:-1]
# RA = RA[1:-1]
# dec = dec[1:-1]
# mesh = np.meshgrid(RA,dec)

# f_med = medfilt(cube_p[:,:,0]*(np.pi/cube_p[:,:,2])**0.5, kernel)

# plt.figure()
# plt.imshow(v_map, vmin=-1e5, vmax=1e5)
# plt.figure()
# plt.imshow(f_med)

def fit_los_velocity_exp_disk(xy, v, p0, p_constraint):
    def f_constrained(xy, *args):
        params = []
        k = 0
        for cons in p_constraint:
            if cons != None:
                params.append(cons)
            else:
                params.append(args[k])
                k+=1
        return los_velocity_exp_disk(xy, *params)
    
    
    params_0 = []
    k = 0
    for cons in p_constraint:
        if cons == None:
            params_0.append(p0[k])
        k+=1
        
    p, cov = curve_fit(f_constrained, xy, v, p0 = params_0, maxfev=int(1e6))
    
    full_p = []
    k=0
    for cons in p_constraint:
        if cons != None:
            full_p.append(cons)
        else:
            full_p.append(p[k])
            k+=1
    # plt.figure()
    # plt.imshow(np.reshape(v, (int(len(xy[0])**0.5), int(len(xy[1])**0.5))))
    # plt.figure()
    # plt.imshow(np.reshape(f_constrained(xy, *p), (int(len(xy[0])**0.5), int(len(xy[1])**0.5))))
    return p, cov, full_p

# xy = [mesh[1].flatten()*pc_per_sec*parsec, mesh[0].flatten()*pc_per_sec*parsec]
# v = v_map.flatten()

# p0 = [1e3,1*parsec,1*parsec,-np.pi/4,np.pi/4,1e2*parsec,1e5*M_sun.value/parsec**2]
# p_constraint = [0,0, 0, None, np.pi/4, 1e2*parsec,1e5*M_sun.value/parsec**2]
# p1, cov1 = fit_los_velocity_exp_disk(xy, v, p0, p_constraint)


# p_constraint = [None, 0, 0, p1[0], np.pi/4, None, None]
# p2, cov2 = fit_los_velocity_exp_disk(xy, v, p0, p_constraint)



# print(p1, np.sqrt(np.diag(cov1)))
# print(p2, np.sqrt(np.diag(cov2)))

    
    
    
    
# x = np.arange(-1e3*parsec, 1e3*parsec, 10*parsec)
# y = np.arange(-1e3*parsec, 1e3*parsec, 10*parsec)
# mesh = np.meshgrid(x,y)
# sig0 = 1e5*2e30/parsec**2

# map_los = los_velocity_exp_disk(mesh, 5e2*parsec, 5e2*parsec, 0.5*np.pi/4, np.pi/2, 500*parsec, sig0)
# plt.imshow(map_los)
# # cube_p, cube_err, im_snr = fit_line_cube(l, data, mode='gauss')
# pos = cube_p[:,:,1]
# v_map = c*(pos-np.median(pos))/np.median(pos)
# mask = np.ones(np.shape(im_snr), dtype='bool')
# mask[6:62,6:62]=False

# v = np.ma.masked_array(v, mask)
# plt.imshow(v, vmin=-1e5, vmax=1e5)

# mesh = np.meshgrid(obj.RA*50*parsec, obj.dec*50*parsec)

# x0, y0 = np.unravel_index(np.argmax(medfilt(cube_p[:,:,0],(5,5))), np.shape(im_snr)) 
# ra0 = obj.RA[x0]
# dec0 = obj.dec[y0]


# p0=[0,0, -np.pi, np.pi/4, 1e3*parsec, 100]
# p0 =[0,0, 5*np.pi/4, np.pi/4, 1e2*parsec, 100]

# def g(mesh, PA):
#     p0=[1e1*parsec,0, PA, np.pi/4, 1e2*parsec, 100]
#     return los_velocity_exp_disk(mesh, *p0)

# p1, c = curve_fit(g, np.array([xp, yp]), zp)

# def h(mesh, x0, y0, i, rd, sigma_0):
#     p0=[x0,y0, p1[0], i, rd, sigma_0]
#     return los_velocity_exp_disk(mesh, *p0)

# p2, c = curve_fit(h, np.array([xp, yp]), zp, p0=[ra0*50*parsec, dec0*50*parsec, np.pi/4, 1e2, 100])

# def k(mesh, rd, sigma0):
#     p0 = [ra0*50*parsec, dec0*50*parsec, p1[0], np.pi/4, rd, sigma0]
#     return los_velocity_exp_disk(mesh, *p0)

# p3, c = curve_fit(k, np.array([xp, yp]), zp, p0=[1e1*parsec, 1e2])

# p, c = curve_fit(los_velocity_exp_disk, np.array([xp, yp]), zp, p0=p0, maxfev=int(1e6))


# ranges = [[-2e2*parsec, 2e2*parsec], [-2e2*parsec, 2e2*parsec], [0, 2*np.pi], [0, np.pi], [0,1e3*parsec], [0, 1e4]]

# ranges = [[-2e0*parsec, 2e0*parsec], [-2e0*parsec, 2e0*parsec], [0, 2*np.pi], [0, np.pi], [0,1e3*parsec], [0, 1e4]]

# xp = mesh[0][mask]
# yp = mesh[1][mask]
# zp = v_map[mask]-np.median(v_map[mask])
# khi2_min=1e50
# t0 = time.time()
# for k in range(int(1e5)):
#     p = get_random_params(ranges)
#     mod = los_velocity_exp_disk(np.array([xp, yp]), *p)
#     khi2 = np.sqrt(np.sum((zp-mod)**2))
#     if khi2 < khi2_min:
#         print(khi2, p)
#         best_params = p
#         khi2_min = khi2
# print(time.time()-t0)
# plt.figure()
# plt.subplot(121)
# im_v_mod = np.ma.masked_array(los_velocity_exp_disk(mesh, *best_params), mask)
# plt.imshow(im_v_mod, vmin=-1e5, vmax=1e5)
# plt.subplot(122)
# im_v = np.ma.masked_array(v_map-np.median(v_map[mask]), mask)
# plt.imshow(im_v, vmin=-1e5, vmax=1e5)
    
