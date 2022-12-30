#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:51:01 2021

@author: pierre
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
# J band
F_0 = 3.129e-13*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
bwidth = 0.162
filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/J.txt'
filt = np.loadtxt(filter_name).T

# NGC 1808
mag = 10.1
coadd_name = '/media/pierre/Disque_2/SNR/SORTED_COADD/NGC1808/J/0.25/NGC1808-KHJ_COADD_OBJ_53463.00876717.fits'
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T00:52:17.366_tpl/HIP-028060-telluric-for-NGC1808_STD_STAR_SPECTRA.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/J.fits'


# NGC 5253
mag = 12.61
coadd_name = glob.glob('/media/pierre/Disque_2/SNR/SORTED_COADD/NGC5253/J/0.25/*COADD_OBJ*.fits')[1]
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T08:22:55.934_tpl/HIP-070243-telluric-for-NGC5253_STD_STAR_SPECTRA.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_5253/J.fits'

# NGC 7552
mag = 10.65
coadd_name = glob.glob('/media/pierre/Disque_2/SNR/SORTED_COADD/NGC7552/J/0.25/*COADD_OBJ*.fits')[2]
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T08:22:55.934_tpl/HIP-070243-telluric-for-NGC5253_STD_STAR_SPECTRA.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_7552/J.fits'



# H band
F_0 = 1.133e-13*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
bwidth = 0.251
filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/H.txt'
filt = np.loadtxt(filter_name).T


# NGC 1808
mag = 9.2
coadd_name = '/media/pierre/Disque_2/SNR/SORTED_COADD/NGC1808/H/0.25/NGC1808-KHJ_COADD_OBJ_53462.99568248.fits'
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T00:45:55.796_tpl/HIP-028060-telluric-for-NGC1808_STD_STAR_SPECTRA.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/H.fits'

# NGC 5253
mag = 12.09
coadd_name = glob.glob('/media/pierre/Disque_2/SNR/SORTED_COADD/NGC5253/H/0.25/*COADD_OBJ*.fits')[1]
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T08:18:04.338_tpl/HIP-070243-telluric-for-NGC5253_STD_STAR_SPECTRA.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_5253/H.fits'

# NGC 7552
mag = 9.87
coadd_name = glob.glob('/media/pierre/Disque_2/SNR/SORTED_COADD/NGC7552/H/0.25/*COADD_OBJ*.fits')[2]
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T08:18:04.338_tpl/HIP-070243-telluric-for-NGC5253_STD_STAR_SPECTRA.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_7552/H.fits'





# K band
F_0 = 4.283e-14*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
bwidth = 0.261
filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/K.txt'
filt = np.loadtxt(filter_name).T


# NGC 1808
mag = 8.8
coadd_name = '/media/pierre/Disque_2/SNR/SORTED_COADD/NGC1808/K/0.25/NGC1808-KHJ_COADD_OBJ_53462.98342928.fits'
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T00:41:54.979_tpl/HIP-028060-telluric-for-NGC1808_STD_STAR_SPECTRA.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/K.fits'

# NGC 5253
mag = 11.55
coadd_name = glob.glob('/media/pierre/Disque_2/SNR/SORTED_COADD/NGC5253/K/0.25/*COADD_OBJ*.fits')[0]
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T08:25:57.761_tpl/HIP-070243-telluric-for-NGC5253_STD_STAR_SPECTRA.fits'
spec_skyname = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T08:25:57.761_tpl/HIP-070243-telluric-for-NGC5253_OBS_SKY.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_5253/K.fits'

# NGC 7552
mag = 9.47
coadd_name = glob.glob('/media/pierre/Disque_2/SNR/SORTED_COADD/NGC7552/K/0.25/*COADD_OBJ*.fits')[0]
spec_star_filename = '/media/pierre/Disque_2/SNR/archive/075B0648A/reflex_end_products/2021-04-15T23:48:58/SINFO.2005-04-03T08:25:57.761_tpl/HIP-070243-telluric-for-NGC5253_STD_STAR_SPECTRA.fits'
typ = 'b8v'
savename='/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_7552/K.fits'



# pixs = xmax-28:xmax+28,ymax-28:ymax+28

def get_obj(filename):
    hdu = fits.open(filename)
    obj = hdu[0].data
    return np.nan_to_num(obj)
    

def get_info(filename):
    hdu = fits.open(filename)
    prog_id = hdu[0].header['HIERARCH ESO OBS PROG ID'].replace('.','').replace('-','').replace('(','').replace(')','')
    filt = hdu[0].header['HIERARCH ESO INS GRAT1 NAME'].replace(' ','')
    pix_scale = hdu[0].header['HIERARCH ESO INS OPTI1 NAME'].replace(' ','')
    date = hdu[0].header['MJD-OBS']
    return prog_id, date, filt, pix_scale

def get_obj_info(filename):
    hdu = fits.open(filename)
    prog_id = hdu[0].header['HIERARCH ESO OBS PROG ID'].replace('.','').replace('-','').replace('(','').replace(')','')
    obj_name = hdu[0].header['HIERARCH ESO OBS TARG NAME']
    date = hdu[0].header['MJD-OBS']
    return obj_name, date


def find_std_stars(filename, rep='/media/pierre/Disque_2/SNR/archive'):
    prog_id, date, filt, pix_scale = get_info(filename)
    rep_of_interest = rep+'/'+prog_id+'/reflex_end_products'
    for rep in glob.glob(rep_of_interest+'/*'):
        print(rep)
        for re in glob.glob(rep+'/*'):
            star_spectra = glob.glob(re+'/*STD_STAR_SPECTRA.fits')
            for spec in star_spectra:
                std_id, date, std_filt, std_pix = get_info(spec)
                if std_id == prog_id and std_filt == filt and std_pix == pix_scale:
                    print(star_spectra, get_obj_info(spec), get_info(spec))

def load_th_spec(typ, rep='/home/pierre/Documents/2021/SINFONI_analysis/pickles_atlas/'):
    full = np.loadtxt(rep+'/uk'+typ+'.dat')
    return full[:,0], full[:,1], full[:,2]

find_std_stars(coadd_name)


spec_star_wl = fits.open(spec_star_filename)[1].data['wavelength']
spec_star_spectrum = fits.open(spec_star_filename)[1].data['counts_bkg']
spec_star_bg = fits.open(spec_star_filename)[1].data['bkg_tot']
spec_star_tot = fits.open(spec_star_filename)[1].data['counts_tot']


w, f, e = load_th_spec(typ)
spec_th = np.interp(spec_star_wl, w/1e4, f)

transmi = np.nan_to_num(spec_star_spectrum/spec_th)
transmi /= np.median(transmi)
transmi[transmi<5e-2]=5e-2


obj = get_obj(coadd_name)
xmax, ymax = np.unravel_index(np.argmax(np.median(np.nan_to_num(obj),0)), np.shape(obj[0]))
obj_detransmitted = obj/transmi[:, np.newaxis, np.newaxis]

filter_interp = np.interp(spec_star_wl, filt[0], filt[1])
obj_2MASSed = obj_detransmitted*filter_interp[:, np.newaxis, np.newaxis]
tot_adu = np.sum(np.nan_to_num(obj_2MASSed[:,xmax-28:xmax+28,ymax-28:ymax+28]))
tot_flux = F_0*10**(-0.4*mag)*bwidth
adu = tot_flux/tot_adu
dwl = spec_star_wl[2]-spec_star_wl[1]
obj_final = obj_detransmitted*adu

# hdu = fits.open(coadd_name)
# hdu[0].data = obj_final

# hdu.writeto(savename)