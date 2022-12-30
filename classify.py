#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:29:56 2021

@author: pierre
"""

import os
import glob
import shutil
from astropy.io import fits

def makedir(path):
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except:
            print('Error creating : '+path)
            
def make_sorted_dir(obj, grating, pix_scale, path = './SORTED_FITS/'):
    makedir(path)
    makedir(path+'/'+obj)
    makedir(path+'/'+obj+'/'+grating)
    makedir(path+'/'+obj+'/'+grating+'/'+pix_scale)
    return path+'/'+obj+'/'+grating+'/'+pix_scale+'/'

rep = '/media/pierre/Disque_2/SNR/archive/'

dirs = glob.glob(rep+'/0*')

hdus = []
dates =  []
progs_ids = []
for di in dirs:
    di = di+'/reflex_end_products/'
    dib = glob.glob(di+'/*')
    subdirs = glob.glob(dib[-1]+'/*')
    for dia in subdirs:
        files = glob.glob(dia+'/*COADD_OBJ*')
        for file in files:
            filename = file.split('/')[-1]
            hdu = fits.open(file)
            header = hdu[0].header
            mjd = str(header['MJD-OBS']).replace(' ','')
            obj = header['OBJECT'].replace(' ','')
            grating = header['HIERARCH ESO INS FILT1 NAME'].replace(' ','')
            pix_scale = header['HIERARCH ESO INS OPTI1 NAME'].replace(' ','')
            path = make_sorted_dir(obj, grating, pix_scale, path='/media/pierre/Disque_2/SNR/SORTED_COADD')
            dest_file = path+filename.replace(".fits", '_'+mjd+'.fits')
            if not os.path.exists(path+filename.replace(".fits", '_'+mjd+'.fits')):
                shutil.copy(file, dest_file)
            else:
                print('File already exists')
            
rep = '/media/pierre/Disque_2/SNR/SORTED_COADD/'
dirs = glob.glob(rep+'/*')
names = []
for di in dirs:
    name = di.split('/')[-1]
    names.append(name)
print(sorted(list(set(names))))

    # prog_ids = []
    # prog_dates = []
    # fits_of_interest = glob.glob(di+'/*.fits')
    # for fit in fits_of_interest:
    #     hdu = fits.open(fit)
    #     try:
    #         dates.append(hdu[0].header['DATE-OBS'][:10])
    #         prog_dates.append(hdu[0].header['DATE-OBS'][:10])
    #         progs_ids.append(hdu[0].header['HIERARCH ESO OBS PROG ID'])
    #         prog_ids.append(hdu[0].header['HIERARCH ESO OBS PROG ID'])
            
    #     except:
    #         b = 0
    # print(list(set(prog_ids)))
    # print(sorted(list(set(prog_dates))))
            
        
# dates =  []
# for hdu in hdus:
#     dates.append(hdu[0].header['DATE-OBS'][:10])
    
print(list(set(dates)))