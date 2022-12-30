#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 19:20:06 2021

@author: pierre
"""
from copy import deepcopy
from sinfobj import sinfobj
import numpy as np
import glob
import os
from astropy.io import fits

filenames = glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/*Circi*/K_50/*/*.fits')

def flatten(t):
    return [item for sublist in t for item in sublist]

def shift_2_objs(obj1, obj2):
    
    obj1.recenter_on_max()
    obj2.recenter_on_max()
    new_obj = deepcopy(obj1)
    
    data1 = obj1.data
    RA1 = obj1.RA
    dec1 = obj1.dec
    lam1 = obj1.lam
    pix_scale_1 = np.around(obj1.pix_scale, 3)
    pix_bw_1 = np.around(obj1.pix_bw, 6)
    mesh1 = np.meshgrid(dec1, lam1, RA1)
    
    data2 = obj2.data
    RA2 = obj2.RA
    dec2 = obj2.dec
    lam2 = obj2.lam
    pix_scale_2 = np.around(obj2.pix_scale, 3)
    pix_bw_2 = np.around(obj2.pix_bw, 6)
    mesh2 = np.meshgrid(dec2, lam2, RA2)
    
    if (pix_bw_1 != pix_bw_2) or (pix_scale_1 != pix_scale_2):
        print('Attention, les deux fichiers semblent correspondre à deux modes différents : ')
        print(pix_bw_1, pix_bw_2, pix_scale_1, pix_scale_2)
    
    # print(np.arange(np.min([np.min(obj1.dec), np.min(obj2.dec)]), np.max([np.max(obj1.dec), np.max(obj2.dec)]), pix_scale_1))
    new_RA = np.arange(np.max([np.max(obj1.RA), np.max(obj2.RA)]), np.min([np.min(obj1.RA), np.min(obj2.RA)])-pix_scale_1, -pix_scale_1)
    new_dec = np.arange(np.min([np.min(obj1.dec), np.min(obj2.dec)]), np.max([np.max(obj1.dec), np.max(obj2.dec)])+pix_scale_1, pix_scale_1)
    new_lam = np.arange(np.min([np.min(obj1.lam), np.min(obj2.lam)]), np.max([np.max(obj1.lam), np.max(obj2.lam)])+pix_bw_1, pix_bw_1)
    
    
    x0_RA_1 = np.argmin((new_RA-RA1[0])**2)
    x0_RA_2 = np.argmin((new_RA-RA2[0])**2)
    x0_dec_1 = np.argmin((new_dec-dec1[0])**2)
    x0_dec_2 = np.argmin((new_dec-dec2[0])**2)
    
    
    x0_lam_1 = np.argmin((new_lam-lam1[0])**2)
    x0_lam_2 = np.argmin((new_lam-lam2[0])**2)
    
    new_mesh = np.array(np.meshgrid(new_dec, new_lam, new_RA))
    new_cube_1 = np.zeros(np.shape(new_mesh[0]))
    new_cube_2 = np.zeros(np.shape(new_mesh[0]))
    new_mask_1 = np.ones(np.shape(new_mesh[0]), dtype='boolean')
    new_mask_2 = np.ones(np.shape(new_mesh[0]), dtype='boolean')
    
    new_cube_1[x0_lam_1:x0_lam_1+len(lam1), x0_dec_1:x0_dec_1+len(dec1), x0_RA_1:x0_RA_1+len(RA1)] = data1
    new_cube_2[x0_lam_2:x0_lam_2+len(lam2), x0_dec_2:x0_dec_2+len(dec2), x0_RA_2:x0_RA_2+len(RA2)] = data2
    new_mask_1[x0_lam_1:x0_lam_1+len(lam1), x0_dec_1:x0_dec_1+len(dec1), x0_RA_1:x0_RA_1+len(RA1)] = False 
    new_mask_2[x0_lam_2:x0_lam_2+len(lam2), x0_dec_2:x0_dec_2+len(dec2), x0_RA_2:x0_RA_2+len(RA2)] = False
    
    masked_cube_1 = np.ma.masked_array(new_cube_1, new_mask_1)
    masked_cube_2 = np.ma.masked_array(new_cube_2, new_mask_2)

    return masked_cube_1, masked_cube_2

def shift(objs):
    
    datas = []
    RAs = []
    decs = []
    lams = []
    pix_scales = []
    pix_bws = []
    
    
    for obj in objs:
        obj.recenter_on_max()
        datas.append(obj.data)
        decs.append(obj.dec)
        RAs.append(obj.RA)
        lams.append(obj.lam)
        pix_scales.append(np.around(obj.pix_scale, 6))
        pix_bws.append(np.around(obj.pix_bw, 6))
    
    suspected_problem = False
    suspected_bad_objs = []
    for k in range(len(pix_scales)-1):
        if (pix_bws[0] != pix_bws[k+1]) or (pix_scales[0] != pix_scales[k+1]):
            suspected_problem = True
            suspected_bad_objs.append(k)
            
    if suspected_problem:
        print('Attention, les fichiers suivants semblent correspondre à deux modes différents : ')
        suspected_bad_objs = np.array(suspected_bad_objs)
        print(suspected_bad_objs)
        # print(suspected_bad_objs, pix_bws[0], pix_bws[suspected_bad_objs], pix_scales[0], pix_scales[suspected_bad_objs])
   
    new_RA = np.arange(np.max(flatten(RAs)), np.min(flatten(RAs))-pix_scales[0], -pix_scales[0])
    new_dec = np.arange(np.min(flatten(decs)), np.max(flatten(decs))+pix_scales[0], pix_scales[0])
    new_lam = np.arange(np.min(flatten(lams)), np.max(flatten(lams))+pix_bws[0], pix_bws[0])
    new_mesh = np.array(np.meshgrid(new_dec, new_lam, new_RA))
    
    masked_cubes = []
    
    for k in range(len(pix_scales)):
        x0_RA = np.argmin((new_RA-RAs[k][0])**2)
        x0_dec = np.argmin((new_dec-decs[k][0])**2)
        x0_lam = np.argmin((new_lam-lams[k][0])**2)

        new_cube = np.zeros(np.shape(new_mesh[0]))
        new_mask = np.ones(np.shape(new_mesh[0]), dtype='bool')
        new_cube[x0_lam:x0_lam+len(lams[k]), x0_dec:x0_dec+len(decs[k]), x0_RA:x0_RA+len(RAs[k])] = datas[k]
        new_mask[x0_lam:x0_lam+len(lams[k]), x0_dec:x0_dec+len(decs[k]), x0_RA:x0_RA+len(RAs[k])] = False

        masked_cube = np.ma.masked_array(new_cube, new_mask)
        masked_cubes.append(masked_cube)
        
    return masked_cubes, new_RA, new_dec, new_lam

def make_new_sinfobj(filename_ref, data, RA, dec, lam):
    hdu = fits.open(filename_ref)
    obj = sinfobj(hdu)
    hdu.close()
    
    obj.data = data
    obj.RA = RA
    obj.dec = dec
    obj.lam = lam
    return obj
    
    
def make_new_fits(filename_save, filename_ref, data, header=None):
    hdu = fits.open(filename_ref)
    hdu[0].data = data
    if header != None:
        hdu[0].header = header
    if os.path.isfile(filename_save):
        os.remove(filename_save)
    hdu.writeto(filename_save)
    hdu.close()

# filenames = glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/*1808*/H*/*/*.fits')
# for filename in filenames:
#     if "N_2" in filename:
#         filenames.remove(filename)
# filenames = glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/*Circi*/K_50/*/*.fits')

# objs = []

# for filename in filenames:
#     hdu = fits.open(filename)
#     objs.append(sinfobj(hdu))
#     hdu.close()
    
# ab, RA, dec, lam = shift(objs)

# test = np.median(ab, 0)

# obj = make_new_sinfobj(filenames[1], test, RA, dec, lam)
# make_new_fits("/home/pierre/test.fits", filenames[1], test)




