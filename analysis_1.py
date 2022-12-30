# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from sinfobj import sinfobj

def plot_spec(cube, label = ''):
    plt.plot(np.median(cube, (1,2)), label=label)
    

def plot_all(sinfo, rep):
    try:
        os.mkdir(rep)
    except:
        k = 123
    sinfo.plot_spec_tot_wl(1.6, 1.75, mode='median')
    plt.savefig(rep+'/full_spectrum.png')
    if 'H' in sinfo.header["HIERARCH ESO INS GRAT1 NAME"]:
        max_lam = sinfo.get_max_lam(lam0 = 1.644, lam1 = 1.700)
    else:
        max_lam = sinfo.get_max_lam()
    sinfo.plot_im_lam(max_lam-0.02)
    plt.savefig(rep+'/im_conti.png')
    sinfo.plot_im_lam(max_lam)
    plt.savefig(rep+'/im_em_line.png')
    sinfo.plot_gif(max_lam)
    plt.savefig(rep+'/im_em_line_gif.png')
    sinfo.plot_em_line(max_lam, 0.01)
    plt.savefig(rep+'/im_em_line_intrinsic.png')
    sinfo.imshow_speed(max_lam, 0.01)
    plt.savefig(rep+'/im_em_line_speed.png')
    plt.close('all')
    
    
rep = '/media/pierre/Disque_2/SNR/archive/SINFO_4/reflex_end_products/2021-04-15T23:48:58/'
rep = '/media/pierre/Disque_2/SNR/archive/SINFONI_CenA/reflex_end_products/2021-04-14T20:38:09/'
rep = '/media/pierre/Disque_2/SNR/archive/077B0079A/reflex_end_products/2021-04-22T10:57:05/'
rep = '/media/pierre/Disque_2/SNR/archive/079B0576A/reflex_end_products/2021-04-22T11:41:33/'


meta_rep = glob.glob("/media/pierre/Disque_2/SNR/archive/0*")
for repo in meta_rep:
    print(repo)
    rep = repo+'/reflex_end_products/2021*/'
    
    dirs = glob.glob(rep+'/*')
    
    #%%
    
    hdus = []
    for di in dirs:
        fits_of_interest = glob.glob(di+'/*OBS_OBJ.fits')
        for fit in fits_of_interest:
            hdu = fits.open(fit)
            hdus.append(hdu)
            
    #%%
    
    cubes = {}
    lams = {}
    objs={}
    
    for hdu in hdus:
        header = hdu[0].header
        objs[header['OBJECT']+'_'+header["HIERARCH ESO INS GRAT1 NAME"]] = sinfobj(hdu)
    #%%
    
    #%%
    
    plt.ioff()
    for obj_name in objs:
        o = objs[obj_name]
        rep = '/home/pierre/Documents/2021/SINFONI_analysis/'+o.header['OBJECT']+'_'+o.header["HIERARCH ESO INS GRAT1 NAME"]+'/'
        plot_all(o, rep)
    # objs['NGC253_H'].imshow_speed(1.6445, 0.01)
    # objs['NGC253_H'].plot_gif(1.6445)
