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
    if ('H' in sinfo.header["HIERARCH ESO INS GRAT1 NAME"]) and ('K' in sinfo.header["HIERARCH ESO INS GRAT1 NAME"]):
        sinfo.plot_spec_tot(mode='median')
        plt.savefig(rep+'/full_spectrum.png')
        max_lam = sinfo.get_max_lam(lam0 = 1.644, lam1 = 1.700)
    elif 'H' in sinfo.header["HIERARCH ESO INS GRAT1 NAME"]:
        sinfo.plot_spec_tot_wl(1.6, 1.75, mode='median')
        plt.savefig(rep+'/full_spectrum.png')
        max_lam = sinfo.get_max_lam(lam0 = 1.644, lam1 = 1.700)
    else:
        max_lam = sinfo.get_max_lam()
        sinfo.plot_spec_tot(mode='median')
        plt.savefig(rep+'/full_spectrum.png')
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

def makedir(path):
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except:
            print('Error creating : '+path)
            
def make_sorted_dir(obj, grating, pix_scale, suffix='', path = './SORTED_FITS/'):
    makedir(path)
    makedir(path+'/'+obj)
    makedir(path+'/'+obj+'/'+grating)
    makedir(path+'/'+obj+'/'+grating+'/'+pix_scale)
    makedir(path+'/'+obj+'/'+grating+'/'+pix_scale+'/'+suffix)
    return path+'/'+obj+'/'+grating+'/'+pix_scale+'/'+suffix+'/'

rep_source = '/media/pierre/Disque_2/SNR/SORTED_COADD/'
rep_dest = '/media/pierre/Disque_2/SNR/SORTED_PLOTS_COADD/'

reps_source = glob.glob(rep_source+'/*')
for rep in reps_source:
    subreps_grat = glob.glob(rep+'/*')
    for subrep_grat in subreps_grat:
        subreps_pixscale = glob.glob(subrep_grat+'/*')
        for subrep_pixscale in subreps_pixscale:
            files = glob.glob(subrep_pixscale+'/*.fits')
            for file in files:
                if not (('MASK' in file) or ('_MED_' in file)): 
                    hdu = fits.open(file)
                    sinfobject = sinfobj(hdu)
                    header = sinfobject.header
                    mjd = str(header['MJD-OBS']).replace(' ','')
                    obj = header['OBJECT'].replace(' ','')
                    grating = header['HIERARCH ESO INS FILT1 NAME'].replace(' ','')
                    pix_scale = header['HIERARCH ESO INS OPTI1 NAME'].replace(' ','')
                    print(obj, grating, pix_scale, mjd)
                    path = make_sorted_dir(obj, grating, pix_scale, suffix=mjd, path='/media/pierre/Disque_2/SNR/SORTED_PLOTS_COADD')
                    plot_all(sinfobject, path)
    
                    
                
                
