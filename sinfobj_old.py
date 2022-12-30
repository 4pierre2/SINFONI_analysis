#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 10:25:26 2021

@author: pierre
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

def gauss(x, a, x0, b, c):
    return a*np.exp(-b*(x-x0)**2)+c

class sinfobj():
    def __init__(self, hdu, hdu_mask=''):
        self.header = hdu[0].header
        self.data = np.nan_to_num(hdu[0].data)
        if hdu_mask:
            self.mask = hdu_mask[0].data
        else:
            self.mask = np.zeros(np.shape(self.data))
        self.lam = (np.arange(len(self.data))-self.header['CRPIX3'])*self.header['CDELT3']+self.header['CRVAL3']
        self.pix_bw = self.lam[1]-self.lam[0]
        # self.lam = np.arange(self.header['CRVAL3'], self.header['CRVAL3']+len(self.data)*self.header['CDELT3'], self.header['CDELT3'])
        self.true_RA = (np.arange(len(self.data[0]))-self.header['CRPIX1'])*self.header['CDELT1']+self.header['CRVAL1']
        self.true_dec = (np.arange(len(self.data[0]))-self.header['CRPIX2'])*self.header['CDELT2']+self.header['CRVAL2']
        self.RA = (np.arange(len(self.data[0]))-self.header['CRPIX1'])*self.header['CDELT1']*3600
        self.dec = (np.arange(len(self.data[0]))-self.header['CRPIX2'])*self.header['CDELT2']*3600
        self.pix_scale = self.dec[1]-self.dec[1]
        
    def get_spec_tot(self, mode='sum'):
        if mode == 'sum':
            spec = np.sum(self.data, (1,2))/self.pix_bw
        if mode == 'mean':
            spec = np.mean(self.data, (1,2))
        if mode == 'median':
            spec = np.median(self.data, (1,2))
        return self.lam, spec
    
    def get_max_lam(self, lam0=0, lam1=1e8):
        lam, spec = self.get_spec_tot(mode='median')
        spec[lam<lam0]=-1e5
        spec[lam>lam1]=-1e5
        max_lam = self.lam[np.argmax(spec)]
        return max_lam
            
    def plot_spec_tot(self, newfig=True, legend=True, mode='sum'):
        if newfig:
            plt.figure()
        if mode == 'sum':
            spec = np.sum(self.data, (1,2))
        if mode == 'mean':
            spec = np.mean(self.data, (1,2))
        if mode == 'median':
            spec = np.median(self.data, (1,2))
        if legend:
            plt.plot(self.lam, spec, label=self.header['OBJECT'])
            plt.legend()
        else:
            plt.plot(self.lam, spec)
        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel('Flux (Arbitrary units)')
        
    def plot_spec_tot_wl(self, l0, l1, newfig=True, legend=True, mode='sum'):
        x0 = np.argmin((self.lam-l0)**2)
        x1 = np.argmin((self.lam-l1)**2)
        if newfig:
            plt.figure()
        if mode == 'sum':
            spec = np.sum(self.data, (1,2))
        if mode == 'mean':
            spec = np.mean(self.data, (1,2))
        if mode == 'median':
            spec = np.median(self.data, (1,2))
        if legend:
            plt.plot(self.lam[x0:x1], spec[x0:x1], label=self.header['OBJECT'])
            plt.legend()
        else:
            plt.plot(self.lam[x0:x1], spec[x0:x1])
        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel('Flux (Arbitrary units)')
        
    def getspec_int(self, ra0=-1, ra1=1, dec0=-1, dec1=1, mode='sum'):
        x0 = np.argmin((self.RA-ra0)**2)
        x1 = np.argmin((self.RA-ra1)**2)
        y0 = np.argmin((self.dec-dec0)**2)
        y1 = np.argmin((self.dec-dec1)**2)
        if x0 != x1 and y0!=y1:
            small_cube = self.data[:,np.min([x0,x1]):np.max([x0,x1]),np.min([y0,y1]):np.max([y0,y1])]
            if mode == 'sum':
                spec = np.sum(small_cube, (1,2))/self.pix_bw
            if mode == 'mean':
                spec = np.mean(small_cube, (1,2))
            if mode == 'median':
                spec = np.median(small_cube, (1,2))
        else:
            spec = self.data[0, x0, y0]
        return self.lam, spec
            
    def plot_spec_int(self, ra0=-1, ra1=1, dec0=-1, dec1=1, newfig=True, legend=True, mode='sum'):
        if newfig:
            plt.figure()
        x0 = np.argmin((self.RA-ra0)**2)
        x1 = np.argmin((self.RA-ra1)**2)
        y0 = np.argmin((self.dec-dec0)**2)
        y1 = np.argmin((self.dec-dec1)**2)
        if x0 != x1 and y0!=y1:
            small_cube = self.data[:,np.min([x0,x1]):np.max([x0,x1]),np.min([y0,y1]):np.max([y0,y1])]
            if mode == 'sum':
                spec = np.sum(small_cube, (1,2))
            if mode == 'mean':
                spec = np.mean(small_cube, (1,2))
            if mode == 'median':
                spec = np.median(small_cube, (1,2))
        else:
            spec = self.data[0, x0, y0]
        if legend:
            plt.plot(self.lam, spec, label=self.header['OBJECT'])
            plt.legend()
        else:
            plt.plot(self.lam, spec)
            
        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel('Flux (Arbitrary units)')
        
        
    def get_spec_pos(self, ra=0, dec=0):
        x = np.argmin((self.RA-ra)**2)
        y = np.argmin((self.dec-dec)**2)
        spec = self.data[:,x,y]
        return self.lam, spec
    
    def plot_spec_pos(self, ra=0, dec=0, newfig=True, legend=True):
        if newfig:
            plt.figure()
        x = np.argmin((self.RA-ra)**2)
        y = np.argmin((self.dec-dec)**2)
        spec = self.data[:,x,y]
        if legend:
            plt.plot(self.lam, spec, label=self.header['OBJECT'])
            plt.legend()
        else:
            plt.plot(self.lam, spec)
            
        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel('Flux (Arbitrary units)')
            
    def get_im_lam(self, lam):
        z = np.argmin((self.lam-lam)**2)
        im = self.data[z,:,:]
        return im

    def plot_im_lam(self, lam='default', newfig=True, log=False, title='object'):
        if newfig:
            plt.figure()
            
        if lam == 'default':
            lam = np.mean(self.lam)
        z = np.argmin((self.lam-lam)**2)
        im = self.data[z,:,:]
        if log:
            plt.imshow(im, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)], norm=LogNorm())
        else:
            plt.imshow(im, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)])
        if title == 'object':
            plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        if title == 'lambda':
            plt.title('{:.2f}'.format(lam*1000)+' nm')
            
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        
    def get_im_int(self, lam0='default', lam1='default', mode='sum'):
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        z0 = np.argmin((self.lam-lam0)**2)
        z1 = np.argmin((self.lam-lam1)**2)
        small_cube = self.data[z0:z1,:,:]
        if mode == 'sum':
            im = np.sum(small_cube, 0)
        if mode == 'mean':
            im = np.mean(small_cube, 0)
        if mode == 'median':
            im = np.median(small_cube, 0)
        return im

    def plot_im_int(self, lam0='default', lam1='default', newfig=True, log=False, mode='sum'):
        if newfig:
            plt.figure()
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        z0 = np.argmin((self.lam-lam0)**2)
        z1 = np.argmin((self.lam-lam1)**2)

        small_cube = self.data[z0:z1,:,:]
        if mode == 'sum':
            im = np.sum(small_cube, 0)
        if mode == 'mean':
            im = np.mean(small_cube, 0)
        if mode == 'median':
            im = np.median(small_cube, 0)
            
        if log:
            plt.imshow(im, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)], norm=LogNorm())
        else:
            plt.imshow(im, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)])

        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
    
    def get_db_im(self, lam0a, lam0b, lam1a, lam1b, mode='mean'):
        im0 = self.get_im_int(lam0a, lam0b, mode=mode)
        im1 = self.get_im_int(lam1a, lam1b, mode=mode)
        return im0-im1
    
    def plot_db_im(self, lam0a, lam0b, lam1a, lam1b, mode='mean', newfig=True, log=False):
        im = self.get_db_im(lam0a, lam0b, lam1a, lam1b, mode=mode)
        if newfig:
            plt.figure()
        if log:
            plt.imshow(im, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)], norm=LogNorm())
        else:
            plt.imshow(im, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)])
        
        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
    
    
    def plot_gif(self, lam0, newfig=True):
        # plt.figure(figsize=(15,5))
        fig, axes = plt.subplots(3, 5, figsize=(15,5))
        ks = np.arange(-7,8)
        l = 0
        for l in range(len(axes.flatten())):
            i, j = np.unravel_index(l, np.shape(axes))
            plt.sca(axes[i, j])
            k = ks[l]
            self.plot_im_lam(lam0-k*float(self.header['CDELT3']), newfig=False, title='lambda')
        plt.tight_layout()
    
    
    def get_lininterp_im(self, lam0, bwidth, mode='mean'):
        lam0a = lam0-bwidth/2
        lam0b = lam0+bwidth/2
        lam1 = lam0a-bwidth
        lam2 = lam0b+bwidth
        
        im0 = self.get_im_int(lam0a, lam0b, mode=mode)
        im1 = self.get_im_int(lam1, lam0a, mode=mode)
        im2 = self.get_im_int(lam0b, lam2, mode=mode)
        return im0-(im1+im2)/2, (im1+im2)/2
    
    def plot_em_line(self,lam0, bwidth, mode='mean', newfig=True, log=False):
        im_line, im_conti = self.get_lininterp_im(lam0, bwidth, mode=mode)
        plt.figure(figsize=(14, 7))
        plt.subplot(121)
        if log:
            plt.imshow(im_conti, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)], norm=LogNorm())
        else:
            plt.imshow(im_conti, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)])
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        plt.title('Continuum emission')
        plt.subplot(122)
        if log:
            plt.imshow(im_line, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)], norm=LogNorm())
        else:
            plt.imshow(im_line, extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)])
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        plt.title('Emission line at '+'{:.1f}'.format(lam0*1000)+' nm')
        
    def get_interp_ifu(self, lam0, bwidth, mode='mean'):
        lam0a = lam0-bwidth/2
        lam0b = lam0+bwidth/2
        lam1 = lam0a-bwidth
        lam2 = lam0b+bwidth
        
        z0a = np.argmin((self.lam-lam0a)**2)
        z0b = np.argmin((self.lam-lam0b)**2)
        
        im1 = self.get_im_int(lam1, lam0a, mode=mode)
        im2 = self.get_im_int(lam0b, lam2, mode=mode)
        
        interp_ifu = self.data[z0a:z0b]-(im1+im2)/2
        
        return interp_ifu, (im1+im2)/2
        
    
    def imshow_speed(self, lam0, bwidth, mode='mean', threshold='default', lam_ref='default', newfig=True):
        ifu, conti = self.get_interp_ifu(lam0, bwidth, mode='mean')
        if newfig:
            plt.figure()
        if threshold=='default':
            threshold = np.mean(ifu)
        mean_ifu = np.mean(ifu, 0)
        mask = mean_ifu < threshold
        speed = np.array(np.argmax(ifu, 0), dtype='float')
        imforshow = np.ma.masked_array(speed, mask)
        imforshow -= imforshow[np.unravel_index(np.argmax(mean_ifu), np.shape(mean_ifu))]
        imforshow *= 3e5*float(self.header['CDELT3'])/lam0
        plt.imshow(imforshow, cmap='rainbow', extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)], vmin=-100, vmax=100)
        cbar = plt.colorbar()
        cbar.set_label('Doppler shift ($km.s^{-1}$)', rotation=90)
        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        
    def imshow_speed_fit(self, lam0, bwidth, mode='mean', threshold='default', lam_ref='default'):
        ifu, conti = self.get_interp_ifu(lam0, bwidth, mode='mean')
        if threshold=='default':
            threshold = 1.5*np.mean(ifu)
        mean_ifu = np.mean(ifu, 0)
        mask = mean_ifu < threshold
        speed = np.array(np.argmax(ifu, 0), dtype='float')
        imforshow = np.ma.masked_array(speed, mask)
        imforshow -= imforshow[np.unravel_index(np.argmax(mean_ifu), np.shape(mean_ifu))]
        imforshow *= 3e5*float(self.header['CDELT3'])/lam0
        plt.imshow(imforshow, cmap='rainbow', extent=[np.min(self.RA), np.max(self.RA), np.min(self.dec), np.max(self.dec)])
        cbar = plt.colorbar()
        cbar.set_label('Doppler shift ($km.s^{-1}$)', rotation=90)
        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
    
    def fit_all_gauss(self, lam0, bwidth = 10):
        for i in range(np.shape(self.data)[0]):
            for j in range(np.shape(self.data)[0]):
                p0 = []
        
    