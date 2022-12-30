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
from scipy.optimize import curve_fit
from scipy.signal import medfilt
from scipy.ndimage import median_filter, maximum_filter
from astropy.io import fits
import glob
import scipy.ndimage
from scipy.constants import G, parsec, c
from scipy.special import iv, kn
from scipy.signal import find_peaks



AL_AK = np.loadtxt('/home/pierre/Documents/CLOUDY/AL_AK.csv', delimiter=';')
AL_AK[:,0] = 1/AL_AK[::-1,0]
AL_AK[:,1] = AL_AK[::-1,1]

def extinct(lam, AK):
    alak = np.interp(lam, AL_AK[:, 0], AL_AK[:, 1])
    ext = 10**(0.4*alak*AK)
    return ext


def gauss(x, a, x0, b, c):
    return a*np.exp(-b*(x-x0)**2)+c

def gauss_abs(x, a, x0, b, c):
    if a>0:
        return x*1e5
    else:
        return a*np.exp(-b*(x-x0)**2)+c
    
def flux_gauss(p):
    flux_gauss = p[0]*(np.pi/p[2])**0.5
    return flux_gauss
    
def err_gauss(p, err):
    err_gauss = ((err[0]*(np.pi/p[2])**0.5)**2+(err[2]*p[0]*np.pi**0.5*p[2]**(-3/2)/2)**2)**0.5
    return err_gauss
    
def amplitude(flux, amplitude):
    return np.pi*amplitude**2/flux**2

def erreur_amplitude(flux, err_flux, a, err_a):
    return ((err_a*a*2*np.pi/flux**2)**2+(err_flux*2*np.pi*a**2/flux**3)**2)**0.5  
    
def fit_line_raw(l, sp, kernel = 'default'):
    if kernel =='default':
        kernel = 15
    sp_med = medfilt(sp, kernel)
    pix_scale = l[1]-l[0]
    c = np.min(sp_med)
    a = np.max(sp)-c
    x0 = l[np.argmax(sp)]
    flux = (np.sum(sp)-len(l)*c)*pix_scale
    b = amplitude(flux, a)
    p = np.array([a, x0, b, c])
    
    var = np.var(sp-sp_med)
    # var = np.var(sp-gauss(l,*p))
    err_c = var**0.5
    err_a = (2*var)**0.5
    err_x0 = pix_scale/2
    err_flux = np.sqrt(2*len(l)*var)*pix_scale
    err_b = erreur_amplitude(flux, err_flux, a, err_a)
    err = np.array([err_a, err_x0, err_b, err_c])
    snr = flux/err_flux
    return p, err, snr

def fit_line_gauss(l, sp, maxfev=1e4):
    p0, err_p0, snr_temp = fit_line_raw(l,sp)
    try:
        p, c = curve_fit(gauss, l, sp, p0=p0, maxfev=int(maxfev))
        err = np.sqrt(np.diag(c))
        flux = flux_gauss(p)
        err_flux = err_gauss(p, err)
        snr = flux/err_flux
        if np.isnan(np.sum(err_flux)):
            argh = 1/0
    except Exception as e:
        p = p0
        err = err_p0
        snr = -1
    return p, err, snr

def fit_line_cube(l, cube, mode='raw'):
    if mode == 'raw':
        f = fit_line_raw
    elif mode == 'gauss':
        f = fit_line_gauss
    sh = np.shape(cube)
    cube_p = []
    cube_err = []
    im_snr = []
    for i in range(sh[1]):
        cube_p.append([])
        cube_err.append([])
        im_snr.append([])
        for j in range(sh[2]):
            p, err, snr = f(l, cube[:,i,j])
            cube_p[i].append(p)
            cube_err[i].append(err)
            im_snr[i].append(snr)
    cube_p = np.array(cube_p)
    cube_err = np.array(cube_err)
    im_snr = np.array(im_snr)
    return cube_p, cube_err, im_snr


class sinfobj():
    def __init__(self, hdu, z=0., hdu_mask=''):
        self.header = hdu[0].header
        self.data = np.nan_to_num(hdu[0].data)#/float(self.header['HIERARCH ESO DET DIT'])
        if hdu_mask:
            self.mask = hdu_mask[0].data
        else:
            self.mask = np.zeros(np.shape(self.data))
        self.pix_bw = self.header['CDELT3']/(1.+z)
        # self.lam = (np.arange(len(self.data))-self.header['CRPIX3'])*self.header['CDELT3']+self.header['CRVAL3']
        # self.lam = np.arange(self.header['CRVAL3'], self.header['CRVAL3']+len(self.data)*self.header['CDELT3'], self.header['CDELT3'])
        self.lam = (np.arange(len(self.data))-self.header['CRPIX3'])*self.header['CDELT3']+self.header['CRVAL3']
        self.lam /= 1+z
        # self.true_RA = (np.arange(len(self.data[0]))-self.header['CRPIX1'])*self.header['CDELT1']+self.header['CRVAL1']
        # self.true_dec = (np.arange(len(self.data[0]))-self.header['CRPIX2'])*self.header['CDELT2']+self.header['CRVAL2']
        self.RA = (np.arange(len(self.data[0][0]))-self.header['CRPIX1'])*self.header['CDELT1']*3600
        self.dec = (np.arange(len(self.data[0]))-self.header['CRPIX2'])*self.header['CDELT2']*3600
        self.pix_scale = self.dec[1]-self.dec[0]
        
        wl = self.header['CRVAL3']
        if wl < 1.3:
            self.band ='J'
        elif wl < 1.8:
            self.band='H'
        elif wl < 2.:
            self.band='HK'
        else:
            self.band='K'
            
        if (self.band=='J') or (self.band=='H') or (self.band=='K'):
            self.filter_names=['/home/pierre/Documents/SINFONI_ANALYSIS/v1/2MASS_filters/'+self.band+'.txt']
        
        if self.band=='HK':
            self.filter_names=['/home/pierre/Documents/SINFONI_ANALYSIS/v1/2MASS_filters/H.txt', '/home/pierre/Documents/SINFONI_ANALYSIS/v1/2MASS_filters/K.txt']
    
    def recenter_coordinates(self, RA, dec):
        self.RA = self.RA-RA
        self.dec = self.dec-dec
        
    def recenter_on_max(self):
        im = np.median(self.data,0)
        medim = medfilt(im, (3,3))
        xdec, xRA = np.unravel_index(np.argmax(medim), np.shape(medim))
        RA = self.RA[xRA]
        dec = self.dec[xdec]
        self.recenter_coordinates(RA, dec)
        
    def crop_pix(self, plam0, plam1, pdec0, pdec1, pra0, pra1):
        self.data = self.data[plam0:plam1+1, np.min([pdec0,pdec1]):np.max([pdec0,pdec1])+1, np.min([pra0,pra1]):np.max([pra0,pra1])+1]
        self.lam = self.lam[plam0:plam1+1]
        self.RA = self.RA[np.min([pra0,pra1]):np.max([pra0,pra1])+1]
        # self.true_RA = self.true_RA[np.min([pra0,pra1]):np.max([pra0,pra1])+1]
        self.dec = self.dec[np.min([pdec0,pdec1]):np.max([pdec0,pdec1])+1]
        # self.true_dec = self.true_dec[np.min([pdec0,pdec1]):np.max([pdec0,pdec1])+1]
        
    def crop(self, ra0='default', ra1='default', dec0='default', dec1='default', lam0='default', lam1='default'):
        if ra0 == 'default':
            ra0 = np.min(self.RA)
        if ra1 == 'default':
            ra1 = np.max(self.RA)
        if dec0 == 'default':
            dec0 = np.min(self.dec)
        if dec1 == 'default':
            dec1 = np.max(self.dec)
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        x0, y0 = self.ptx(ra0, dec0)
        x1, y1 = self.ptx(ra1, dec1)
        z0 = self.ltx(lam0)
        z1 = self.ltx(lam1)
        self.crop_pix(z0, z1, y0, y1, x0, x1)
        self.header['CRPIX1'] -= x0
        self.header['NAXIS1'] -= len(self.RA)
        self.header['CRPIX2'] -= y0
        self.header['NAXIS1'] -= len(self.dec)
        self.header['CRPIX3'] -= z0
        self.header['NAXIS1'] -= len(self.lam)

    
    
    def ltx(self, lam):
        return np.argmin((self.lam-lam)**2)
        
    def ptx(self, ra, dec):
        return np.argmin((self.RA-ra)**2), np.argmin((self.dec-dec)**2)
    
    
    def get_small_cube(self, ra0='default', ra1='default', dec0='default', dec1='default', lam0='default', lam1='default'):
        if ra0 == 'default':
            ra0 = np.min(self.RA)
        if ra1 == 'default':
            ra1 = np.max(self.RA)
        if dec0 == 'default':
            dec0 = np.min(self.dec)
        if dec1 == 'default':
            dec1 = np.max(self.dec)
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        x0, y0 = self.ptx(ra0, dec0)
        x1, y1 = self.ptx(ra1, dec1)
        z0 = self.ltx(lam0)
        z1 = self.ltx(lam1)
        if x1 == x0:
            x1 += 1
        if y1 == y0:
            y1 += 1
        if z0 == z1:
            z1 += 1
        small_cube = self.data[z0:z1+1,np.min([y0,y1]):np.max([y0,y1])+1,np.min([x0,x1]):np.max([x0,x1])+1]
        l = self.lam[z0:z1+1]
        r = self.RA[x0:x1+1]
        d = self.dec[y0:y1+1]
        return l, r, d, small_cube
            
    
    def get_spec(self, ra0='default', ra1='default', dec0='default', dec1='default', lam0='default', lam1='default', filt=None):
        if ra0 == 'default':
            ra0 = np.min(self.RA)
        if ra1 == 'default':
            ra1 = np.max(self.RA)
        if dec0 == 'default':
            dec0 = np.min(self.dec)
        if dec1 == 'default':
            dec1 = np.max(self.dec)
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        x0, y0 = self.ptx(ra0, dec0)
        x1, y1 = self.ptx(ra1, dec1)
        z0 = self.ltx(lam0)
        z1 = self.ltx(lam1)
        if x1 == x0:
            x1 += 1
        if y1 == y0:
            y1 += 1
        if z0 == z1:
            z1 += 1
        small_cube = self.data[z0:z1,np.min([y0,y1]):np.max([y0,y1]),np.min([x0,x1]):np.max([x0,x1])]
        
        if not (filt is None):
            filter_interp = np.interp(self.lam[z0:z1], filt[0], filt[1])
            small_cube = small_cube*filter_interp[:, np.newaxis, np.newaxis]
        
        lam = self.lam[z0:z1]
        spec = np.sum(small_cube, (1,2))/self.pix_bw
        return lam, spec
            
    
    def get_conti(self, ra0='default', ra1='default', dec0='default', dec1='default', lam0='default', lam1='default', filt=None, kernel=21):
        lam, spec = self.get_spec(ra0, ra1, dec0, dec1, filt=filt)
        conti = medfilt(spec,kernel)
        
        
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        z0 = self.ltx(lam0)
        z1 = self.ltx(lam1)
        
        return self.lam[z0:z1], conti[z0:z1]
        
    def get_im(self, lam0='default', lam1='default', ra0='default', ra1='default', dec0='default', dec1='default', filt=None):
        if ra0 == 'default':
            ra0 = np.min(self.RA)
        if ra1 == 'default':
            ra1 = np.max(self.RA)
        if dec0 == 'default':
            dec0 = np.min(self.dec)
        if dec1 == 'default':
            dec1 = np.max(self.dec)
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        x0, y0 = self.ptx(ra0, dec0)
        x1, y1 = self.ptx(ra1, dec1)
        z0 = self.ltx(lam0)
        z1 = self.ltx(lam1)
        if x1 == x0:
            x1 += 1
        if y1 == y0:
            y1 += 1
        if z0 == z1:
            z1 += 1
        small_cube = self.data[z0:z1,np.min([y0,y1]):np.max([y0,y1]),np.min([x0,x1]):np.max([x0,x1])]
        
        if not (filt is None):
            filter_interp = np.interp(self.lam[z0:z1], filt[0], filt[1])
            small_cube = small_cube*filter_interp[:, np.newaxis, np.newaxis]
        im = np.sum(small_cube, 0)/self.pix_scale**2
        return im        

            
    def plot_spec(self, ra0='default', ra1='default', dec0='default', dec1='default', lam0='default', lam1='default', title=False, newfig=True, legend=False):
        if newfig:
            plt.figure()
        lam, spec = self.get_spec(ra0, ra1, dec0, dec1, lam0, lam1)
        if legend:
            plt.plot(lam, spec, label=self.header['OBJECT'])
            plt.legend()
        else:
            plt.plot(lam, spec)
        if title:
            plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel('Flux ($W.m^{-2}.\mu m^{-1}$)')
        

    def plot_im(self, lam0='default', lam1='default', ra0='default', ra1='default', dec0='default', dec1='default', newfig=True, log=False, mode='sum', origin='lower', filt=None):
        if ra0 == 'default':
            ra0 = np.min(self.RA)
        if ra1 == 'default':
            ra1 = np.max(self.RA)
        if dec0 == 'default':
            dec0 = np.min(self.dec)
        if dec1 == 'default':
            dec1 = np.max(self.dec)
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
            
        if newfig:
            plt.figure()
        im = self.get_im(lam0, lam1, ra0=ra0, ra1=ra1, dec0=dec0, dec1=dec1, filt=filt)
        
        if log:
            plt.imshow(im, extent=[np.max([ra0,ra1]), np.min([ra0,ra1]), np.min([dec0, dec1]), np.max([dec0,dec1])], norm=LogNorm(), origin=origin)
        else:
            plt.imshow(im, extent=[np.max([ra0,ra1]), np.min([ra0,ra1]), np.min([dec0, dec1]), np.max([dec0,dec1])], origin=origin)

        plt.title(self.header['OBJECT']+'_'+self.header["HIERARCH ESO INS GRAT1 NAME"])
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        cbar = plt.colorbar()
        cbar.set_label('Flux ($W.m^{-2}.arcsec^{-2}$)')
    
    def test_line(self, lam, half_width=6e-3, ra0='default', ra1='default', dec0='default', dec1='default', lam0='default', lam1='default', threshold = 3):

        if lam0 == 'default':
            lam0 = lam-half_width
        if lam1 == 'default':
            lam1 = lam+half_width
                
        l, spec = self.get_spec(ra0, ra1, dec0, dec1, lam0, lam1)
        p_raw, err_raw, snr_raw = fit_line_raw(l, spec)
        p_gauss, err_gauss, snr_gauss = fit_line_gauss(l, spec)
        
        condition = snr_gauss > threshold
        return condition
    
    def detect_lines(self, half_width=6e-3, step=1e-3, kernel=5, ra0='default', ra1='default', dec0='default', dec1='default', lam0='default', lam1='default', threshold = 3):
        
        lines_final = []
        snr_finals = []
        for filter_name in self.filter_names:
            if lam0 == 'default':
                filt = np.loadtxt(filter_name).T
                filter_interp = np.interp(self.lam, filt[0], filt[1])
                lam_0 = self.lam[filter_interp>0.3][0]
            if lam1 == 'default':
                filt = np.loadtxt(filter_name).T
                filter_interp = np.interp(self.lam, filt[0], filt[1])
                lam_1 = self.lam[filter_interp>0.3][-1]
                
            lams = np.arange(lam_0,lam_1,step)
            lines = np.zeros(np.shape(lams))
            k = 0
            for lam in lams:
                lines[k] = self.test_line(lam, half_width=half_width, ra0=ra0, ra1=ra1,dec0=dec0, dec1=dec1, lam0=lam-half_width, lam1=lam+half_width, threshold=threshold)
                k += 1
                
            detection_map = medfilt(lines,kernel)
            lines_temp = lams[find_peaks(detection_map)[0]]
            for line in lines_temp:
                l, spec = self.get_spec(ra0, ra1, dec0, dec1, line-half_width, line+half_width)
                p_gauss, err_gauss, snr_gauss = fit_line_gauss(l, spec)
                if snr_gauss > threshold:
                    lines_final.append(p_gauss[1])
        return lines_final, snr_finals
            
        
    def get_barycenter(self, l, sp, std=0):
        bc = np.sum(l*sp)/np.sum(sp)
        err = (np.sum(std**2*(l-bc)**2)/(np.sum(sp)**2))**0.5
        return bc, err
    
    def get_flux_pos_em_line(self, lam0='default', lam1='default', ra0='default', ra1='default', dec0='default', dec1='default', radius=8, kernel=51):
        if ra0 == 'default':
            ra0 = np.min(self.RA)
        if ra1 == 'default':
            ra1 = np.max(self.RA)
        if dec0 == 'default':
            dec0 = np.min(self.dec)
        if dec1 == 'default':
            dec1 = np.max(self.dec)
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        x0, y0 = self.ptx(ra0, dec0)
        x1, y1 = self.ptx(ra1, dec1)
        z0 = self.ltx(lam0)
        z1 = self.ltx(lam1)
        if x1 == x0:
            x1 += 1
        if y1 == y0:
            y1 += 1
        if z0 == z1:
            z1 += 1
        data = self.data[:,np.min([x0,x1]):np.max([x0,x1]),np.min([y0,y1]):np.max([y0,y1])]
        fls = np.zeros(np.shape(data[0]))
        fls_err = np.zeros(np.shape(data[0]))
        poss = np.zeros(np.shape(data[0]))
        poss_err = np.zeros(np.shape(data[0]))
        
        l = self.lam[z0:z1]
        for i in range(np.shape(data[0])[0]):
            for j in range(np.shape(data[0])[1]):
                spec = data[z0:z1,i,j]
                limits = [np.max([0,z0-kernel]), np.min([z1+kernel, np.shape(data)[0]])]
                cont = medfilt(data[limits[0]:limits[1],i,j], kernel)[kernel:kernel+z1-z0]
                xm = np.argmax(spec-cont)
                em = (spec-cont)-np.mean(np.concatenate([(spec-cont)[:xm-radius], (spec-cont)[xm+radius:]]))
                std = np.std(np.concatenate([em[:xm-radius], em[xm+radius:]]))
                sp = spec-cont
                fl = np.sum(em)
                std = np.std(np.concatenate([em[:xm-radius], em[xm+radius:]]))
                fl_err = std*len(em)**0.5
                xm = np.argmax(sp)
                std = np.std(sp[xm-radius:xm+radius])
                pos, pos_err = self.get_barycenter(l[xm-radius:xm+radius], em[xm-radius:xm+radius], std)
                fls[i,j] = fl
                fls_err[i,j] = fl_err
                poss[i,j] = pos
                poss_err[i,j] = pos_err
        
        
        return fls, fls_err, poss, poss_err
        
    def get_gauss_flux_pos_em_line(self, lam0='default', lam1='default', ra0='default', ra1='default', dec0='default', dec1='default', radius=8, kernel=51):
        if ra0 == 'default':
            ra0 = np.min(self.RA)
        if ra1 == 'default':
            ra1 = np.max(self.RA)
        if dec0 == 'default':
            dec0 = np.min(self.dec)
        if dec1 == 'default':
            dec1 = np.max(self.dec)
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        x0, y0 = self.ptx(ra0, dec0)
        x1, y1 = self.ptx(ra1, dec1)
        z0 = self.ltx(lam0)
        z1 = self.ltx(lam1)
        if x1 == x0:
            x1 += 1
        if y1 == y0:
            y1 += 1
        if z0 == z1:
            z1 += 1
        data = self.data[:,np.min([x0,x1]):np.max([x0,x1]),np.min([y0,y1]):np.max([y0,y1])]
        fls = np.zeros(np.shape(data[0]))
        fls_err = np.zeros(np.shape(data[0]))
        poss = np.zeros(np.shape(data[0]))
        poss_err = np.zeros(np.shape(data[0]))
        
        l = self.lam[z0:z1]
        
        a = np.zeros(np.shape(data[0]))
        b = np.zeros(np.shape(data[0]))
        c = np.zeros(np.shape(data[0]))
        d = np.zeros(np.shape(data[0]))
        mask = np.zeros(np.shape(data[0]))
        for i in range(np.shape(data[0])[0]):
            for j in range(np.shape(data[0])[1]):
                spec = data[z0:z1+1,i,j]/self.pix_bw
                # spec = self.data[:, i, j]/self.pix_bw
                st = medfilt(spec,kernel_size=(3))
                p0 = [np.max(st)-np.min(st), self.lam[z0:z1+1][np.argmax(st)], 1e6, np.min(abs(st))]

                try:
                    p, cov = curve_fit(gauss, self.lam[z0:z1+1], spec, p0 = p0, maxfev=int(1e4))        
                    a[i,j]=p[0]
                    b[i,j]=p[1]
                    c[i,j]=p[2]
                    d[i,j]=p[3]
                except Exception as e:
                    print(e, p0)
                    try:
                        p, cov = curve_fit(gauss, self.lam[z0:z1+1], spec, p0 = p0, maxfev=int(1e4))        
                        a[i,j]=p[0]
                        b[i,j]=p[1]
                        c[i,j]=p[2]
                        d[i,j]=p[3]
                    except:
                        a[i,j]=p0[0]
                        b[i,j]=p0[1]
                        c[i,j]=p0[2]
                        d[i,j]=p0[3]
                        mask[i,j] = True
        return a*(np.pi/c)**0.5, b
                    
        
    
    def fit_all_gauss(self, lam0, lam_m='default', lam_p='default'):
        if lam_m == 'default':
            lam_m = lam0-5e-3
        if lam_p == 'default':
            lam_p = lam0+5e-3
        lam = self.lam
        x = self.ltx(lam0)
        xm = self.ltx(lam_m)
        xp = self.ltx(lam_p)
        a = np.zeros(np.shape(self.data[0]))
        b = np.zeros(np.shape(self.data[0]))
        c = np.zeros(np.shape(self.data[0]))
        d = np.zeros(np.shape(self.data[0]))
        mask = np.zeros(np.shape(self.data[0]), dtype='bool')
        pls=[]
        for i in range(np.shape(self.data[0])[0]):
            for j in range(np.shape(self.data[0])[1]):
                spec = self.data[:, i, j]/self.pix_bw
                st = medfilt(spec[xm:xp+1],kernel_size=(3))
                p0 = [np.max(st)-np.min(st), lam[xm:xp+1][np.argmax(st)], 1e6, np.min(abs(st))]

                try:
                    p, cov = curve_fit(gauss, lam[xm:xp+1], spec[xm:xp+1], p0 = p0, maxfev=int(1e4))        
                    a[i,j]=p[0]
                    b[i,j]=p[1]
                    c[i,j]=p[2]
                    d[i,j]=p[3]
                except Exception as e:
                    # print(e, p0, pl)
                    try:
                        p, cov = curve_fit(gauss, lam[xm:xp+1], spec[xm:xp+1], p0 = p0, maxfev=int(1e4))        
                        a[i,j]=p[0]
                        b[i,j]=p[1]
                        c[i,j]=p[2]
                        d[i,j]=p[3]
                    except:
                        a[i,j]=p0[0]
                        b[i,j]=p0[1]
                        c[i,j]=p0[2]
                        d[i,j]=p0[3]
                        mask[i,j] = True
        return a, b, c, d, mask
                    
                        
    def fit_all_gauss_abs(self, lam0, lam_m='default', lam_p='default'):
        if lam_m == 'default':
            lam_m = lam0-5e-3
        if lam_p == 'default':
            lam_p = lam0+5e-3
        lam = self.lam
        x = self.ltx(lam0)
        xm = self.ltx(lam_m)
        xp = self.ltx(lam_p)
        a = np.zeros(np.shape(self.data[0]))
        b = np.zeros(np.shape(self.data[0]))
        c = np.zeros(np.shape(self.data[0]))
        d = np.zeros(np.shape(self.data[0]))
        mask = np.zeros(np.shape(self.data[0]), dtype='bool')
        pls=[]
        for i in range(np.shape(self.data[0])[0]):
            for j in range(np.shape(self.data[0])[1]):
                spec = self.data[:, i, j]/self.pix_bw
                st = medfilt(spec[xm:xp+1],kernel_size=(3))
                p0 = [np.min(st)-np.max(st), lam[xm:xp+1][np.argmin(st)], 1e6, np.max(abs(st))]

                try:
                    p, cov = curve_fit(gauss_abs, lam[xm:xp+1], spec[xm:xp+1], p0 = p0, maxfev=int(1e4))        
                    a[i,j]=p[0]
                    b[i,j]=p[1]
                    c[i,j]=p[2]
                    d[i,j]=p[3]
                except Exception as e:
                    # print(e, p0, pl)
                    try:
                        p, cov = curve_fit(gauss_abs, lam[xm:xp+1], spec[xm:xp+1], p0 = p0, maxfev=int(1e4))        
                        a[i,j]=p[0]
                        b[i,j]=p[1]
                        c[i,j]=p[2]
                        d[i,j]=p[3]
                    except:
                        a[i,j]=p0[0]
                        b[i,j]=p0[1]
                        c[i,j]=p0[2]
                        d[i,j]=p0[3]
                        mask[i,j] = True
        return a, b, c, d, mask
                    
           
    def load_star_spectra(self, path = '/home/pierre/Documents/SINFONI_ANALYSIS/v1/GNIRS_spectra'):
        liste = glob.glob(path+'/*')
        specs = []
        wls = []
        for filename in liste:
            try:
                hdu = fits.open(filename)
                spec = hdu[0].data
                wl = ((np.arange(len(hdu[0].data))-hdu[0].header['CRPIX1'])*hdu[0].header['CDELT1']+hdu[0].header['CRVAL1'])/1e4
                hdu.close()
                specs.append(spec)
                wls.append(wl)
            except:
                print('Unable to load '+filename)
        return wls, specs
            
    def fit_stellar_velocity(self, lam0='default', lam1='default', ra0='default', ra1='default', dec0='default', dec1='default', kernel=201, kernel_max=35):
        if ra0 == 'default':
            ra0 = np.min(self.RA)
        if ra1 == 'default':
            ra1 = np.max(self.RA)
        if dec0 == 'default':
            dec0 = np.min(self.dec)
        if dec1 == 'default':
            dec1 = np.max(self.dec)
        if lam0 == 'default':
            lam0 = np.min(self.lam)
        if lam1 == 'default':
            lam1 = np.max(self.lam)
        x0, y0 = self.ptx(ra0, dec0)
        x1, y1 = self.ptx(ra1, dec1)
        z0 = self.ltx(lam0)
        z1 = self.ltx(lam1)
        if x1 == x0:
            x1 += 1
        if y1 == y0:
            y1 += 1
        if z0 == z1:
            z1 += 1

        wls, specs = self.load_star_spectra()
        
        wl = self.lam[z0:z1]
        spis = []
        for k in range(len(wls)):
            # plt.plot(wls[k], specs[k])
            si = np.interp(wl, wls[k], specs[k])
            spis.append(si)
        sp = np.mean(spis, 0)
        
        
        def make_spec(lam, a, delta_p):
            sp_moved = scipy.ndimage.shift(sp, delta_p, mode='constant', cval=np.median(sp))
            spec = (a+sp_moved)/(1+a)
            return spec
        
        def fit_raw_p(s):
            
            k2min = 1e100
            xmin = 0
            xmin_spis = 0
            for x_spis in np.arange(len(spis)):
                sp_star = spis[x_spis]
                for x in np.arange(-40,-29,1):
                    sp_moved = np.roll(s, x)
                    k2 = np.sum(abs(sp_moved-sp_star))
                    if k2<k2min:
                        k2min = k2
                        xmin = x
                        xmin_spis = x_spis
            return xmin, xmin_spis
        
        
        data = self.data[:,np.min([x0,x1]):np.max([x0,x1]),np.min([y0,y1]):np.max([y0,y1])]
        ages = np.zeros(np.shape(data[0]))
        poss = np.zeros(np.shape(data[0]))
        
        for i in range(np.shape(data[0])[0]):
            for j in range(np.shape(data[0])[1]):
                spec = data[z0:z1,i,j]
                cont = median_filter(data[:,i,j], kernel)[z0:z1]
                cont_max = maximum_filter(median_filter(data[:,i,j], 11), kernel_max)[z0:z1]
                cont *= np.median(cont_max)/np.median(cont)
                mask = [cont<=(0.1*np.median(cont))]
                cont[mask] = spec[mask]
                dat = spec/cont
                pos_raw, age_raw = fit_raw_p(dat)
                poss[i,j] = pos_raw
                ages[i,j] = age_raw
        return poss, ages
    
"""            
                
filenamej = "/home/pierre//Documents/CLEANED_1/NGC1808_J_125/cleaned.fits"     
filenameh = "/home/pierre//Documents/CLEANED_1/NGC1808_H_125/cleaned.fits"     
filenamek = "/home/pierre//Documents/CLEANED_1/NGC1808_K_125/cleaned.fits"     
      
hduj = fits.open(filenamej)
objj = sinfobj(hduj, z=0)          
cubej = objj.data
lamj = objj.lam
dlamj = lamj[1]-lamj[0]
imj = objj.get_im(filt=np.loadtxt(objj.filter_names[0]).T)
magj = -2.5*np.log10(imj/3.129e-13/0.162)
hduh = fits.open(filenameh)
objh = sinfobj(hduh, z=0)          
cubeh = objh.data
lamh = objh.lam
dlamh = lamh[1]-lamh[0]
imh = objh.get_im(filt=np.loadtxt(objh.filter_names[0]).T)
magh = -2.5*np.log10(imh/1.133e-13/0.251)
hduk = fits.open(filenamek)
objk = sinfobj(hduk, z=0)          
cubek = objk.data
lamk = objk.lam
dlamk = lamk[1]-lamk[0]
imk = objk.get_im(filt=np.loadtxt(objk.filter_names[0]).T)
magk = -2.5*np.log10(imk/4.283e-14/0.262)


# def fit_conti(l, cube):
    
#%%

ssps = []
ages = []
path = '/home/pierre/Documents/SPISEA/SSPS_1/Z0_1A/'

filenames = glob.glob(path+'/*.txt')


for filename in filenames:
    if 'lam.txt' not in filename:
        age = float(filename.split('/')[-1].split('.tx')[0])
        ages.append(age)
        ssps.append(np.loadtxt(filename))
    else:
        lam_ssp = np.loadtxt(filename)
            
#%%

cube = np.concatenate([objj.data/dlamj, objh.data/dlamh, objk.data/dlamk])
lam = np.concatenate([lamj, lamh, lamk])

# cube = cube[::5]
# lam = lam[::5]
l0j = np.argmin(abs(lam_ssp-np.min(lamj)))
l1j = np.argmin(abs(lam_ssp-np.max(lamj)))
lamj_ssp = lam_ssp[l0j:l1j]
l0h = np.argmin(abs(lam_ssp-np.min(lamh)))
l1h = np.argmin(abs(lam_ssp-np.max(lamh)))
lamh_ssp = lam_ssp[l0h:l1h]
l0k = np.argmin(abs(lam_ssp-np.min(lamk)))
l1k = np.argmin(abs(lam_ssp-np.max(lamk)))
lamk_ssp = lam_ssp[l0k:l1k]

ssps_convolved = []
for ssp in ssps:
    sspj = scipy.ndimage.zoom(ssp[l0j:l1j], len(lamj)/len(lamj_ssp), order=1)
    sspj_convolved = scipy.ndimage.gaussian_filter(sspj, (np.mean(lamj)/2000)/(lamj[1]-lamj[0]))
    ssph = scipy.ndimage.zoom(ssp[l0h:l1h], len(lamh)/len(lamh_ssp), order=1)
    ssph_convolved = scipy.ndimage.gaussian_filter(ssph, (np.mean(lamh)/3000)/(lamh[1]-lamh[0]))
    sspk = scipy.ndimage.zoom(ssp[l0k:l1k], len(lamk)/len(lamk_ssp), order=1)
    sspk_convolved = scipy.ndimage.gaussian_filter(sspk, (np.mean(lamk)/4000)/(lamk[1]-lamk[0]))
    ssps_convolved.append(np.concatenate([sspj_convolved, ssph_convolved, sspk_convolved]))


#%%
from scipy.optimize import curve_fit
from scipy.ndimage import shift
from tqdm import tqdm

plt.ioff()

AL_AK = np.loadtxt('/home/pierre/Documents/CLOUDY/AL_AK.csv', delimiter=';')
AL_AK[:,0] = 1/AL_AK[::-1,0]
AL_AK[:,1] = AL_AK[::-1,1]

def extinct(lam, AK):
    alak = np.interp(lam, AL_AK[:, 0], AL_AK[:, 1])
    ext = 10**(0.4*alak*AK)
    return ext

ages, ssps_convolved = zip(*sorted(zip(ages, ssps_convolved)))


rebin_fac = 5

cube_t = scipy.ndimage.zoom(cube, [1/rebin_fac, 0.25, 0.25], order=1)
ssps_t = scipy.ndimage.zoom(ssps_convolved, [1, 1/rebin_fac], order=1)
lam_t = scipy.ndimage.zoom(lam, 1/rebin_fac)
lam_min = 0.2
x_min = np.argmin(abs(lam_t-lam_min))
lam_max = 2.5
x_max = np.argmin(abs(lam_t-lam_max))
shape = np.shape(cube_t[x_min:x_max])
age_map = np.zeros(shape[1:])
AK_map = np.zeros(shape[1:])
M_map = np.zeros(shape[1:])
roll_map = np.zeros(shape[1:])
model_cube = np.zeros(shape)

for i in tqdm(range(shape[1])):
    for j in tqdm(range(shape[2]), leave=False):
        r2min = 1e100
        spec_t = cube_t[:, i, j][x_min:x_max]
        for N in range(len(ssps_t)):
            ssp_t = ssps_t[N][x_min:x_max]
            age_t = ages[N]
            def f_for_fit(_, a, AK, roll):
                return shift(a*ssp_t*extinct(lam_t, AK)[x_min:x_max], roll, mode='nearest', order=1)
            p, cov = curve_fit(f_for_fit, lam, spec_t, p0 = [np.mean(spec_t)/np.mean(ssp_t), 1, 1], maxfev=100000)
            model_t = f_for_fit(lam_t[x_min:x_max], *p)
            # plt.plot(lam_t[x_min:x_max], model_t, alpha=0.1, c='b')
            # plt.plot(lam_t[x_min:x_max], spec_t, c='r')
            r2 = np.sum((model_t-spec_t)**2)
            if r2 < r2min:
                r2min = r2
                best_model = model_t
                best_age = age_t
                best_roll = p[2]
                best_AK = p[1]
                best_M = p[0]
        roll_map[i, j] = best_roll
        age_map[i, j] = best_age
        AK_map[i, j] = best_AK
        M_map[i, j] = best_M
        model_cube[:, i, j] = best_model
plt.figure()
plt.imshow(age_map)
plt.colorbar()
plt.savefig('/home/pierre/fit0.png')
plt.close()
plt.figure()
plt.imshow(AK_map)
plt.colorbar()
plt.savefig('/home/pierre/fit1.png')
plt.close()
plt.figure()
plt.imshow(M_map)
plt.colorbar()
plt.savefig('/home/pierre/fit2.png')
plt.close()
                
   
    #%%
    
for i in tqdm(range(shape[1])):
    for j in tqdm(range(shape[2]), leave=False):
        plt.plot(lam_t[x_min:x_max], cube_t[:, i, j][x_min:x_max])         
        plt.plot(lam_t[x_min:x_max], model_cube[:, i, j])
        plt.show()         
        plt.pause(0.2)
        plt.clf()
# def prop(spec, a):
#     return a*spec


# def prop(spec, a, shift, sigma):
#     return scipy.ndimage.gaussian_filter(scipy.ndimage.gaussian_filter((a*scipy.ndimage.shift(spec, shift, order=1, mode='nearest')), sigma), 1.48e-3/(lamk[1]-lamk[0]))

# ages, ssps_convolved = zip(*sorted(zip(ages, ssps_convolved)))

# shape = np.shape(cube)
# age1_map = np.zeros(shape[1:])
# age2_map = np.zeros(shape[1:])
# roll_map = np.zeros(shape[1:])
# std_map = np.zeros(shape[1:])
# model_ifu = np.zeros(shape)
# for i in tqdm(range(shape[1])):
#     for j in tqdm(range(shape[2]), leave=False):
#         r2min = 1e100
#         spec = cube[:, i, j]
#         for N1 in np.arange(0, len(ssps_convolved[:200]), 4):
#             for N2 in np.arange(0, len(ssps_convolved[:]), 10):
#                 if (N2 > N1):
#                     ssp1 = ssps_convolved[N1]
#                     age1 = ages[N1]
#                     ssp2 = ssps_convolved[N2]
#                     age2 = ages[N2]
#                     def prop(spec, a, b, shift, sigma):
#                         spect = scipy.ndimage.shift(a*ssp1+b*ssp2, shift, order=1, mode='nearest')
#                         return scipy.ndimage.gaussian_filter(scipy.ndimage.gaussian_filter(spect, sigma), 1.48e-3/(lamk[1]-lamk[0]))

#                     p, c = curve_fit(prop, ssp, spec, p0=[np.median(ssp1)/np.median(ssp)/2, np.median(ssp2)/np.median(ssp)/2, 15, 4], maxfev=100000)
#                     model = prop(ssp, *p)
#                     r2 = np.sum((model-spec)**2)
#                     if r2 < r2min:
#                         best_model = model
#                         best_age1 = age1
#                         best_age2 = age2
#                         best_roll = p[2]
#                         best_std = p[3]
#                         r2min = r2
#         age1_map[i, j] = best_age1
#         age2_map[i, j] = best_age2
#         roll_map[i, j] = best_roll
#         std_map[i, j] = best_std
#         model_ifu[:, i, j] = best_model
#         plt.figure()
#         plt.imshow(age1_map)
#         plt.colorbar()
#         plt.savefig('/home/pierre/fit0.png')
#         plt.close()
#         plt.figure()
#         plt.imshow(age2_map)
#         plt.colorbar()
#         plt.savefig('/home/pierre/fit1.png')
#         plt.close()
#         plt.figure()
#         plt.imshow(roll_map)
#         plt.colorbar()
#         plt.savefig('/home/pierre/fit2.png')
#         plt.close()
#         plt.figure()
#         plt.imshow(std_map)
#         plt.colorbar()
#         plt.savefig('/home/pierre/fit3.png')
#         plt.close()
        
        
        
                    
"""