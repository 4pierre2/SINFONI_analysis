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
from sinfobj import sinfobj, fit_line_gauss, fit_line_cube
from fit_speed import fit_los_velocity_exp_disk, los_velocity_exp_disk
from copy import deepcopy
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from scipy.constants import c, pi, G, parsec
from scipy.optimize import curve_fit
from astropy.constants import M_sun
from matplotlib.backends.backend_pdf import PdfPages
import scipy.signal 
from astropy.cosmology import FlatLambdaCDM
import shutil
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

G2 = G*parsec**3/M_sun.value
def get_coord(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    RA = header['RA']
    dec = header['DEC']
    hdu.close()
    return RA, dec

def get_band(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    wl = header['CRVAL3']
    if wl < 1.3:
        band ='j'
    elif wl < 1.8:
        band='h'
    elif wl < 2.:
        band='hk'
    else:
        band='k'
    return band


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
        try:
            name = result_table["MAIN_ID"][x].decode()
        except:
            name = result_table["MAIN_ID"][x]
        z = result_table["Z_VALUE"][x]
        J = result_table["FLUX_J"][x]
        H = result_table["FLUX_H"][x]
        K = result_table["FLUX_K"][x]
    except Exception as e:
        print(e)
    return name, z, J, H, K

def test_line(lam, spec, rad=5, threshold=3):
    filtered = medfilt(spec,3)
    xmax = np.argmax(filtered)
    not_max = np.concatenate([spec[:xmax-rad], spec[xmax+rad:]])
    if np.max(spec-np.median(spec))>threshold*np.std(not_max-np.median(spec)):
        return True
    else:
        return False
    

def load_star_spectra(path = '/home/pierre/Documents/SINFONI_ANALYSIS/v1/GNIRS_spectra'):
    liste = glob.glob(path+'/*')
    specs = []
    wls = []
    for filename in liste:
        try:
            hdu = fits.open(filename)
            spec = hdu[0].data
            wl = (np.arange(len(hdu[0].data))-hdu[0].header['CRPIX1'])*hdu[0].header['CDELT1']+hdu[0].header['CRVAL1']
            hdu.close()
            specs.append(spec)
            wls.append(wl)
        except:
            print('Unable to load '+filename)
    return wls, specs

        # hdu.close()


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

filenames = np.concatenate([glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/N*/*/N_0/*.fits'),glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/M*/*/*/*.fits')])

filenames = glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/N*/*/*/*.fits')

# filenames = glob.glob('/media/pierre/Disque_2/SNR/CALIBRATED_1/NGC1808*/*/*/*.fits')


x = np.arange(1.48,2.22,1e-3)
objs = []
def automatic_analysis(filename, em_lines = 'default', kernel=31, radius=8, extent=[-1,1,-1,1], half_width=6e-3, mode='gauss'):
    name, z, J, H, K = get_simbad_info(filename)
    luminosity_distance = cosmo.luminosity_distance(z)
    band = get_band(filename).upper()
    hdu = fits.open(filename)
    obj = sinfobj(hdu, z=z)
    obj.data[~np.isfinite(obj.data)]=0
    hdu.close()
    obj.recenter_on_max()
    pc_per_sec = cosmo.kpc_proper_per_arcmin(z).value*1e3/60
    pix_scale = obj.pix_scale*pc_per_sec
    
    filter_names = obj.filter_names
    for filter_name in filter_names:
        filt = np.loadtxt(filter_name).T
        filter_interp = np.interp(obj.lam, filt[0], filt[1])
        
        lam0 = obj.lam[filter_interp>0.3][0]
        lam1 = obj.lam[filter_interp>0.3][-1]
        
        obj.plot_im( filt=filt, log=False)
        plt.title(name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
        obj.plot_im(filt=filt, log=True)
        plt.title(name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
        obj.plot_spec(lam0=lam0, lam1=lam1)
        plt.title(name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
        
        obj.plot_im(ra0=-1, ra1=1, dec0=-1, dec1=1,filt=filt, log=False)
        plt.title('Central 2"x2" of '+name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
        obj.plot_im(ra0=-1, ra1=1, dec0=-1, dec1=1,filt=filt, log=True)
        plt.title('Central 2"x2" of '+name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
        obj.plot_spec(ra0=-1, ra1=1, dec0=-1, dec1=1, lam0=lam0, lam1=lam1)
        plt.title('Central 2"x2" of '+name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
        
        obj.plot_im(ra0=-0.5, ra1=0.5, dec0=-0.5, dec1=0.5,filt=filt, log=False)
        plt.title('Central 1"x1" of '+name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
        obj.plot_im(ra0=-0.5, ra1=0.5, dec0=-0.5, dec1=0.5,filt=filt, log=True)
        plt.title('Central 1"x1" of '+name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
        obj.plot_spec(ra0=-0.5, ra1=0.5, dec0=-0.5, dec1=0.5, lam0=lam0, lam1=lam1)
        plt.title('Central 1"x1" of '+name+' '+filter_name.split('/')[-1].replace('.txt',''))
        
    if em_lines == 'default':
        em_lines = obj.detect_lines(ra0=extent[0], ra1=extent[1], dec0=extent[2], dec1=extent[3], threshold=5, kernel=11)[0]
    print(em_lines)
    k = 0
    for em_line in em_lines:
         
        l, sp_full = obj.get_spec(lam0=em_line-half_width, lam1=em_line+half_width)
        l, sp_extent = obj.get_spec(ra0=extent[0], ra1=extent[1], dec0=extent[2], dec1=extent[3], lam0=em_line-half_width, lam1=em_line+half_width)
 
        l, RA_full, dec_full, data_full = obj.get_small_cube(lam0=em_line-half_width, lam1=em_line+half_width)
        l, RA_ext, dec_ext, data_extent = obj.get_small_cube(lam0=em_line-half_width, lam1=em_line+half_width, ra0=extent[0], ra1=extent[1], dec0=extent[2], dec1=extent[3])
 
        print('before fit')
        cube_p_full, cube_err_full, im_snr_full = fit_line_cube(l, data_full/obj.pix_bw, mode=mode)
        cube_p_extent, cube_err_extent, im_snr_extent = fit_line_cube(l, data_extent/obj.pix_bw, mode=mode)
        print('after_fit')
        ra0 = np.min(obj.RA)
        ra1 = np.max(obj.RA)
        dec0 = np.min(obj.RA)
        dec1 = np.max(obj.RA)
        ex_full = [np.max([ra0,ra1]), np.min([ra0,ra1]), np.min([dec0, dec1]), np.max([dec0,dec1])]
        ex_extent = [np.max([extent[0],extent[1]]), np.min([extent[0],extent[1]]), np.min([extent[2], extent[3]]), np.max([extent[2], extent[3]])]
        
        
        
        
        
        
        
        
        obj.plot_spec(lam0=em_line-half_width, lam1=em_line+half_width)
        plt.title("Integrated emission line at $"+"{:.3f}".format(em_line)+"\ \mu m$")
        
        obj.plot_spec(ra0=extent[0], ra1=extent[1], dec0=extent[2], dec1=extent[3], lam0=em_line-half_width, lam1=em_line+half_width)
        plt.title("Central region emission line at $"+"{:.3f}".format(em_line)+"\ \mu m$")
        
        
        
        
        
        
        
        im_flux_full = np.nan_to_num(cube_p_full[:,:,0]*(np.pi/cube_p_full[:,:,2])**0.5)/obj.pix_scale**2
        medim_flux_full = medfilt(im_flux_full, (3,3))
        verymedim_flux_full = medfilt(im_flux_full, (7,7))
        im_flux_extent = np.nan_to_num(cube_p_extent[:,:, 0]*(np.pi/cube_p_extent[:,:,2])**0.5)/obj.pix_scale**2
        medim_flux_extent = medfilt(im_flux_extent, (3,3))
        verymedim_flux_extent = medfilt(im_flux_extent, (7,7))
        
        plt.figure(); plt.tight_layout()
        plt.title("Full $"+"{:.3f}".format(em_line)+"\ \mu m$")
        vmax = np.max(verymedim_flux_full)
        plt.imshow(im_flux_full, extent=ex_full, vmin = 0, vmax= vmax)
        cbar = plt.colorbar()
        cbar.set_label('Flux ($W.m^{-2}.arcsec^{-2}$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        
        plt.figure(); plt.tight_layout()
        plt.title("Central $"+"{:.3f}".format(em_line)+"\ \mu m$")
        vmax = np.max(verymedim_flux_extent)
        plt.imshow(im_flux_extent, extent=ex_extent, vmin = 0, vmax= vmax)
        cbar = plt.colorbar()
        cbar.set_label('Flux ($W.m^{-2}.arcsec^{-2}$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        
        plt.figure(); plt.tight_layout()
        plt.title("Full $"+"{:.3f}".format(em_line)+"\ \mu m$")
        vmax = np.max(verymedim_flux_full)
        plt.imshow(medim_flux_full, extent=ex_full, vmin = 0, vmax= vmax)
        cbar = plt.colorbar()
        cbar.set_label('Flux ($W.m^{-2}.arcsec^{-2}$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        
        plt.figure(); plt.tight_layout()
        plt.title("Central $"+"{:.3f}".format(em_line)+"\ \mu m$")
        vmax = np.max(verymedim_flux_extent)
        plt.imshow(medim_flux_extent, extent=ex_extent, vmin = 0, vmax= vmax)
        cbar = plt.colorbar()
        cbar.set_label('Flux ($W.m^{-2}.arcsec^{-2}$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        
        if (em_line < 2.14) and (em_line > 2.08):  
            
            factor = 5.0776e16*luminosity_distance.value**2
            special_map_full = factor*im_flux_full/pc_per_sec**2
            special_map_extent = factor*im_flux_extent/pc_per_sec**2
            med_special_map_full = factor*medim_flux_full/pc_per_sec**2
            med_special_map_extent = factor*medim_flux_extent/pc_per_sec**2
            
            ex_full_pc = np.array(ex_full)*pc_per_sec
            ex_extent_pc = np.array(ex_extent)*pc_per_sec
            
            vmax = np.max(medfilt(med_special_map_full,(5,5)))
            
            plt.figure(); plt.tight_layout()
            plt.imshow(special_map_full, cmap='jet', origin='lower', extent=ex_full_pc, vmin = 0, vmax= vmax)
            plt.xlabel('Relative Right Ascension (pc)')
            plt.ylabel('Relative Declination (pc)')
            cbar = plt.colorbar()
            cbar.set_label('Mass ($M_\odot.pc^{-2}$)')
            plt.title('Mass = '+str(np.sum(special_map_full)*pix_scale**2)+r' $M_\odot$')
            
            plt.figure(); plt.tight_layout()
            plt.imshow(special_map_extent, cmap='jet', origin='lower', extent=ex_extent_pc, vmin = 0, vmax= vmax)
            plt.xlabel('Relative Right Ascension (pc)')
            plt.ylabel('Relative Declination (pc)')
            cbar = plt.colorbar()
            cbar.set_label('Mass ($M_\odot.pc^{-2}$)')
            plt.title('Mass = '+str(np.sum(special_map_extent)*pix_scale**2)+r' $M_\odot$')
            
            plt.figure(); plt.tight_layout()
            plt.imshow(med_special_map_full, cmap='jet', origin='lower', extent=ex_full_pc, vmin = 0, vmax= vmax)
            plt.xlabel('Relative Right Ascension (pc)')
            plt.ylabel('Relative Declination (pc)')
            cbar = plt.colorbar()
            cbar.set_label('Mass ($M_\odot.pc^{-2}$)')
            plt.title('Mass = '+str(np.sum(med_special_map_full)*pix_scale**2)+r' $M_\odot$')
            
            plt.figure(); plt.tight_layout()
            plt.imshow(med_special_map_extent, cmap='jet', origin='lower', extent=ex_extent_pc, vmin = 0, vmax= vmax)
            plt.xlabel('Relative Right Ascension (pc)')
            plt.ylabel('Relative Declination (pc)')
            cbar = plt.colorbar()
            cbar.set_label('Mass ($M_\odot.pc^{-2}$)')
            plt.title('Mass = '+str(np.sum(med_special_map_extent)*pix_scale**2)+r' $M_\odot$')
        
        if (em_line < 1.67) and (em_line > 1.64):  
            pir2 = 4*np.pi*(luminosity_distance.value*1e6*parsec)**2
            factor = 1e7*8*pir2/1e35
            
            special_map_full = factor*im_flux_full/pc_per_sec**2
            special_map_extent = factor*im_flux_extent/pc_per_sec**2
            med_special_map_full = factor*medim_flux_full/pc_per_sec**2
            med_special_map_extent = factor*medim_flux_extent/pc_per_sec**2
            
            ex_full_pc = np.array(ex_full)*pc_per_sec
            ex_extent_pc = np.array(ex_extent)*pc_per_sec
            
            vmax = np.max(scipy.signal.medfilt(med_special_map_full, (5,5)))
            
            plt.figure(); plt.tight_layout()
            plt.imshow(special_map_full, cmap='jet', origin='lower', extent=ex_full_pc, vmin = 0, vmax= vmax)
            plt.xlabel('Relative Right Ascension (pc)')
            plt.ylabel('Relative Declination (pc)')
            cbar = plt.colorbar()
            cbar.set_label('Supernova rate ($10^{-7}\ yr^{-1}.pc^{-2}$)')
            plt.title('SNR = '+str(np.sum(special_map_full)*pix_scale**2/1e7)+r' SNR/yr')
            
            plt.figure(); plt.tight_layout()
            plt.imshow(special_map_extent, cmap='jet', origin='lower', extent=ex_extent_pc, vmin = 0, vmax= vmax)
            plt.xlabel('Relative Right Ascension (pc)')
            plt.ylabel('Relative Declination (pc)')
            cbar = plt.colorbar()
            cbar.set_label('Supernova rate ($10^{-7}\ yr^{-1}.pc^{-2}$)')
            plt.title('SNR = '+str(np.sum(special_map_extent)*pix_scale**2/1e7)+r' SNR/yr')
            
            plt.figure(); plt.tight_layout()
            plt.imshow(med_special_map_full, cmap='jet', origin='lower', extent=ex_full_pc, vmin = 0, vmax= vmax)
            plt.xlabel('Relative Right Ascension (pc)')
            plt.ylabel('Relative Declination (pc)')
            cbar = plt.colorbar()
            cbar.set_label('Supernova rate ($10^{-7}\ yr^{-1}.pc^{-2}$)')
            plt.title('SNR = '+str(np.sum(med_special_map_full)*pix_scale**2/1e7)+r' SNR/yr')
            
            plt.figure(); plt.tight_layout()
            plt.imshow(med_special_map_extent, cmap='jet', origin='lower', extent=ex_extent_pc, vmin = 0, vmax= vmax)
            plt.xlabel('Relative Right Ascension (pc)')
            plt.ylabel('Relative Declination (pc)')
            cbar = plt.colorbar()
            cbar.set_label('Supernova rate ($10^{-7}\ yr^{-1}.pc^{-2}$)')
            plt.title('SNR = '+str(np.sum(med_special_map_extent)*pix_scale**2/1e7)+r' SNR/yr')
            
            
        
        im_pos_full = np.nan_to_num(cube_p_full[:,:,1])
        medim_pos_full = medfilt(im_pos_full, (3,3))
        verymedim_pos_full = medfilt(im_pos_full, (7,7))
        veryverymedim_pos_full = medfilt(im_pos_full, (21,21))
        im_pos_extent = np.nan_to_num(cube_p_extent[:,:,1])
        medim_pos_extent = medfilt(im_pos_extent, (3,3))
        verymedim_pos_extent = medfilt(im_pos_extent, (7,7))
        
        
        
        pos_med_full = np.median(medim_pos_full)
        pos_med_extent = np.median(medim_pos_extent)
        pos_verymed_full = np.median(medim_pos_full)
        pos_verymed_extent = np.median(medim_pos_extent)
        
        v_full = c*(im_pos_full-pos_med_full)/pos_med_full
        medv_full = c*(verymedim_pos_full-pos_med_full)/pos_med_full
        v_extent = c*(im_pos_extent-pos_verymed_extent)/pos_verymed_extent
        medv_extent = c*(verymedim_pos_extent-pos_verymed_extent)/pos_verymed_extent
        
        delta_speed = c*(np.max(veryverymedim_pos_full)-pos_med_full)/pos_med_full
        delta_pos = delta_speed/c
        plt.figure(); plt.tight_layout()
        plt.imshow(im_pos_full, extent=ex_full, vmin = pos_med_full-delta_pos, vmax = pos_med_full+delta_pos, cmap = 'jet')
        cbar = plt.colorbar()
        cbar.set_label('Wavelength shift($\mu m$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        plt.title("Central wavelength fit around $"+"{:.3f}".format(pos_med_full)+"\ \mu m$")
        
        
        plt.figure(); plt.tight_layout()
        plt.imshow(im_pos_extent, extent=ex_extent, vmin = pos_med_extent-delta_pos, vmax = pos_med_extent+delta_pos, cmap = 'jet')
        cbar = plt.colorbar()
        cbar.set_label('Wavelength shift ($\mu m$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        plt.title("Central wavelength fit around $"+"{:.3f}".format(pos_med_extent)+"\ \mu m$")
        
        
        plt.figure(); plt.tight_layout()
        plt.imshow(v_full*1e-3, extent=ex_full, vmin = -delta_speed*1e-3, vmax = delta_speed*1e-3, cmap = 'jet')
        cbar = plt.colorbar()
        cbar.set_label('Speed ($km.s^{-1}$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        plt.title("Central wavelength fit around $"+"{:.3f}".format(pos_med_full)+"\ \mu m$")
        
        plt.figure(); plt.tight_layout()
        plt.imshow(medv_full*1e-3, extent=ex_full, vmin = -delta_speed*1e-3, vmax = delta_speed*1e-3, cmap = 'jet')
        cbar = plt.colorbar()
        cbar.set_label('Speed ($km.s^{-1}$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        plt.title("Central wavelength fit around $"+"{:.3f}".format(pos_med_full)+"\ \mu m$")
        
        plt.figure(); plt.tight_layout()
        plt.imshow(v_extent*1e-3, extent=ex_extent, vmin = -delta_speed*1e-3, vmax = delta_speed*1e-3, cmap = 'jet')
        cbar = plt.colorbar()
        cbar.set_label('Speed ($km.s^{-1}$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        plt.title("Central wavelength fit around $"+"{:.3f}".format(pos_med_extent)+"\ \mu m$")
        
        plt.figure(); plt.tight_layout()
        plt.imshow(medv_extent*1e-3, extent=ex_extent, vmin = -delta_speed*1e-3, vmax = delta_speed*1e-3, cmap = 'jet')
        cbar = plt.colorbar()
        cbar.set_label('Speed ($km.s^{-1}$)')
        plt.xlabel('Relative Right Ascension (")')
        plt.ylabel('Relative Declination (")')
        plt.title("Central wavelength fit around $"+"{:.3f}".format(pos_med_extent)+"\ \mu m$")
        
        mask = 12
        v_map = medv_full[mask:-mask, mask:-mask]
        ra = obj.RA[mask:-mask]
        dec = obj.dec[mask:-mask]
        mesh = np.meshgrid(ra, dec)
        
        xy = [mesh[1].flatten()*pc_per_sec*parsec, mesh[0].flatten()*pc_per_sec*parsec]
        v = v_map.flatten()
        
        p0 = [1e3,1*parsec,1*parsec,-np.pi/4,np.pi/4,1e2*parsec,1e2*M_sun.value/parsec**2]
        p_constraint = [0,0, 0, None, np.pi/4, 1e2*parsec,1e2*M_sun.value/parsec**2]
        p1, cov1, full_p1 = fit_los_velocity_exp_disk(xy, v, p0, p_constraint)
        
        
        p_constraint = [0, 0, 0, p1[0], np.pi/4, None, None]
        print('before vel_fit')
        p2, cov2, full_p2 = fit_los_velocity_exp_disk(xy, v, p0, p_constraint)
        print('after_vel_fit')
        plt.figure(); plt.tight_layout()
        plt.subplot(121)
        plt.imshow(v_map*1e-3, extent=ex_full, vmin = -delta_speed*1e-3, vmax = delta_speed*1e-3, cmap = 'jet')
        plt.subplot(122)
        plt.imshow(los_velocity_exp_disk(np.array([mesh[1], mesh[0]])*pc_per_sec*parsec, *full_p2)*1e-3, extent=ex_full, vmin = -delta_speed*1e-3, vmax = delta_speed*1e-3, cmap = 'jet')
        plt.title('$r_{d}='+"{:.0f}".format(full_p2[-2]/parsec)+'\ pc$, $\Sigma_{0}='+"{:.0f}".format(full_p2[-1])+'\ M_{\odot}.pc^{-2}$')
        cbar = plt.colorbar()
        cbar.set_label('Speed ($km.s^{-1}$)')
        k += 1

        # if (em_line < 1.67) and (em_line > 1.64):     

        #     pir2 = 4*np.pi*(luminosity_distance.value*1e6*parsec)**2
        #     snr = 10**(1.01*np.log10(medfilt(fl, (1,1))*pir2*1e7)-41.17)
        #     snr = 8*medfilt(fl, (1,1))*pir2/1e35
        #     plt.figure(); plt.tight_layout()
        #     im = snr/pix_scale**2
        #     vmin = 0
        #     vmax = np.max(np.nan_to_num(medfilt(im,(1,1))))
        #     plt.imshow(im, vmin=vmin, vmax=vmax, origin='lower',  extent=extent_pc)
        #     plt.xlabel('Relative Right Ascension (pc)')
        #     plt.ylabel('Relative Declination (pc)')
        #     plt.title('SNR = '+str(np.sum(snr))+r'SNR/yr')
        #     cbar = plt.colorbar()
        #     cbar.set_label('Supernova rate ($10^{-7}\ yr^{-1}.pc^{-2}$)')
      
        
        
        
        
        # fl, pos = obj.get_gauss_flux_pos_em_line(lam0=em_line-half_width, lam1=em_line+half_width, ra0=-1, ra1=1, dec0=-1, dec1=1, radius=radius, kernel=kernel)

        # plt.figure(); plt.tight_layout()
        # obj.plot_im(ra0=-1, ra1=1, dec0=-1, dec1=1, lam0=em_line-half_width, lam1=em_line+half_width, newfig=False)
        # plt.gca().images[-1].colorbar.remove()
        # plt.imshow(medfilt(fl, (3,3)), vmin = 0, vmax = np.max(medfilt(fl, (1,1))), extent=[extent[0], extent[1], extent[2], extent[3]], origin='lower')
        # plt.colorbar()
        # plt.title("Emission l*ne at $"+"{:.3f}".format(em_line)+"\ \mu m$")
        # plt.figure(); plt.tight_layout()
        # obj.plot_im(ra0=-1, ra1=1, dec0=-1, dec1=1, lam0=em_line-half_width, lam1=em_line+half_width, newfig=False)
        # plt.gca().images[-1].colorbar.remove()
        # vmin = np.min(pos)*1000
        # vmax = np.max(pos)*1000
        # plt.imshow(1000*medfilt(pos, (3,3)), vmin = vmin, vmax=vmax, cmap = 'jet', extent=[extent[0], extent[1], extent[2], extent[3]], origin='lower')
        # plt.colorbar()
        # plt.title("Emission line at $"+"{:.3f}".format(em_line)+"\ \mu m$")
        
        # plt.figure(); plt.tight_layout()
        # obj.plot_im(ra0=-1, ra1=1, dec0=-1, dec1=1, lam0=em_line-half_width, lam1=em_line+half_width, newfig=False)
        # plt.gca().images[-1].colorbar.remove()
        # vmax = np.min([500, 1e-3*c*(np.median(pos)-np.min(pos))/np.median(pos)])
        # vmin = np.max([-500, 1e-3*c*(np.median(pos)-np.max(pos))/np.median(pos)])
        # plt.imshow(1e-3*c*(np.median(pos)-medfilt(pos, (3,3)))/np.median(pos), vmin = vmin, vmax=vmax, cmap = 'jet', extent=[extent[0], extent[1], extent[2], extent[3]], origin='lower')
        # plt.colorbar()
        # plt.title("Emission line at $"+"{:.3f}".format(em_line)+"\ \mu m$")
        
    
        # extent_pc = [extent[0]*pc_per_sec, extent[1]*pc_per_sec, extent[2]*pc_per_sec, extent[3]*pc_per_sec]
        # print(extent)
        # if (em_line < 2.14) and (em_line > 2.08):     
        #     plt.figure(); plt.tight_layout()
        #     map_h2hotmass = 5.0776e16*luminosity_distance.value**2*medfilt(fl, (1,1))/pix_scale**2
        #     plt.imshow(map_h2hotmass, cmap='jet', origin='lower', extent=extent_pc)
        #     plt.xlabel('Relative Right Ascension (pc)')
        #     plt.ylabel('Relative Declination (pc)')
        #     cbar = plt.colorbar()
        #     cbar.set_label('Mass ($M_\odot.pc^{-2}$)')

        # if (em_line < 1.67) and (em_line > 1.64):     

        #     pir2 = 4*np.pi*(luminosity_distance.value*1e6*parsec)**2
        #     snr = 10**(1.01*np.log10(medfilt(fl, (1,1))*pir2*1e7)-41.17)
        #     snr = 8*medfilt(fl, (1,1))*pir2/1e35
        #     plt.figure(); plt.tight_layout()
        #     im = snr/pix_scale**2
        #     vmin = 0
        #     vmax = np.max(np.nan_to_num(medfilt(im,(1,1))))
        #     plt.imshow(im, vmin=vmin, vmax=vmax, origin='lower',  extent=extent_pc)
        #     plt.xlabel('Relative Right Ascension (pc)')
        #     plt.ylabel('Relative Declination (pc)')
        #     plt.title('SNR = '+str(np.sum(snr))+r'SNR/yr')
        #     cbar = plt.colorbar()
        #     cbar.set_label('Supernova rate ($10^{-7}\ yr^{-1}.pc^{-2}$)')
      
        # k+=1
                


    # fl, fl_err, pos, pos_err = obj.fit_stellar_velocity(lam0=2.25, lam1=2.42, ra0=-1, ra1=1, dec0=-1, dec1=1, kernel=kernel)
    # plt.figure(); plt.tight_layout()
    # plt.imshow(medfilt(pos, (3,3)))
    # plt.colorbar()


            
        
    
    # spec = obj.get_spec(-1, 1, -1, 1)[1]
    # conti = obj.get_conti(-1, 1, -1, 1, kernel=101)[1]
    # plt.plot((spec-conti)/np.median(conti))
    # em_line = (spec-conti)/np.median(conti)
    
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
        
    
# em_lines=[]
for filename in filenames:
    print(30*[filename])
    plt.ioff()
    # print(filename)
    try:
        automatic_analysis(filename, mode='gauss')
    except Exception as e:
        print(e)
    rep = filename.replace('/calibrated.fits', '')
    rep_to_plot = rep+'/plots_new_gauss'
    if not os.path.isdir(rep_to_plot):
        os.mkdir(rep_to_plot)
    else:
        shutil.rmtree(rep_to_plot)
        os.mkdir(rep_to_plot)
    plt.figure()
    plt.plot(np.arange(-np.pi, np.pi, 1e-2), np.sin(np.arange(-np.pi, np.pi, 1e-2)))
    multipage(rep_to_plot+'/plots.pdf')
    manypdf(rep_to_plot+'/plot')
    manypng(rep_to_plot+'/plot')
    plt.close('all')
    plt.ion()
    # name, z, J, H, K = get_simbad_info(filename)
    # hdu = fits.open(filename)
    # obj = sinfobj(hdu, z=z)
    # header=hdu[0].header
    # print(get_band(filename), header['CRVAL3'],header['CDELT3'],header['CRPIX3'])
    # hdu.close()
    
    
    # spec = obj.get_spec(0.5, 1, 0.5, 1)[1][220:1800]
    # conti = obj.get_conti(0.5, 1, 0.5, 1, kernel=21)[1][220:1800]
    # lam = obj.lam[220:1800]
    # em_line = (spec-conti)/np.median(conti)
    # em_line_interp = np.interp(x, lam, em_line)
    
    # if np.isfinite(em_line_interp).all():
    #     if test_line(x[153:175], em_line_interp[153:175], rad = 3):
    #         em_lines.append(em_line_interp)
            # plt.figure(); plt.tight_layout()
            # plt.plot(x[153:175], em_line_interp[153:175])
            # plt.pause(0.1)
    # plt.plot((spec-conti)/np.median(conti))
    # plt.ylim(-0.1,1)
    # plt.pause(0.1)

# def A_l(wl, Av, Rv=4.05):
#     return Av/Rv*(2.659*(-1.857+1.040/wl)+Rv)

# def multipage(filename, figs=None, dpi=200):
#     pp = PdfPages(filename)
#     if figs is None:
#         figs = [plt.figure(n) for n in plt.get_fignums()]
#     for fig in figs:
#         fig.savefig(pp, format='pdf')
#     pp.close()
    

# def manypdf(filename, figs=None, dpi=200, forma='pdf'):
#     if figs is None:
#         figs = [plt.figure(n) for n in plt.get_fignums()]
#     k=0
#     for fig in figs:
#         pp = PdfPages(filename+'_'+str(k)+forma)
#         fig.savefig(pp, format=forma)
#         pp.close()
#         k+=1
    
# def manypng(filename, figs=None, dpi=200):
#     forma = 'png'
#     if figs is None:
#         figs = [plt.figure(n) for n in plt.get_fignums()]
#     k=0
#     for fig in figs:
#         fig.savefig(filename+'_'+str(k)+'.'+forma, format=forma)
#         k+=1
        
# hdu_1808_j= fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/J.fits')
# hdu_1808_h = fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/H.fits')
# hdu_1808_k = fits.open('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/K.fits')
# ngc1808_j = sinfobj(hdu_1808_j)
# ngc1808_h = sinfobj(hdu_1808_h)
# ngc1808_k = sinfobj(hdu_1808_k)


# plt.figure(figsize=(12,6))
# ngc1808_j.plot_spec(-1,1,-1,1,newfig=False)
# ngc1808_h.plot_spec(-1,1,-1,1, newfig=False)
# ngc1808_k.plot_spec(-1,1,-1,1, newfig=False)
# plt.axis([1.11,2.45,0,9e-14])



# a18, b18, c18, d18, mask18 = ngc1808_j.fit_all_gauss(1.286,1.285,1.29)
# ab18, bb18, cb18, db18, maskb18= ngc1808_k.fit_all_gauss(2.173,2.170,2.178)
# ext1808 =  11.52*np.log10(5.89/(a18/ab18))


# mask1808 = np.zeros(np.shape(ext1808))

# mask1808[:,:6]=1
# mask1808[:,61:]=1
# # inds = np.array([[59, 17],[58, 17],[51, 43],[50, 43],[43, 55],[42, 55],[41, 56],[40, 56],[17, 59],[16, 59]])
# # for ind in inds:
# #     print(ind[0], ind[1])
# #     mask1808[ind[0], ind[1]]=1
# # inds = np.array([[23, 29], [22, 29], [63,  9],[62,  9],[53, 12],[52, 12],[55, 17],[54, 17],[17, 42],[16, 42],[17, 51],[16, 51],[11, 52],[10, 52],[ 5, 48],[ 4, 48],[ 3, 26],[ 2, 26],[ 3, 27],[ 2, 27],[ 7, 38],[ 6, 38],[ 5, 48],[ 4, 48]])
# # for ind in inds:
# #     print(ind[0], ind[1])
# #     mask1808[ind[0], ind[1]]=1
    
    
    
# map1808 = np.ma.masked_array(ext1808, mask1808)

# map1808 = scipy.signal.medfilt(map1808,(3,3))
# map1808 = np.ma.masked_array(map1808, mask1808)
# # plt.imshow(map1808, vmin=0, vmax=5, cmap='jet')
# plt.figure(); plt.tight_layout()
# plt.imshow(map1808, origin='lower', vmin=0, vmax=6, cmap='jet')
# plt.xlabel('Relative Right Ascension (")')
# plt.ylabel('Relative Declination (")')
# cbar = plt.colorbar()
# cbar.set_label('$A_v$')

# # plt.colorbar()



# a_ab = 5.89/10**(7/11.52)

 
# ah218, bh218, ch218, dh218, maskh218 = ngc1808_k.fit_all_gauss(2.130,2.120,2.132)

# plt.figure(); plt.tight_layout()
# map_h2flux = np.ma.masked_array(10**(0.4*A_l(2.12,ext1808))*ah218*(np.pi/ch218)**0.5, mask1808)
# plt.imshow(map_h2flux, cmap='jet', vmin = 8e-21, vmax=1e-19, origin='lower', norm=LogNorm())
# plt.xlabel('Relative Right Ascension (")')
# plt.ylabel('Relative Declination (")')
# cbar = plt.colorbar()
# cbar.set_label('Flux ($W.m^{-2}.pixel^{-1}$)')


# plt.figure(); plt.tight_layout()
# map_h2hotmass = 5.0776e16*Dl**2*map_h2flux/pix_scale**2
# plt.imshow(map_h2hotmass, cmap='jet', origin='lower')
# plt.xlabel('Relative Right Ascension (")')
# plt.ylabel('Relative Declination (")')
# cbar = plt.colorbar()
# cbar.set_label('Mass ($M_\odot.arcsec^{-1}$)')

# # plt.figure(); plt.tight_layout()
# map_h2lam = np.ma.masked_array(bh218, mask1808)
# # img = plt.imshow(map_h2lam, cmap='jet', vmin = 2.1282, vmax=2.1302, origin='lower')
# map_h2speed = c*(map_h2lam/z-2.122)/2.122
# plt.figure(); plt.tight_layout()
# plt.imshow(map_h2speed*1e-3, cmap='jet', origin='lower', vmin=-200, vmax=200, extent=[np.min(ngc1808_k.RA), np.max(ngc1808_k.RA), np.min(ngc1808_k.dec), np.max(ngc1808_k.dec)])
# plt.xlabel('Relative Right Ascension (")')
# plt.ylabel('Relative Declination (")')
# cbar = plt.colorbar()
# cbar.set_label('Velocity ($km.s^{-1}$)')
# # as18, bs18, cs18, ds18, masks18 = ngc1808_k.fit_all_gauss_abs(2.330,2.325,2.345)

# # map_slam = np.ma.masked_array(bs18, mask1808)
# # map_sspeed = c*(map_slam/z-2.3227)/2.3227
# # plt.figure(); plt.tight_layout()
# # plt.imshow(map_sspeed*1e-3, cmap='jet', origin='lower', vmin=-200, vmax=200)
# # plt.colorbar()


# def plummer_full(mesh, Vs, M, A, i0, psi0, x0, y0):
#     if i0>2*np.pi:
#         return 1e25*i0*mesh[0]
#     if i0<0*np.pi:
#         return 1e25*i0*mesh[0]
#     if psi0>2*np.pi:
#         return 1e25*i0*mesh[0]
#     if psi0<0*np.pi:
#         return 1e25*i0*mesh[0]
#     xs = mesh[0]
#     ys = mesh[1]
#     R = ((xs-x0)**2+(ys-y0)**2)**0.5
#     psi = np.arctan2(ys, xs)
#     psi[(xs==0)*(ys==0)] = 0
#     T1 = Vs
#     T2 = ((R**2*G2*M)/(R**2+A**2)**1.5)**0.5
#     T3 = (np.sin(i0)*np.cos(psi-psi0))/((np.cos(psi-psi0)**2+np.sin(psi-psi0)**2/np.cos(i0)**2)**0.75)
#     resulting_map = np.array(T1+T2*T3, dtype='float')
#     return resulting_map.flatten()#*0+psi.flatten()


# def plummer(mesh, Vs, M, A, i0, psi0):
#     if i0>2*np.pi:
#         return 1e25*i0*mesh[0]
#     if i0<0*np.pi:
#         return 1e25*i0*mesh[0]
#     if psi0>2*np.pi:
#         return 1e25*i0*mesh[0]
#     if psi0<0*np.pi:
#         return 1e25*i0*mesh[0]
#     xs = mesh[0]
#     ys = mesh[1]
#     R = ((xs)**2+(ys)**2)**0.5
#     psi = np.arctan2(ys, xs)
#     psi[(xs==0)*(ys==0)] = 0
#     T1 = Vs
#     T2 = ((R**2*G2*M)/(R**2+A**2)**1.5)**0.5
#     T3 = (np.sin(i0)*np.cos(psi-psi0))/((np.cos(psi-psi0)**2+np.sin(psi-psi0)**2/np.cos(i0)**2)**0.75)
#     resulting_map = np.array(T1+T2*T3, dtype='float')
#     return resulting_map.flatten()#*0+psi.flatten()


# mesh = np.meshgrid(ngc1808_k.RA*62, ngc1808_k.dec*62 )
# m = mask1808
# pf, cf = curve_fit(plummer_full, [mesh[0][m==0],mesh[1][m==0]], map_h2speed[m==0], p0 = [0, 1e5, 10, 50*np.pi/180, 130*np.pi/180, 0, 0], maxfev=10000)

# p, c = curve_fit(plummer, [mesh[0][m==0],mesh[1][m==0]], map_h2speed[m==0], p0 = [0, 1e5, 100, 50*np.pi/180, 130*np.pi/180], maxfev=10000)

# plt.figure(); plt.tight_layout()
# plt.imshow(np.reshape(plummer_full(mesh, *pf)*1e-3, np.shape(mesh[0])), cmap='jet', origin='lower', vmin=-200, vmax=200, extent=[np.min(ngc1808_k.RA), np.max(ngc1808_k.RA), np.min(ngc1808_k.dec), np.max(ngc1808_k.dec)])
# plt.xlabel('Relative Right Ascension (")')
# plt.ylabel('Relative Declination (")')
# cbar = plt.colorbar()
# cbar.set_label('Velocity ($km.s^{-1}$)')
# # res_min=1e30
# # for M in 10**np.arange(10):
# #     for A in 5**np.arange(5):
# #         for i0 in np.arange(0,2*np.pi,np.pi/8):
# #             for psi0 in np.arange(0,2*np.pi,np.pi/8):
# #                 rm = plummer([mesh[0][m==0],mesh[1][m==0]], 0, M, A, i0, psi0)
# #                 res = np.sum((rm-map_h2speed[m==0])**2)
# #                 if res < res_min:
# #                     res_min = res
# #                     print(M, A, i0, psi0)
    


# mask1808_h = np.zeros(np.shape(ext1808))

# mask1808_h[:,:6]=1
# mask1808_h[:,61:]=1
# # inds = np.array([[59, 17],[58, 17],[51, 43],[50, 43],[43, 55],[42, 55],[41, 56],[40, 56],[17, 59],[16, 59]])
# # for ind in inds:
# #     print(ind[0], ind[1])
# #     mask1808[ind[0], ind[1]]=1

# afe218, bfe218, cfe218, dfe218, maskfe18 = ngc1808_h.fit_all_gauss(1.65, 1.645, 1.655)

# pir2 = 4*np.pi*(Dl*1e6*parsec)**2

# snr = 8.08*10**(0.4*A_l(1.644,ext1808))*afe218*(np.pi/cfe218)**0.5*pir2*1e-35

# snr = np.ma.masked_array(snr, mask1808_h)
# plt.figure(); plt.tight_layout()
# plt.imshow(1e5*snr, origin='lower', vmin=0, vmax=10, extent=[np.min(ngc1808_h.RA), np.max(ngc1808_h.RA), np.min(ngc1808_h.dec), np.max(ngc1808_h.dec)])
# plt.xlabel('Relative Right Ascension (")')
# plt.ylabel('Relative Declination (")')
# cbar = plt.colorbar()
# cbar.set_label('Supernova rate ($10^{-5}\ yr^{-1}.pixel^{-1}$)')


# def multipage(filename, figs=None, dpi=200):
#     pp = PdfPages(filename)
#     if figs is None:
#         figs = [plt.figure(n) for n in plt.get_fignums()]
#     for fig in figs:
#         fig.savefig(pp, format='pdf')
#     pp.close()
    

# def manypdf(filename, figs=None, dpi=200, forma='pdf'):
#     if figs is None:
#         figs = [plt.figure(n) for n in plt.get_fignums()]
#     k=0
#     for fig in figs:
#         pp = PdfPages(filename+'_'+str(k)+'.'+forma)
#         fig.savefig(pp, format=forma)
#         pp.close()
#         k+=1
    
# def manypng(filename, figs=None, dpi=200):
#     forma = 'png'
#     if figs is None:
#         figs = [plt.figure(n) for n in plt.get_fignums()]
#     k=0
#     for fig in figs:
#         fig.savefig(filename+'_'+str(k)+'.'+forma, format=forma)
#         k+=1
# multipage('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/ANALYSIS/plots.pdf')
# manypdf('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/ANALYSIS/plot')
# manypng('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1808/ANALYSIS/plot')
# plt.close('all')
