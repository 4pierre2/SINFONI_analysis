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
from get_mag import get_coord_sinfoni, get_2MASS_mags
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
import scipy.signal
from difflib import SequenceMatcher
import glob
import re
# J band

def get_obj(filename):
    hdu = fits.open(filename)
    obj = hdu[0].data
    hdu.close()
    return np.nan_to_num(obj)

def get_coord(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    RA = header['RA']
    dec = header['DEC']
    hdu.close()
    return RA, dec

def get_header_name(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    name = header['OBJECT']
    hdu.close()
    return name


def get_header_pixscale(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    pixscale = header['CDELT2']
    hdu.close()
    return pixscale
    

def get_info(filename):
    hdu = fits.open(filename)
    prog_id = hdu[0].header['HIERARCH ESO OBS PROG ID'].replace('.','').replace('-','').replace('(','').replace(')','')
    filt = hdu[0].header['HIERARCH ESO INS GRAT1 NAME'].replace(' ','')
    pix_scale = hdu[0].header['HIERARCH ESO INS OPTI1 NAME'].replace(' ','')
    date = hdu[0].header['MJD-OBS']
    hdu.close()
    return prog_id, date, filt, pix_scale

def get_obj_info(filename):
    hdu = fits.open(filename)
    prog_id = hdu[0].header['HIERARCH ESO OBS PROG ID'].replace('.','').replace('-','').replace('(','').replace(')','')
    obj_name = hdu[0].header['HIERARCH ESO OBS TARG NAME']
    date = hdu[0].header['MJD-OBS']
    hdu.close()
    return obj_name, date

def get_coadd_name(di):
    possible_filenames = glob.glob(di+'/*COADD*.fits')
    filenames = []
    for filename in possible_filenames:
        if not (('MED' in filename) or ('MASK' in filename)):
            filenames.append(filename)
    if len(filenames) == 1:
        filename = filenames[0]
    return filename


def get_obsobj_name(di):
    possible_filenames = glob.glob(di+'/*OBS_OBJ*.fits')
    filenames = []
    for filename in possible_filenames:
        if not (('MED' in filename) or ('MASK' in filename)):
            filenames.append(filename)
    filename = filenames[0]
    return filename


def find_std_stars(filename, rep='/media/pierre/Disque_2/SNR/archive'):
    prog_id, date, filt, pix_scale = get_info(filename)
    rep_of_interest = rep+'/'+prog_id+'/reflex_end_products'
    for rep in glob.glob(rep_of_interest+'/*'):
        for re in glob.glob(rep+'/*'):
            star_spectra = glob.glob(re+'/*STD_STAR_SPECTRA.fits')
            for spec in star_spectra:
                std_id, date, std_filt, std_pix = get_info(spec)
                if std_id == prog_id and std_filt == filt and std_pix == pix_scale:
                    print(star_spectra, get_obj_info(spec), get_info(spec))
                    
def find_sky(filename, rep='/media/pierre/Disque_2/SNR/archive'):
    prog_id, date, filt, pix_scale = get_info(filename)
    rep_of_interest = rep+'/'+prog_id+'/reflex_end_products'
    for rep in glob.glob(rep_of_interest+'/*'):
        for re in glob.glob(rep+'/*'):
            sky_spectra = glob.glob(re+'/*_SKY_MED.fits')
            for spec in sky_spectra:
                std_id, date, std_filt, std_pix = get_info(spec)
                if std_id == prog_id and std_filt == filt and std_pix == pix_scale:
                    print(sky_spectra, get_obj_info(spec), get_info(spec))
                    
def find_filenames(names, rep='/media/pierre/Disque_2/SNR/archive'):
    # prog_id, date, filt, pix_scale = get_info(filename)
    rep_of_interest = rep#+'/'+prog_id+'/reflex_end_products'
    for rep1 in glob.glob(rep_of_interest+'/*'):
        for rep in glob.glob(rep1+'/*_end*'):
            for rez in glob.glob(rep+'/*'):
                for r in glob.glob(rez+'/*'):
                    files = glob.glob(r+'/*.fits')
                    for file in files:
                        k = True
                        for name in names:
                            if not (name in file):
                                k = False
                        if k:
                            print(file)

def find_dirs(names, rep='/media/pierre/Disque_2/SNR/archive'):
    # prog_id, date, filt, pix_scale = get_info(filename)
    rep_of_interest = rep#+'/'+prog_id+'/reflex_end_products'
    dirs = []
    for rep1 in glob.glob(rep_of_interest+'/*'):
        for rep in glob.glob(rep1+'/*_end*'):
            for rez in glob.glob(rep+'/*'):
                for r in glob.glob(rez+'/*'):
                    files = glob.glob(r+'/*.fits')
                    for file in files:
                        k = True
                        for name in names:
                            if not (name in file):
                                k = False
                        if k:
                            dirs.append(r)
    return set(dirs)

def load_th_spec(filename):
    full = np.loadtxt(filename)
    return full[:,0], full[:,1], full[:,2]

def get_mag_type(std_dir, band='h'):
    spec_star_filenames = glob.glob(std_dir+'/*COADD_*.fits')
    filenames = []
    for filename in spec_star_filenames:
        if not (('MED_' in filename) or ('MASK_' in filename)):
            filenames.append(filename)
    if len(filenames) == 1:
        spec_star_filename = filenames[0]
    ra, dec = get_coord(spec_star_filename)
    ra, dec = get_coord_path(std_dir)
    coords_obj = SkyCoord(ra, dec, unit='deg')
    customSimbad = Simbad()
    customSimbad.add_votable_fields('sptype')
    customSimbad.add_votable_fields('flux('+str(band.upper())+')')
    result_table = customSimbad.query_region(SkyCoord(ra, dec, unit='deg'))
    dists = []
    for r in result_table:
        rastd = r['RA'].split(' ')
        rastd = rastd[0]+'h'+rastd[1]+'m'+rastd[2]+'s'
        decstd = r['DEC'].split(' ')
        decstd = decstd[0]+'d'+decstd[1]+'m'+decstd[2]+'s'
        coords_std = SkyCoord(rastd, decstd)
        dist = coords_obj.separation(coords_std)
        dists.append(dist.value)
    
    x = np.argmin(result_table["FLUX_"+band.upper()])
    x = np.argmin(dists)
    mag = result_table["FLUX_"+band.upper()][x]
    sptype = re.sub('[\W_]+', '', result_table["SP_TYPE"][x].decode().lower())
    return mag, sptype

def query_simbad(ra, dec, radius=10/3600):
    sptype = None
    try:
        customSimbad = Simbad()
        customSimbad.add_votable_fields('sptype')
        customSimbad.add_votable_fields('flux(H)')
        # customSimbad.add_votable_fields('maintype')
        criteria ='ra > '+str(ra-radius)+' & '+'ra < '+str(ra+radius)+' & dec > '+str(dec-radius)+' & '+'dec < '+str(dec+radius)+"& maintypes='Star'"
        # result_table = customSimbad.query_region(SkyCoord(ra, dec, unit='deg'), radius=radius)
        result_table = customSimbad.query_criteria(criteria)
        x = np.argmin(result_table["FLUX_H"])
        mag = result_table["FLUX_H"][x]
        sptype = re.sub('[\W_]+', '', result_table["SP_TYPE"][x].decode().lower())
    except Exception as e:
        a = 1
    return sptype

def get_simbad_name(ra, dec, radius='0d0m20s'):
    sptype = None
    try:
        customSimbad = Simbad()
        customSimbad.add_votable_fields('sptype')
        customSimbad.add_votable_fields('flux(V)')
        # customSimbad.add_votable_fields('maintype')
        # criteria ='ra > '+str(ra-radius)+' & '+'ra < '+str(ra+radius)+' & dec > '+str(dec-radius)+' & '+'dec < '+str(dec+radius)
        result_table = customSimbad.query_region(SkyCoord(ra, dec, unit='deg'), radius=radius)
        # result_table = customSimbad.query_criteria(criteria)
        x = np.argmin(result_table["FLUX_V"])
        mag = result_table["FLUX_V"][x]
        sptype = re.sub('[\W_]+', '', result_table["SP_TYPE"][x].decode().lower())
        name = result_table["MAIN_ID"][x].decode()
    except Exception as e:
        print(e)
        name = ' '
    return name


def get_type(filename):
    ra, dec = get_coord_sinfoni(filename)
    customSimbad = Simbad()
    customSimbad.add_votable_fields('sptype')
    customSimbad.add_votable_fields('flux(H)')
    result_table = customSimbad.query_region(SkyCoord(ra, dec, unit='deg'), radius='0d1m00s')
    x = np.argmin(result_table["FLUX_H"])
    mag = result_table["FLUX_H"][x]
    sptype = re.sub('[\W_]+', '', result_table["SP_TYPE"][x].decode().lower())
    return sptype
 
def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()
   
def find_std_th_spec(sptype, path='/home/pierre/Documents/2021/SINFONI_analysis/pickles_atlas/'):
    names = glob.glob(path+'/*')
    xmax=0
    for name in names:
        x = similar('uk'+sptype, name)
        if x>xmax:
            xmax = x
            name_max = name
    return name_max

def get_fov(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    pscale = header['CDELT2']*3600
    fov = pscale*np.shape(hdu[0].data)[1]
    hdu.close()
    return fov

def get_coord_path(path):
    ra = 0
    dec = 0
    if os.path.isdir(path):
        files = os.listdir(path)
        found = False
        k = 0
        while not (found or k==len(files)):
            file = files[k]
            # print(file, file[-8])
            if "OBJ.fits" == file[-8:]:
                ra, dec = get_coord(path+'/'+file)
                found=True
            elif "STD.fits" == file[-8:]:
                ra, dec = get_coord(path+'/'+file)
                found=True
            elif "PSF.fits" == file[-8:]:
                ra, dec = get_coord(path+'/'+file)
                found=True
            k+=1
    return ra, dec


def get_coadd_name_path(path):
    fi = None
    if os.path.isdir(path):
        files = os.listdir(path)
        found = False
        k = 0
        while not (found or k==len(files)):
            file = files[k]
            # print(file, file[-8])
            if "COADD_OBJ.fits" == file[-14:]:
                fi = file
                found=True
            elif "COADD_STD.fits" == file[-14:]:
                fi = file
                found=True
            elif "COADD_PSF.fits" == file[-14:]:
                fi = file
                found=True
            k+=1
    return fi


    
def get_type_path(path):
    sptype = None
    if os.path.isdir(path):
        files = os.listdir(path)
        found = False
        k = 0
        while not (found or k==len(files)):
            file = files[k]
            # print(file, file[-8])
            if "OBJ.fits" == file[-8:]:
                ra, dec = get_coord(path+'/'+file)
                sptype = query_simbad(ra, dec)
                if sptype != None:
                    found=True
            elif "STD.fits" == file[-8:]:
                ra, dec = get_coord(path+'/'+file)
                sptype = query_simbad(ra, dec)
                if sptype != None:
                    found=True
            elif "PSF.fits" == file[-8:]:
                ra, dec = get_coord(path+'/'+file)
                sptype = query_simbad(ra, dec)
                if sptype != None:
                    found=True
            k+=1
    return sptype



def find_std_dir(obj_dir, rep_ori='/media/pierre/Disque_2/SNR/archive'):
    coadd_names = glob.glob(obj_dir+'/*COADD_OBJ*')
    for coadd in coadd_names:
        if not (('MED_' in coadd) or ('MASK_' in coadd)):
            coadd_name = coadd
    ra_obj, dec_obj = get_coord(coadd_name)
    prog_id, date, filt, pix_scale = get_info(coadd_name)
    fov = get_fov(coadd_name)
    rep_of_interest = rep_ori+'/*/reflex_end_products'
    possible_dirs = []
    radecs = []
    dists = []
    for rep in glob.glob(rep_of_interest+'/*'):
        for rez in glob.glob(rep+'/*'):
            ra, dec = get_coord_path(rez)
            radecs.append([ra,dec])
            dist = ((ra_obj-ra)**2+(dec_obj-dec)**2)**0.5
            if dist*3600 > 5*fov/2**0.5:
                dists.append(dist)
                possible_dirs.append(rez)
    
    sorted_dirs = [x for _, x in sorted(zip(dists, possible_dirs))]
    sorted_dists = sorted(dists)
    found = False
    for delta_t in [0.5,3,15,100]:
        k = 0
        while not (found or k>=len(sorted_dirs)):
            dist = sorted_dists[k]
            di = sorted_dirs[k]
            coadd_name_std = get_coadd_name_path(di)
            if coadd_name_std != None:
                sptype = get_type_path(di)
                prog_id_std, date_std, filt_std, pix_scale_std = get_info(di+'/'+coadd_name_std)
                same_mode = (filt_std==filt) and (pix_scale_std==pix_scale) and (abs(date_std-date) < delta_t)
                if sptype != None and same_mode:
                    found = True
            k+=1
            

    return di

def find_std_dir_old(obj_dir, rep_ori='/media/pierre/Disque_2/SNR/archive'):
    coadd_names = glob.glob(obj_dir+'/*COADD_OBJ*')
    for coadd in coadd_names:
        if not (('MED_' in coadd) or ('MASK_' in coadd)):
            coadd_name = coadd
    ra_obj, dec_obj = get_coord_sinfoni(coadd_name)
    prog_id, date, filt, pix_scale = get_info(coadd_name)
    fov = get_fov(coadd_name)
    rep_of_interest = rep_ori+'/*/reflex_end_products'
    possible_dirs = []
    for rep in glob.glob(rep_of_interest+'/*'):
        for rez in glob.glob(rep+'/*'):
            star_spectra = glob.glob(rez+'/*.fits')
            for spec in star_spectra:
                if not (('MED_' in star_spectra) or ('MASK_' in star_spectra) or ('SKY_' in star_spectra) or ('RESAMPLED' in star_spectra)or ('FLAT' in star_spectra)):
                    std_id, date, std_filt, std_pix = get_info(spec)
                    if std_id == prog_id and std_filt == filt and std_pix == pix_scale:
                        if not (rez.replace('//','/') == obj_dir):
                            possible_dirs.append(rez)
    possible_dirs = list(set(possible_dirs))
    print(possible_dirs)
    # print(possible_dirs)
    dists = []
    full_names = []
    for di in possible_dirs:
        names = glob.glob(di+'/*.fits')
        for name in names:
            ra_std, dec_std = get_coord_sinfoni(name)
            dist = (abs((ra_obj-ra_std)**2+(dec_obj-dec_std)**2))**0.5
            if dist > (fov/2**0.5):
                dists.append(dist)
            else:
                dists.append(1e6)  
            full_names.append(name)
    return possible_dirs[np.argmin(dists)]
                    
def get_band(filename):
    hdu = fits.open(filename)
    header = hdu[0].header
    wl = header['CRVAL3']
    if wl < 1.3:
        band ='j'
    elif wl < 1.8:
        band='h'
    else:
        band='k'
    return band


def get_std_name(di):
    possible_filenames = glob.glob(di+'/*_STD_STAR_SPECTRA.fits')
    filenames = []
    for filename in possible_filenames:
        if not (('MED' in filename) or ('MASK' in filename)):
            filenames.append(filename)
    if len(filenames) == 1:
        filename = filenames[0]
    return filename

def get_std_spectrum(std_dir):
    try:
        name = get_std_name(std_dir)
        spec_star_wl = fits.open(name)[1].data['wavelength']
        spec_star_spectrum = fits.open(name)[1].data['counts_bkg']
        print('ok')
    except:
        name = get_coadd_name(std_dir)
        hdu = fits.open(name)
        cube = np.nan_to_num(hdu[0].data)
        xmax, ymax = np.unravel_index(np.argmax(np.median(cube,0)), np.shape(np.median(cube,0)))
        cube = cube[:,xmax-5:xmax+6,ymax-5:ymax+6]
        spec_star_spectrum = np.median(cube,(1,2))
        spec_star_wl = (np.arange(len(spec_star_spectrum))-hdu[0].header['CRPIX3'])*hdu[0].header['CDELT3']+hdu[0].header['CRVAL3']
        print('ok, but the other one')
        if np.nan_to_num(np.sum(hdu[0].data)) <= 0:
            name = get_obsobj_name(std_dir)
            hdu = fits.open(name)
            cube = np.nan_to_num(hdu[0].data)
            xmax, ymax = np.unravel_index(np.argmax(np.median(cube,0)), np.shape(np.median(cube,0)))
            cube = cube[:,xmax-5:xmax+6,ymax-5:ymax+6]
            spec_star_spectrum = np.median(cube,(1,2))
            spec_star_wl = (np.arange(len(spec_star_spectrum))-hdu[0].header['CRPIX3'])*hdu[0].header['CDELT3']+hdu[0].header['CRVAL3']
            print('ok, but the other other one')
    return spec_star_wl, spec_star_spectrum

#%%

def calibrate(obj_dir):
    print('      ')
    print('      ')
    print('##################')
    print('# Calibrating '+obj_dir)
    print('##################')
    print('      ')
    print('      ')
    
    std_dir = find_std_dir(obj_dir)
    obj_name = get_coadd_name(obj_dir)
    std_name = get_coadd_name(std_dir)
    std_coadd_name = get_coadd_name(std_dir)
    band = get_band(obj_name)
    
    print('Associating '+obj_dir+' with '+str(std_dir)+' : ')
    print(obj_name, ' : ', get_info(obj_name))
    print(std_name, ' : ', get_info(std_name))
        
    if band == 'j':
        F_0 = 3.129e-13*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
        bwidth = 0.162
        filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/J.txt'
        filt = np.loadtxt(filter_name).T
    
    if band == 'h':
        F_0 = 1.133e-13*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
        bwidth = 0.251
        filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/H.txt'
        filt = np.loadtxt(filter_name).T
        
    if band =='k':
        F_0 = 4.283e-14*1e4# *1e4 pour passer de W/cm2/um à W/m2/um
        bwidth = 0.262
        filter_name = '/home/pierre/Documents/2021/SINFONI_analysis/2MASS_filters/K.txt'
        filt = np.loadtxt(filter_name).T
        
    print('      ')
    print('Getting 2MASS survey images and computing magnitude for the object')
    raobj, decobj = get_coord_path(obj_dir)
    fov = get_fov(obj_name)
    mags_obj = np.array(get_2MASS_mags(raobj, decobj, fov=fov, band=band, verbose=False))
    mag_obj = np.mean(mags_obj[~np.isnan(mags_obj)])
    rastd, decstd = get_coord_path(std_dir)
    print(band.upper()+' : ',mag_obj)
    
    print('      ')
    print('Query to SIMBAD to get the spectral type of the STD:')
    mag_cat_std, sptype = get_mag_type(std_dir, band=band)
    spec_star_wl, spec_star_spectrum = get_std_spectrum(std_dir)
    print(sptype.upper())
    
    print('Loading the Pickles spectrum:')
    print(find_std_th_spec(sptype))
    w, f, e = load_th_spec(find_std_th_spec(sptype))
    spec_th = np.interp(spec_star_wl, w/1e4, f)
    print("Computing Transmission")
    transmi = np.nan_to_num(spec_star_spectrum/spec_th)
    transmi /= np.median(transmi)
    transmi[transmi<5e-2]=5e-2


    
    print("Opening ", obj_name)
    obj = get_obj(obj_name)
    xmax, ymax = np.unravel_index(np.argmax(scipy.signal.medfilt(np.median(np.nan_to_num(obj),0),(5,5))), np.shape(obj[0]))
    print("Photocenter estimated at x, y =",xmax,',', ymax)
    obj_detransmitted = obj/transmi[:, np.newaxis, np.newaxis]
    print("Division by transmission")
    
    filter_interp = np.interp(spec_star_wl, filt[0], filt[1])
    obj_2MASSed = obj_detransmitted*filter_interp[:, np.newaxis, np.newaxis]
    limits = [np.max([0,xmax-28]), np.min([xmax+28, np.shape(obj[0])[0]]),np.max([0,ymax-28]), np.min([ymax+28, np.shape(obj[0])[1]])]
    tot_adu = np.sum(np.nan_to_num(obj_2MASSed[:,limits[0]:limits[1],limits[2]:limits[3]]))
    print(tot_adu, 'ADU')
    tot_flux = F_0*10**(-0.4*mag_obj)*bwidth
    print(tot_flux, ' W.m-2$')
    adu = tot_flux/tot_adu
    print('=>   1 ADU = ',adu,' W.m-2')
    obj_final = obj_detransmitted*adu
    print("Calibrating in flux, from ADU to W.m-2")
    return obj_final


#%%


rep = '/media/pierre/Disque_2/SNR/archive'
dirs = list(find_dirs('COADD', rep))

    #%%
    
import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")
    
    


rep_ori='/media/pierre/Disque_2/SNR/archive'
rep_of_interest = rep_ori+'/*/reflex_end_products'
possible_dirs = []
radecs = []
n1 = len(glob.glob(rep_of_interest+'/*'))
i = 0
j = 0
tot = 0
for rep in glob.glob(rep_of_interest+'/*'):
    n2 = len( glob.glob(rep+'/*'))
    j=0
    for rez in glob.glob(rep+'/*'):
        old_stdout = sys.stdout
        # print(tot, ' : ', i, '/', n1, '    ', j, '/', n2)
        # print(rez)
        names = glob.glob(rez+'/*COADD_OBJ*.fits')
        ok = True
        desired_name = 'NGC1808'
        for n in names:
            name = get_header_name(n)
            if desired_name in name:
                ok = True
                print(desired_name, name)
        print(rez)
        if (len(names) > 0) and ok:
            try:
                if not os.path.isdir(rez+'/calibrated/'):
                    os.mkdir(rez+'/calibrated/')
                log_file = open(rez+'/calibrated/'+"calibration.log","w")
                # sys.stdout = log_file

                obj_calibrated = calibrate(rez)
                
                hdu = fits.open(get_coadd_name(rez))
                
                hdu[0].data = obj_calibrated
                
                hdu.writeto(rez+'/calibrated/calibrated.fits', overwrite=True)
                print("Saving to "+rez+'/calibrated/calibrated.fits', name)
            except Exception as e:
                print(e)
                print('ERREUR')
                print('ERREUR')
                print('ERREUR')
        #     log_file.close()
        # if os.path.isfile(rez+'/calibrated/'+"calibration.log"):
        #     log_file = open(rez+'/calibrated/'+"calibration.log","r")
        #     if ok:
        #         print(log_file.read())
        #     log_file.close()
            
            
        

        sys.stdout = old_stdout
        

        j+=1
        tot+=1
    i+=1
radecs = np.array(radecs)

#%%
import shutil

rep_ori='/media/pierre/Disque_2/SNR/archive'
rep_of_interest = rep_ori+'/*/reflex_end_products'

rep_save = '/media/pierre/Disque_2/SNR/CALIBRATED_1/'
if not os.path.isdir(rep_save):
    os.mkdir(rep_save)

n1 = len(glob.glob(rep_of_interest+'/*'))
i = 0
j = 0
tot = 0
for rep in glob.glob(rep_of_interest+'/*'):
    n2 = len( glob.glob(rep+'/*'))
    j=0
    for rez in glob.glob(rep+'/*'):
        old_stdout = sys.stdout
        print(tot, ' : ', i, '/', n1, '    ', j, '/', n2)
        # print(rez)
        # print(glob.glob(rez))
        # print(glob.glob(rez+'/calibrated/*'))
        names = glob.glob(rez+'/*COADD_OBJ*.fits')
        ok = True
        desired_name = 'NGC1808'
        for n in names:
            name = get_header_name(n)
            if desired_name in name:
                ok = True
                print(desired_name, name)
        if (len(glob.glob(rez+'/*COADD_OBJ*.fits')) > 0) and ok:
            if os.path.isfile(rez+'/calibrated/calibrated.fits') and os.path.isfile(rez+'/calibrated/calibration.log'):
                print(glob.glob(rez+'/calibrated/calibrated.fits'))
                ra, dec = get_coord(rez+'/calibrated/calibrated.fits')
                simbad_name = get_simbad_name(ra, dec)
                header_name = get_header_name(rez+'/calibrated/calibrated.fits')
                print(simbad_name)
                print(header_name)
                di = rep_save+'/'+re.sub('[\W_]+', '', simbad_name)
                band = get_band(rez+'/calibrated/calibrated.fits').upper()+'_'+str(int(3600*1e3*get_header_pixscale(rez+'/calibrated/calibrated.fits')))
                num_path = '0'
                if not os.path.isdir(di):
                    os.mkdir(di)
                if not os.path.isdir(di+'/'+band):
                    os.mkdir(di+'/'+band)
                    os.mkdir(di+'/'+band+'/N_0')
                else:
                    last_num = '0'
                    last_num_int = int(last_num)
                    for f in glob.glob(di+'/'+band+'/N_*'):
                        num = f.replace(di+'/'+band+'/N_','').replace('/','')
                        if int(num)>last_num_int:
                            last_num = num
                            last_num_int = int(num)
                    num_path = str(last_num_int+1)
                    print(last_num_int)
                    print(num_path)
                    os.mkdir(di+'/'+band+'/N_'+num_path)
                    print('')
                # last_num_int = int(last_num)
                shutil.copy(rez+'/calibrated/calibrated.fits', di+'/'+band+'/N_'+num_path+'/calibrated.fits')
                shutil.copy(rez+'/calibrated/calibration.log', di+'/'+band+'/N_'+num_path+'/calibration.log')
                
                if not os.path.isfile(di+'/'+band+'/'+re.sub('[\W_]+', '', header_name)):
                    np.savetxt(di+'/'+band+'/'+re.sub('[\W_]+', '', header_name), [ra, dec])
                
                        
                    
                    
                # log_file = open(rez+'/calibrated/'+"calibration.log","w")
        j+=1
        tot+=1
    i+=1