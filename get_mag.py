#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 13:39:27 2021

@author: pierre
"""

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from glob import glob
import urllib
import pandas as pd
import requests
import gzip
import os
import pickle

def get_metadata(ra, dec):
    webUrl = urllib.request.urlopen('https://irsa.ipac.caltech.edu/ibe/search/twomass/allsky/allsky?POS='+str(ra)+','+str(dec))
    data = webUrl.read()
    return data.decode()

def clean_metadata(metadata):
    liste = metadata.split('\n')
    liste_finale = []
    k = 0
    for l in liste:
        if len(l)>0:
            if not l[0]=='\\':
                if not (('double|' in l) or ('deg|' in l) or ('null' in l)):
                    if k == 0:
                        l = l.replace(' ','')
                        l = l.replace('|',',')
                        k+=1
                    else:
                        while '  ' in l:
                            l = l.replace('  ', ' ')
                        l = l.replace(' ', ',')
                    while l[0]==',':
                        l = l[1:]
                    while l[-1]==',':
                        l = l[:-1]
                    liste_finale.append(l.split(',')[:15])
    return liste_finale

def get_url(ordate, hemisphere, scanno, fname):
    ordate = '{:06d}'.format(ordate)
    hemisphere = '{:1s}'.format(hemisphere)
    scanno = '{:03d}'.format(scanno)
    fname = '{:20s}'.format(fname)
    path = ordate+hemisphere+'/s'+scanno+'/image/'+fname
    url = 'https://irsa.ipac.caltech.edu/ibe/data/twomass/allsky/allsky/' + path
    return url

def get_urls(clean_metadata):
    ds = pd.DataFrame(clean_metadata[1:], columns=clean_metadata[0])
    ordates = [int(x) for x in ds['ordate']]
    hemispheres = [str(x) for x in ds['hemisphere']]
    scannos = [int(x) for x in ds['scanno']]
    fnames = [str(x) for x in ds['fname']]
    urls = []
    for k in range(len(hemispheres)):
        urls.append(get_url(ordates[k], hemispheres[k], scannos[k], fnames[k]))
    return urls
            
def download_url(url, filename):   
    r = requests.get(url.replace(' ',''))
    with open(filename, 'wb') as f:
        f.write(gzip.decompress(r.content))

def phot_sinf(im,x0,y0,zp,fov=3):
    y = np.arange(len(im[:,0]))
    x = np.arange(len(im[0,:]))
    xx, yy = np.meshgrid(x,y)
    mask = ((((xx-x0)**2)**0.5)<=fov/2) *((((yy-y0)**2)**0.5)<=fov/2)
    im_no_bg = im-np.median(np.nan_to_num(im))
    phot = np.sum(im_no_bg[mask])
    
    mag = -2.5*np.log10(phot)+zp
    return mag


def get_mag(hdu, ra, dec, fov=3):
    im = hdu[0].data
    w = WCS(hdu[0].header)
    zp = hdu[0].header['MAGZP']
    coords = SkyCoord(ra, dec, unit='deg')
    x, y = w.world_to_pixel(coords)
    pix_scale = 3600*hdu[0].header['CDELT1']
    fov_pix = fov*1./abs(pix_scale)
    mag = phot_sinf(im, x, y, zp, fov=fov_pix)
    return mag

def get_coord_sinfoni(filename):
    hdu = fits.open(filename)
    im = scipy.signal.medfilt(np.median(hdu[0].data, 0), (5,5))
    xmax, ymax = np.unravel_index(np.argmax(im), np.shape(im))
    w = WCS(hdu[0].header)
    wl = hdu[0].header['CRPIX3']
    coords = w.pixel_to_world(xmax, ymax, wl)
    return coords[0].ra.value, coords[0].dec.value
    

def get_2MASS_mags(ra, dec, band='h', path='./2MASS_ims/', fov=3, verbose=True):
    
    rastr = '{:09.5f}'.format(ra)
    decstr = '{:09.5f}'.format(dec)
    
    if not os.path.isdir(path+'/metadatas/'):
        os.mkdir(path+'/metadatas')
    
    if verbose:print("Checking if these RA and DEC values have already been downloaded...")
    if not os.path.isfile(path+'/metadatas/'+rastr+'_'+decstr):
        if verbose:print("No, downloading now")
        metadata_raw = get_metadata(ra, dec)
        pickle.dump(metadata_raw, open(path+'/metadatas/'+rastr+'_'+decstr, 'wb'))
    else:
        if verbose:print("Yes, loading the metadata")
        metadata_raw = pickle.load(open(path+'/metadatas/'+rastr+'_'+decstr, 'rb'))
        
    if verbose:print("Cleaning metadata")
    metadata = clean_metadata(metadata_raw)
    if verbose:print("Building urls")
    urls = get_urls(metadata)
    
    if not os.path.isdir(path):
        os.mkdir(path)
        
    filenames = []
    for url in urls:
        
        if verbose:print("Dealing with "+url.split('/')[-1]+" :")
        filename = path+url.split('/')[-1].split('.gz')[0]
        
        if verbose:print("   Checking wavelength")
        if filename.split('/')[-1][0] == band:
            if verbose:print("   Checking if already downloaded")
            if not os.path.isfile(filename):
                if verbose:print("   Downloading")
                download_url(url, filename)
            filenames.append(filename)
    if verbose:print(' ')
    if verbose:print('All files loaded, calculating mags...') 
    mags = []
    for filename in filenames:
        if verbose:print("   Computing mag for "+filename.split('/')[-1])
        hdu = fits.open(filename)
        mag = get_mag(hdu, ra, dec, fov=fov)
        mags.append(mag)    
        
    if verbose:print('Magnitudes : ', mags)
    return mags
    
# ra, dec = get_coord_sinfoni('/home/pierre/Documents/2021/SINFONI_analysis/SINFONI_CALIB/NGC_1068/K.fits')
# print(get_2MASS_mags(ra,dec, fov=4))
# ra =  148.968458 
# dec =  69.679703 
# print(get_2MASS_mags(ra,dec, band='k', fov=150))


