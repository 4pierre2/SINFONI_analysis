#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 08:53:40 2022

@author: pierre
"""

import numpy as np
from astropy.io import fits
from sinfobj import sinfobj

hdu_J = fits.open('/media/pierre/Disque_2/SNR/CLEANED_1/NGC1808_J_125/cleaned.fits')
hdu_H = fits.open('/media/pierre/Disque_2/SNR/CLEANED_1/NGC1808_H_125/cleaned.fits')
hdu_K = fits.open('/media/pierre/Disque_2/SNR/CLEANED_1/NGC1808_K_125/cleaned.fits')

hdus = [hdu_J, hdu_H, hdu_K]

sinfobjs = []

for hdu in hdus:
    sinfobjs.append(sinfobj(hdu))
    