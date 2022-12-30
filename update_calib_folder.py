#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 11:18:19 2021

@author: pierre
"""
import glob
import os
from shutil import copyfile
rep1 = '/home/pierre/Downloads/076B0098A/'
rep1 = '/home/pierre/Downloads/099B0551A/'
rep2 = '/home/pierre/Downloads/install/share/esopipes/datastatic/sinfo-3.3.1/'

f1 = glob.glob(rep1+'/*')
f2 = glob.glob(rep2+'/*')

k=0
for f in f1:
    filename = f.split('/')[-1]
    isthere=False
    for fb in f2:
        filenameb = fb.split('/')[-1]
        if filename == filenameb:
            os.remove(fb)
            k+=1
            isthere=True
    # if not isthere and ".fits" in f:
    #     print(f, f.replace('sinfo-3.3.0', 'sinfo-3.3.1'))
    #     copyfile(f, f.replace('sinfo-3.3.0', 'sinfo-3.3.1'))
    #     k+=1
            
        
    print(k)
        
    