#!/usr/bin/env python
#Loic Daumail 04 20 2023
'''
The bold data was collected with less slices (4 less on top, 4 less at the bottom of the brain, along z axis).
Thus, this script's goal is to trim the corresponding top up nifti files to allow topup correction
'''

import os
import os.path as op
import numpy as np
import pandas as pd
import json
import glob
import sys

def resize_topup_runs(niftiPath):
    path = niftiPath
    os.chdir(path)
    TUnii = [f for f in os.listdir(f"{path}") if f.endswith('oppPE.nii')]
    outPath = f"{path}/out"
    os.makedirs(outPath, exist_ok=True)
    #1 resize nifti file (this will also change nifti header dim 3 size)
    for n in TUnii:
        TUPath = glob.glob(f"{path}/{n}")[0]
        hdrdata = os.popen(f"fslinfo {TUPath}").readlines()
        for line in hdrdata:
            if line.startswith('dim1'):
                 xsize = int(line.split()[-1])
            if line.startswith('dim2'):
                 ysize = int(line.split()[-1])
            if line.startswith('dim3'):
                 zsize = int(line.split()[-1])
            if line.startswith('dim4'):
                 tsize = int(line.split()[-1])
        cmd = f"fslroi {TUPath} {outPath}/{n} 0 {xsize} 0 {ysize} 4 {zsize-8} 1 {tsize-1}"
        os.system(cmd)
        os.system(f"rm {TUPath}")
        os.system(f"mv {outPath}/{n}.gz {TUPath}")

    os.system(f"rm -r {outPath}")
