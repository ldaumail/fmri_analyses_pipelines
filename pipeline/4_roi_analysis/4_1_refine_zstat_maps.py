#!/usr/bin/env python
#Loic Daumail 08 07 2023
'''
This script is aimed at refining localizer maps for the inducers and the phantom patches conditions, in a way that only
the z-stats for the inducer higher than the z-stats of the phantom condition are saved in the inducer z-stat map,
and only the z-stats of the phantom patch that are higher than those of the inducer z-stat map are saved in the phantom z-stat map
'''

import numpy as np
import nibabel as nib
import os
workDir = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/preprocessing.feat/stats'
# Replace 'z_stat_map1.nii.gz' and 'z_stat_map2.nii.gz' with the paths to your actual z-stat maps
z_stat2_path = f'{workDir}/zstat2.nii.gz'
z_stat3_path = f'{workDir}/zstat3.nii.gz'

# Load the z-stat maps using nibabel
z_stat2 = nib.load(z_stat2_path).get_fdata()
z_stat3 = nib.load(z_stat3_path).get_fdata()

# Find the highest z-stat value in each map

output_map1 = np.where(z_stat2 > z_stat3, z_stat2, 0)
output_map2 = np.where(z_stat3 > z_stat2, z_stat3, 0)

# Save the output maps using nibabel
output_map1_img = nib.Nifti1Image(output_map1, nib.load(z_stat2_path).affine)
output_map2_img = nib.Nifti1Image(output_map2, nib.load(z_stat3_path).affine)

# Replace 'output_map1.nii.gz' and 'output_map2.nii.gz' with the desired output file paths
os.makedirs(os.path.join(f'{workDir}/improved_zstat_maps'), exist_ok=True)
output_map1_path = f'{workDir}/improved_zstat_maps/zstat2_sup_zstat3.nii.gz'
output_map2_path = f'{workDir}/improved_zstat_maps/zstat3_sup_zstat2.nii.gz'

nib.save(output_map1_img, output_map1_path)
nib.save(output_map2_img, output_map2_path)

print("Output maps generated successfully.")





