#!/usr/bin/python

'''
Loic - 07/12/2023
This script is developed to generate ROI masks to be used with FEAT query for ROI analysis
Thus, running this script should happen right after FEAT preprocessing and analysis, and right before Featquery
'''
import os
import subprocess
import glob
#Look at multiple areas
workDir = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/preprocessing.feat'
outDir = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/masksNative'
for a, area in enumerate(['V1_lh','V1_rh', 'V2_lh', 'V2_rh', 'V3_lh', 'V3_rh', 'hV4_lh', 'hV4_rh', 'V3a_lh', 'V3a_rh', 'V3b_lh', 'V3b_rh']): #enumerate(['V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b']):
    area_mask_standard = f'/Users/tonglab/Documents/Loic/atlases/Wang_2015/{area}.nii.gz' #input
    area_mask_native = f'{workDir}/{area}_example_func.nii.gz' #output
    refFunc = f'{workDir}/example_func.nii.gz' #reference
    trans_mat = f'{workDir}/reg/standard2example_func.mat' #transformation matrix
    os.system(f'flirt -in {area_mask_standard} -ref {refFunc} -out {area_mask_native} -init {trans_mat} -interp nearestneighbour -applyxfm')
    os.makedirs(f'{outDir}/improved_masks/{area}_masks', exist_ok=True)
    #Now binarize thresholded zstats, then create mask
    for z in [2, 3]: #[1, 2, 3, 4, 5]
        for th in [2.1, 2.5, 3.1, 3.6]:
            # sorted(glob.glob(os.path.join('/Users/tonglab/Documents/Loic/phantom_mri_data/derivatives/FEAT/events_3column', subject,f'ses-1/task-{scan}', f'*run{run + 1:03}_{condition}.txt')))  # make sure conditions are in same order as the event files stored in event_files

            # List all files in the directory
            # map_name_string = f'zstat{z}_sup_'
            # all_files = os.listdir(directory_path)
            # files_with_specific_string = [file for file in all_files if map_name_string in file]

            zstat = glob.glob(f'{workDir}/stats/improved_zstat_maps/zstat{z}_sup_*.nii.gz')
            thrZstat = f'{workDir}/stats/improved_zstat_maps/thresh{th}_zstat{z}.nii.gz'
            os.makedirs(f'{workDir}/stats/improved_zstat_maps', exist_ok=True)
            ##Create thresholded z-stat map
            os.system(f'fslmaths {zstat[0]} -thr {th} {thrZstat}')

            ## Binarize it
            mask = f'{workDir}/stats/improved_zstat_maps/bin_thresh{th}_zstat{z}.nii.gz'
            os.system(f'fslmaths {thrZstat} -bin {mask}')


            final_mask = f'{outDir}/improved_masks/{area}_masks/zstat{z}_{area}_thr{th}_mask.nii.gz'
            os.system(f'fslmaths {mask} -mul {area_mask_native} {final_mask}')

            # The fslstats command to get the number of non-zero voxels (use '-l 0' to exclude background value)
            command = f'fslstats {final_mask} -V -l 0'

            # Run the command and capture the output
            output = subprocess.check_output(command, shell=True, text=True)

            # The output contains two values: volume and number of voxels
            num_voxels, volume = output.strip().split()

            # Convert the number of voxels to an integer and store it in a Python variable
            num_voxels = int(num_voxels)

            # Print or use the value as needed
            print("Volume:", volume)
            print("Number of Voxels:", num_voxels)

            # Get the original file extension (e.g., '.nii.gz')
            # _, file_extension = os.path.splitext(final_mask)

            # Create the new filename with the number of voxels included
            new_filename = f'{outDir}/improved_masks/{area}_masks/zstat{z}_{area}_thr{th}_numvoxels{num_voxels}.nii.gz'

            # Rename the file with the new filename
            os.rename(final_mask, new_filename)


# 1 Make a mask for the LGN #no activity detected in the LGN
#first, create mask on Fsleyes
#Binarize the thresholded z-stat map
#Multiply both masks together
# outDir = '/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019/data/functional/localizer/run01/outputs/magnitude.feat'
# zstat = f'{outDir}/stats/zstat1.nii.gz'
# thrZstat = f'{outDir}/stats/thresh_zstat1.nii.gz'
# #Create thresholded z-stat map
# os.system(f'fslmaths {zstat} -thr 3.1 {thrZstat}')
#
# ## Create thresholded mask (-bin for binarize)
# #First make a thresholded map (previous step), then binarize it
# mask = f'{outDir}/stats/thresh_zstat1_mask.nii.gz'
# os.system(f'fslmaths {thrZstat} -bin {mask}')

## Multiply both masks together
# lgn_mask =f'{outDir}/stats/zstat1_lgn_mask.nii.gz' #mask created with fsleyes
# final_mask = f'{outDir}/stats/zstat1_lgn_final_mask.nii.gz'
# os.system(f'fslmaths {mask} -mul {lgn_mask} {final_mask}')


## 2 Make masks for V1
#zstat1 = both conditions activate voxel
#zstat2 = inducer patches activate voxel
#zstat 3 = phantom patch activates voxel
#ztat 4 = activation for inducer patches is superior to activation for phantom patch
#ztat 5 = activation for inducer patches is inferior to activation for phantom patch

#Get probabilistic map and transform from standard space to individual space
# outDir = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019/data/functional/localizer/run01/outputs/firstLevel.feat'
#
# v1_mask_standard = '/home/tonglab/Documents/Loic/atlases/Wang_2015/V1.nii.gz' #input
# v1_mask_individual = f'{outDir}/V1_example_func.nii.gz' #output
# refFunc = f'{outDir}/example_func.nii.gz' #reference
# trans_mat = f'{outDir}/reg/standard2example_func.mat' #transformation matrix
# os.system(f'flirt -in {v1_mask_standard} -ref {refFunc} -out {v1_mask_individual} -init {trans_mat} -interp nearestneighbour -applyxfm')
#
# #Now binarize thresholded zstats, then create mask
# for z in [1, 2, 3, 4, 5]:
#     for th in [2.1, 2.5, 3.1, 3.6]:
#         zstat = f'{outDir}/stats/zstat{z}.nii.gz'
#         thrZstat = f'{outDir}/stats/thresh{th}_zstat{z}.nii.gz'
#         ##Create thresholded z-stat map
#         os.system(f'fslmaths {zstat} -thr {th} {thrZstat}')
#
#         ## Binarize it
#         mask = f'{outDir}/stats/bin_thresh{th}_zstat{z}.nii.gz'
#         os.system(f'fslmaths {thrZstat} -bin {mask}')
#
#         # Multiply both masks together
#         v1_mask = f'{outDir}/V1_example_func.nii.gz'
#         final_mask = f'{outDir}/stats/zstat{z}_v1_thr{th}_mask.nii.gz'
#         os.system(f'fslmaths {mask} -mul {v1_mask} {final_mask}')








