
'''
Loic - 07/12/2023
This script is developed to generate ROI masks to be used with FEAT query for ROI analysis
Thus, running this script should happen right after FEAT preprocessing and analysis, and right before Featquery
'''
import os
import subprocess
import glob

subject = 'sub-M132'
session = '230725'
#Look at multiple areas
workDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/second_level.gfeat/cope2.feat'
outDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/phantomEnhanceMasksStandard'
outDir_nat = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/phantomEnhanceMasksNative'
os.makedirs(f'{outDir_nat}', exist_ok=True)
# locDirMasks = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/masksNative'
phDirFeat = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/run01/outputs/topUp/noB0/preprocessing.feat'

# localizer_mask_standard = f'{workDir}/cope2.feat/stats/zstat1.nii.gz' #input
# localizer_mask_native = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/preprocessing.feat/reg/phCope2_example_func.nii.gz' #output
# refFunc = f'{workDir}/cope2.feat/example_func.nii.gz' #reference
# trans_mat = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/phantom/run01/outputs/topUp/noB0/preprocessing.feat/reg/standard2example_func.mat' #transformation matrix
# os.system(f'flirt -in {localizer_mask_standard} -ref {refFunc} -out {localizer_mask_native} -init {trans_mat} -interp nearestneighbour -applyxfm')

for a, area in enumerate(['V1_lh','V1_rh', 'V2_lh', 'V2_rh', 'V3_lh', 'V3_rh', 'hV4_lh', 'hV4_rh', 'V3a_lh', 'V3a_rh', 'V3b_lh', 'V3b_rh']): #enumerate(['V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b']):
    area_mask_standard = f'/Users/tonglab/Documents/Loic/atlases/Wang_2015/{area}.nii.gz' #input
    # area_mask_native = f'{workDir}/{area}_example_func.nii.gz' #output
    # refFunc = f'{workDir}/example_func.nii.gz' #reference
    # trans_mat = f'{workDir}/reg/standard2example_func.mat' #transformation matrix
    # os.system(f'flirt -in {area_mask_standard} -ref {refFunc} -out {area_mask_native} -init {trans_mat} -interp nearestneighbour -applyxfm')
    os.makedirs(f'{outDir}/{area}_masks', exist_ok=True)
    os.makedirs(f'{outDir_nat}/{area}_masks', exist_ok=True)

    # for z in [2, 3]: #[1, 2, 3, 4, 5]
    for th in [3.1, 3.6]:
        phantom_zstat = f'{workDir}/stats/zstat1.nii.gz'
        # zstat = glob.glob(f'{workDir}/stats/improved_zstat_maps/zstat{z}_sup_*.nii.gz')
        thrZstat = f'{workDir}/stats/modif_zstat_maps/thresh{th}_zstat1.nii.gz'
        os.makedirs(f'{workDir}/stats/modif_zstat_maps', exist_ok=True)
        # ##Create thresholded z-stat map
        os.system(f'fslmaths {phantom_zstat} -thr {th} {thrZstat}')
        #
        # ## Binarize it
        mask = f'{workDir}/stats/modif_zstat_maps/bin_thresh{th}_zstat1.nii.gz'
        os.system(f'fslmaths {thrZstat} -bin {mask}')
        #
        #
        # mask_native = f'{workDir}/cope2.feat/stats/modif_zstat_maps/native_bin_thresh{th}_zstat1.nii.gz'
        # os.system( f'flirt -in {mask} -ref {refFunc} -out {mask_native} -init {trans_mat} -interp nearestneighbour -applyxfm')
        # localizer_mask = glob.glob(f'{locDir}/improved_masks/{area}_masks/zstat3_{area}_thr3.6_numvoxels*.nii.gz') #zstats3 = phantom region
        final_mask = f'{outDir}/{area}_masks/phzstat2_loczstat3_{area}_thr{th}_mask.nii.gz' #phzstats2 since we use the second contrast of phantom runs (cope2.feat)
        os.system(f'fslmaths {mask} -mul {area_mask_standard} {final_mask}')
        #
        # # The fslstats command to get the number of non-zero voxels (use '-l 0' to exclude background value)
        command = f'fslstats {final_mask} -V -l 0'
        #
        # # Run the command and capture the output
        output = subprocess.check_output(command, shell=True, text=True)
        #
        # # The output contains two values: volume and number of voxels
        num_voxels, volume = output.strip().split()
        #
        # # Convert the number of voxels to an integer and store it in a Python variable
        num_voxels = int(num_voxels)
        #
        # # Print or use the value as needed
        print("Volume:", volume)
        print("Number of Voxels:", num_voxels)
        #
        # # Get the original file extension (e.g., '.nii.gz')
        # # _, file_extension = os.path.splitext(final_mask)

        # Create the new filename with the number of voxels included
        new_filename = f'{outDir}/{area}_masks/phzstat2_loczstat3_{area}_thr{th}_numvoxels{num_voxels}.nii.gz'

        # Rename the file with the new filename
        os.rename(final_mask, new_filename)

        ## Transform mask to native space
        area_mask_native = f'{outDir_nat}/{area}_masks/phzstat2_loczstat3_{area}_thr{th}_numvoxels{num_voxels}.nii.gz' #output
        refFunc = f'{phDirFeat}/example_func.nii.gz' #reference
        trans_mat = f'{phDirFeat}/reg/standard2example_func.mat' #transformation matrix
        os.system(f'flirt -in {new_filename} -ref {refFunc} -out {area_mask_native} -init {trans_mat} -interp nearestneighbour -applyxfm')

        #recount number of voxels
        command = f'fslstats {area_mask_native} -V -l 0'
        #
        # # Run the command and capture the output
        output = subprocess.check_output(command, shell=True, text=True)
        #
        # # The output contains two values: volume and number of voxels
        num_voxels, volume = output.strip().split()
        #
        # # Convert the number of voxels to an integer and store it in a Python variable
        num_voxels = int(num_voxels)
        #
        # # Print or use the value as needed
        print("Volume:", volume)
        print("Number of Voxels:", num_voxels)
        #
        # # Get the original file extension (e.g., '.nii.gz')
        # # _, file_extension = os.path.splitext(final_mask)

        # Create the new filename with the number of voxels included
        new_count_filename = f'{outDir_nat}/{area}_masks/phzstat2_loczstat3_{area}_thr{th}_numvoxels{num_voxels}.nii.gz'

        # Rename the file with the new filename
        os.rename(area_mask_native, new_count_filename)