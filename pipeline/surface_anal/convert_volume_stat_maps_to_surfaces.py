
'''
Loic - 07/12/2023
This script is developed to generate surface map versions of the statistical maps generated in FEAT outputs.
'''
import os
import subprocess
import os.path as op

subject = 'sub-F019'
session = '230321'

utils = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/code/utils'
sys.path.append(op.expanduser(f'{utils}'))

#phantom data dirs
workDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/second_level.gfeat/cope2.feat'
regDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/reg'
outDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/phantomEnhanceMapsStandard'
outDir_nat = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/phantomEnhanceMapsNative'
outDir_nat_surf = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/phantomEnhanceSurfsNative'
os.makedirs(f'{outDir_nat_surf}', exist_ok=True)
os.makedirs(f'{outDir_nat}', exist_ok=True)
os.makedirs(f'{outDir}', exist_ok=True)
phDirFeat = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/run01/outputs/topUp/noB0/preprocessing.feat'
dataDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data'

#create stat map surface of phantom data for whole hemisphere
for h, hemi in enumerate(['lh', 'rh']):

    os.makedirs(f'{outDir}/{hemi}_maps', exist_ok=True)
    os.makedirs(f'{outDir_nat}/{hemi}_maps', exist_ok=True)
    os.makedirs(f'{outDir_nat_surf}/{hemi}_surfs', exist_ok=True)
    # for z in [2, 3]: #[1, 2, 3, 4, 5]
    for th in [3.1, 3.6]:
        phantom_zstat = f'{workDir}/stats/zstat1.nii.gz'
        thrZstat = f'{workDir}/stats/modif_zstat_maps/thresh{th}_zstat1.nii.gz'
        os.makedirs(f'{workDir}/stats/modif_zstat_maps', exist_ok=True)
        # ##Create thresholded z-stat map
        os.system(f'fslmaths {phantom_zstat} -thr {th} {thrZstat}')
        #
        # # The fslstats command to get the number of non-zero voxels (use '-l 0' to exclude background value)
        command = f'fslstats {thrZstat} -V -l 0'
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
        new_filename = f'{outDir}/{hemi}_maps/phzstat2_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz'

        # Rename the file with the new filename
        os.rename(thrZstat, new_filename)

        ## Transform map to native space
        map_native = f'{outDir_nat}/{hemi}_maps/phzstat2_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz' #output
        refFunc = f'{phDirFeat}/example_func.nii.gz' #reference
        trans_mat = f'{phDirFeat}/reg/standard2example_func.mat' #transformation matrix
        os.system(f'flirt -in {new_filename} -ref {refFunc} -out {map_native} -init {trans_mat} -interp nearestneighbour -applyxfm')

        #recount number of voxels
        command = f'fslstats {map_native} -V -l 0'
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

        # Create the new filename with the number of voxels included
        new_count_filename = f'{outDir_nat}/{hemi}_maps/phzstat2_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz'

        # Rename the file with the new filename
        os.rename(map_native, new_count_filename)

        surfOutDir = f'{outDir_nat_surf}/{hemi}_surfs/phzstat2_{hemi}_thr{th}_surf.nii.gz'
        freesurferCommand = f'mri_vol2surf --mov {new_count_filename} --hemi {hemi} --o {surfOutDir} --regheader {subject}_{session} --projfrac 0.5' # --reg {dataDir}/sess01/bold/register.dof6.lta not needed as already in native space {regDir}/fsAnat2origAnat.mat
        os.system(f'bash {utils}/callFreesurferFunction.sh -s "{freesurferCommand}"')

#Localizer dirs
locDir = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/preprocessing.feat'
locOutDir = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/MapsStandard'
locOutDir_nat = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/MapsNative'
locOutDir_nat_surf = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/localizer/run01/outputs/topUp/noB0/MapsNativeSurfs'
os.makedirs(f'{locOutDir}', exist_ok=True)
os.makedirs(f'{locOutDir_nat}', exist_ok=True)
os.makedirs(f'{locOutDir_nat_surf}', exist_ok=True)
#create stat map surface of localizer data for whole hemisphere
for h, hemi in enumerate(['lh', 'rh']):

    os.makedirs(f'{locOutDir}/{hemi}_maps', exist_ok=True)
    os.makedirs(f'{locOutDir_nat}/{hemi}_maps', exist_ok=True)
    os.makedirs(f'{locOutDir_nat_surf}/{hemi}_surfs', exist_ok=True)
    # for z in [2, 3]: #[1, 2, 3, 4, 5]
    for th in [3.1, 3.6]:
        phantom_zstat = f'{locDir}/stats/zstat3.nii.gz'
        thrZstat = f'{locDir}/stats/modif_zstat_maps/thresh{th}_zstat3.nii.gz'
        os.makedirs(f'{locDir}/stats/modif_zstat_maps', exist_ok=True)
        # ##Create thresholded z-stat map
        os.system(f'fslmaths {phantom_zstat} -thr {th} {thrZstat}')
        #
        # # The fslstats command to get the number of non-zero voxels (use '-l 0' to exclude background value)
        command = f'fslstats {thrZstat} -V -l 0'
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
        new_filename = f'{locOutDir}/{hemi}_maps/loczstat3_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz'

        # Rename the file with the new filename
        os.rename(thrZstat, new_filename)

        # ## Transform map from standard to functional native space
        # map_native = f'{locOutDir_nat}/{hemi}_maps/loczstat3_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz' #output
        # refFunc = f'{locDir}/example_func.nii.gz' #reference
        # trans_mat = f'{locDir}/reg/standard2example_func.mat' #transformation matrix
        # os.system(f'flirt -in {new_filename} -ref {refFunc} -out {map_native} -init {trans_mat} -interp nearestneighbour -applyxfm')
        #
        # #recount number of voxels
        # command = f'fslstats {map_native} -V -l 0'
        # #
        # # # Run the command and capture the output
        # output = subprocess.check_output(command, shell=True, text=True)
        # #
        # # # The output contains two values: volume and number of voxels
        # num_voxels, volume = output.strip().split()
        # #
        # # # Convert the number of voxels to an integer and store it in a Python variable
        # num_voxels = int(num_voxels)
        # #
        # # # Print or use the value as needed
        # print("Volume:", volume)
        # print("Number of Voxels:", num_voxels)
        #
        # # Create the new filename with the number of voxels included
        # new_count_filename = f'{locOutDir_nat}/{hemi}_maps/loczstat3_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz'
        #
        # # Rename the file with the new filename
        # os.rename(map_native, new_count_filename)

        surfOutDir = f'{locOutDir_nat_surf}/{hemi}_surfs/loczstat3_{hemi}_thr{th}_surf.nii.gz'
        freesurferCommand = f'mri_vol2surf --mov {new_filename} --hemi {hemi} --o {surfOutDir} --regheader {subject}_{session} --projfrac 0.5' # --reg {dataDir}/sess01/bold/register.dof6.lta not needed as already in native space {regDir}/fsAnat2origAnat.mat
        os.system(f'bash {utils}/callFreesurferFunction.sh -s "{freesurferCommand}"')

#create surface by area
for a, area in enumerate(['V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b']): #enumerate(['V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b']): 'V1_lh','V1_rh', 'V2_lh', 'V2_rh', 'V3_lh', 'V3_rh', 'hV4_lh', 'hV4_rh', 'V3a_lh', 'V3a_rh', 'V3b_lh', 'V3b_rh'
    for h,hemi in enumerate(['lh', 'rh']):
        area_mask_standard = f'/Users/tonglab/Documents/Loic/atlases/Wang_2015/{area}_{hemi}.nii.gz' #input

        os.makedirs(f'{outDir}/{area}_{hemi}_maps', exist_ok=True)
        os.makedirs(f'{outDir_nat}/{area}_{hemi}_maps', exist_ok=True)
        os.makedirs(f'{outDir_nat_surf}/{area}_{hemi}_surfs', exist_ok=True)
        # for z in [2, 3]: #[1, 2, 3, 4, 5]
        for th in [3.1, 3.6]:
            phantom_zstat = f'{workDir}/stats/zstat1.nii.gz'
            # zstat = glob.glob(f'{workDir}/stats/improved_zstat_maps/zstat{z}_sup_*.nii.gz')
            thrZstat = f'{workDir}/stats/modif_zstat_maps/thresh{th}_zstat1.nii.gz'
            os.makedirs(f'{workDir}/stats/modif_zstat_maps', exist_ok=True)
            # ##Create thresholded z-stat map
            os.system(f'fslmaths {phantom_zstat} -thr {th} {thrZstat}')

            final_map = f'{outDir}/{area}_{hemi}_maps/phzstat2_{area}_{hemi}_thr{th}.nii.gz' #phzstats2 since we use the second contrast of phantom runs (cope2.feat)
            os.system(f'fslmaths {thrZstat} -mul {area_mask_standard} {final_map}')
            #
            # # The fslstats command to get the number of non-zero voxels (use '-l 0' to exclude background value)
            command = f'fslstats {final_map} -V -l 0'
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
            new_filename = f'{outDir}/{area}_{hemi}_maps/phzstat2_{area}_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz'

            # Rename the file with the new filename
            os.rename(final_map, new_filename)

            ## Transform map to native space
            area_map_native = f'{outDir_nat}/{area}_{hemi}_maps/phzstat2_{area}_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz' #output
            refFunc = f'{phDirFeat}/example_func.nii.gz' #reference
            trans_mat = f'{phDirFeat}/reg/standard2example_func.mat' #transformation matrix
            os.system(f'flirt -in {new_filename} -ref {refFunc} -out {area_map_native} -init {trans_mat} -interp nearestneighbour -applyxfm')

            #recount number of voxels
            command = f'fslstats {area_map_native} -V -l 0'
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
            new_count_filename = f'{outDir_nat}/{area}_{hemi}_maps/phzstat2_{area}_{hemi}_thr{th}_numvoxels{num_voxels}.nii.gz'

            # Rename the file with the new filename
            os.rename(area_map_native, new_count_filename)

            # try:
            #     subprocess.run(f'source /Applications/freesurfer/7.4.1/SetUpFreeSurfer', shell=True, check=True)
            # except subprocess.CalledProcessError as e:
            #     print(f"Error: {e}")

            surfOutDir = f'{outDir_nat_surf}/{area}_{hemi}_surfs/phzstat2_{area}_{hemi}_thr{th}_surf.nii.gz'
            freesurferCommand = f'mri_vol2surf --mov {new_count_filename} --hemi {hemi} --o {surfOutDir} --projfrac 0.5' # --reg {dataDir}/sess01/bold/register.dof6.lta not needed as already in native space {regDir}/fsAnat2origAnat.mat
            os.system(f'bash {utils}/callFreesurferFunction.sh -s "{freesurferCommand}"')

# surface = f'{out_dir}/FOV_mask_{hemi}.mgh'
# if not op.isfile(surface) or overwrite:
#     cmd = f'mri_vol2surf --mov {FOV_mask} --out {surface} --regheader {subject} --hemi {hemi}'
#     os.system(cmd)
# # surface to label
# label = f'{surface[:-4]}.label'
# if not op.isfile(label) or overwrite:
#     cmd = f'mri_cor2label --i {surface} --surf {subject} {hemi} --id 1 --l {label}'
#     os.system(cmd)
