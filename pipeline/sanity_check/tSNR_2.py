import numpy as np
import os
import os.path as op
# 20230620
# in this version, the preprocessing is more and done with FSL
reg_dir = '/home/tonglab/Miao/fMRI/figureGround/fMRI/data/individual/F013/temp_tsnr_restingState/230804'
ref_func_t= '/home/tonglab/Miao/fMRI/figureGround/fMRI/data/individual/F013/230804/functional/resting_state/run01/inputs/rawData'

# reg_dir = '/home/tonglab/Miao/fMRI/figureGround/fMRI/data/individual/F013/temp_tsnr_restingState/2211_V4'
# ref_func_t = '/home/tonglab/Miao/fMRI/figureGround/fMRI/data/individual/F013/2211/sub-F013_ses-1_task-restingState_dir-AP_run-1_bold'

ref_func_Preprocessed = reg_dir+'/filtered_func_data'
ref_func = reg_dir+'/rawData'

### remove the skull
overwrite = 0
#### std_anatomical = '/home/tonglab/fsl/data/standard/MNI152_T1_2mm_brain'
ref_anat = '/home/tonglab/freesurfer/subjects/F013/mri/orig/anatomical_brain'
fs_dir = '/home/tonglab/freesurfer/subjects/F013'
subject = 'F013'
ref_anat_brain = '/home/tonglab/Miao/fMRI/figureGround/fMRI/data/individual/F013/temp_tsnr_restingState/anatomical_brain_betExtracted'
if not op.isfile(ref_anat_brain):
    os.system(f'bet {ref_anat} {ref_anat_brain}')

### do preprocessing and average resting state in fsl and then preceed
# if not os.path.exists(ref_func_Preprocessed) or overwrite:
#     os.system(f'mcflirt -in {ref_func_t} -out {ref_func_Preprocessed}')
os.system(f'fslmaths {ref_func_Preprocessed} -Tmean {ref_func}')

### from function to native anatomical
xform_func_to_anat = f'{reg_dir}/func_to_anat.mat'
if not op.isfile(xform_func_to_anat):
    os.system(
        f'epi_reg --epi={ref_func} --t1={ref_anat} --t1brain={ref_anat_brain} --out={xform_func_to_anat[:-4]}')  # BBR
xform_anat_to_func = f'{reg_dir}/anat_to_func.mat'
if not op.isfile(xform_anat_to_func):
    os.system(f'convert_xfm -omat {xform_anat_to_func} -inverse {xform_func_to_anat}')

### from native anatomical to standard anatomical
xform_std_to_anat = f'{fs_dir}/mri/transforms/reg.mni152.2mm.lta'
if not op.isfile(xform_std_to_anat):
    os.system(f'mni152reg --s {subject}')

### TSNR
tSNRdir = reg_dir +'/tsnr'
pathTmean = os.path.join(tSNRdir, "Tmean.nii.gz")
pathTstd = os.path.join(tSNRdir, "Tstd.nii.gz")
pathTSNR = os.path.join(tSNRdir, "tSNR.nii.gz")
if not os.path.exists(pathTmean) or overwrite:
    print('Creating mean of timeseries...')
    os.system(f'fslmaths {ref_func_Preprocessed} -Tmean {pathTmean}')
if not os.path.exists(pathTstd) or overwrite:
    print('Creating std of timeseries...')
    os.system(f'fslmaths {ref_func_Preprocessed} -Tstd {pathTstd}')
if not os.path.exists(pathTSNR) or overwrite:
    print('Calculating tSNR map...')
    os.system(f'fslmaths {pathTmean} -div {pathTstd} {pathTSNR}')

resultsFile_name = tSNRdir + '/results_name.csv'
resultsFile_mean = tSNRdir + '/results_mean.csv'
resultsFile_std = tSNRdir + '/results_std.csv'
resultsFile_number = tSNRdir + '/results_number.csv'


### invert mask
maskName =['V1','V2','V3a','V3b','V3','hV4','LO1','LO2','VO1','VO2','V1-V3','V1-V4']
results = open(resultsFile_name, 'a+')
for i in range(len(maskName)):
    results.write(f'{maskName[i]}\n')
    mask_path_stds = '/home/tonglab/Miao/fMRI/masks/Wang_2015'
    mask_path_std = mask_path_stds +'/'+ maskName[i]+'.nii.gz'
    mask_path_anat = reg_dir +'/anatomicalMask/'+ maskName[i]+'.nii.gz'
    mask_path_func  = reg_dir +'/functionalMask/'+ maskName[i]+'.nii.gz'
    os.system(
        f'mri_vol2vol --mov {mask_path_std} --targ {ref_anat}.nii.gz --lta {xform_std_to_anat} --o {mask_path_anat} --nearest')

    os.system(
        f'flirt -in {mask_path_anat} -ref {ref_func}.nii.gz -applyxfm -init {xform_anat_to_func} -out {mask_path_func} -interp nearestneighbour')

    ### save by region

    os.system(f'fslstats {pathTSNR} -k {mask_path_func} -M >> {resultsFile_mean}')
    os.system(f'fslstats {pathTSNR} -k {mask_path_func} -S >> {resultsFile_std}')
    os.system(f'fslstats {pathTSNR} -k {mask_path_func} -V | tr \' \' \',\' >> {resultsFile_number}')

results.close()

#fslstats /home/tonglab/Miao/fMRI/figureGround/fMRI/data/individual/F013/temp_tsnr_restingState/230804/tsnr/tSNR.nii.gz -k /home/tonglab/Miao/fMRI/figureGround/fMRI/data/individual/F013/temp_tsnr_restingState/230804/functionalMask/V1.nii.gz -M > /home/tonglab/Miao/fMRI/figureGround/fMRI/data/individual/F013/temp_tsnr_restingState/230804/tsnr/results.txt