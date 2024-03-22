#!/usr/bin/python
"""
Loic - 07/17/2023
get data for timeseries of each condition (in % signal change) and place in csv file

-convert filtered_func_data to % signal change
-separate by condition based on event files
-get timeseries for each condition * mask

"""

import os
import os.path as op
import sys
import glob
import pickle
import datetime

overwrite = True

# get scan info from experiment file
utils = '/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/code/utils'
sys.path.append(op.expanduser(f'{utils}'))
from experiment import experiment

# set system to process NIFTI_GZ
os.system('export FSLOUTPUTTYPE=NIFTI_GZ')


for topup in ['topUp']:#, 'noTopUp']:
	for b0 in ['noB0']:#,'b0']:
		for HRFmodel in ['doubleGamma']:#, 'singleGamma', 'singleGamma034']:
			fieldStrength = list(experiment['scanInfo'].keys())[0]
			for subject in experiment['scanInfo'][fieldStrength].keys():
				for session in experiment['scanInfo'][fieldStrength][subject].keys():
					sessDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}'
					#for scan in experiment['scanInfo'][fieldStrength][subject][session]['funcScans'].keys():
						#if scan == 'phantom':
					scan = 'phantom'
					params = experiment['design'][scan]['params']
					for run in range(len(experiment['scanInfo'][fieldStrength][subject][session]['funcScans'][scan])):


						# for debugging
						# topup = 'topUp'
						# b0 = 'noB0'
						# HRFmodel = 'doubleGamma'
						# scan = 'phantom'
						# params = experiment['design'][scan]['params']
						# subject = 'sub-F019'
						# session = '230321'
						# run=0


						print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | topup : {topup} | b0: {b0} | HRF model: {HRFmodel} | Scan: {scan} | Subject: {subject} | Session: {session} | Run: {run+1} ')
						subject_sess = f'{subject}_{session}'
						runDir = os.path.join('/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal', subject_sess, 'data/functional', scan, f'run{run+1:02}/outputs/topUp/noB0', 'preprocessing.feat')

						os.makedirs(os.path.join(runDir, 'timeseries'), exist_ok=True)
						if not os.path.exists(f'{runDir}/timeseries/filtered_func_percentSC.nii.gz'): #or overwrite:
							print('calculating mean...')
							os.system(f'fslmaths {runDir}/filtered_func_data -Tmean {runDir}/timeseries/Tmean')
							print('demeaning...')
							os.system(f'fslmaths {runDir}/filtered_func_data -sub {runDir}/timeseries/Tmean {runDir}/timeseries/filtered_func_demeaned')
							print('converting to percent MR signal...')
							os.system(f'fslmaths {runDir}/timeseries/filtered_func_demeaned -div {runDir}/timeseries/Tmean -mul 100 {runDir}/timeseries/filtered_func_percentSC')

# create timeseries
for area in ['V1_lh','V1_rh', 'V2_lh', 'V2_rh', 'V3_lh', 'V3_rh', 'hV4_lh', 'hV4_rh', 'V3a_lh', 'V3a_rh', 'V3b_lh', 'V3b_rh']: #['V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b']:  # ,
	for zstat in ['zstat2']: #'zstat4'
		print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Area: {area}| zstat: {zstat}')
		os.makedirs(os.path.join('/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/phantom/allRuns/outputs/topUp/noB0/psc/phEnhancement_roi'), exist_ok=True)
		resultsFile = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/phantom/allRuns/topUp/noB0/psc/phEnhancement_roi/timeseries_{area}_ph{zstat}.csv'
		results = open(resultsFile, 'w')
		results.write('topup,b0,HRFmodel,scan,subject,session,region,condition,run,rep,TR,PSC\n')
		for topup in ['topUp']:  # , 'noTopUp']:
			for b0 in ['noB0']:  # ,'b0']:
				for HRFmodel in ['doubleGamma']:  # , 'singleGamma', 'singleGamma034']:
					fieldStrength = list(experiment['scanInfo'].keys())[0]
					for subject in experiment['scanInfo'][fieldStrength].keys():
						for session in experiment['scanInfo'][fieldStrength][subject].keys():
							sessDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}'
							#for scan in experiment['scanInfo'][fieldStrength][subject][session]['funcScans'].keys():
								#if scan == 'phantom':
							scan = 'phantom'
							params = experiment['design'][scan]['params']
							for run in range(len(experiment['scanInfo'][fieldStrength][subject][session]['funcScans'][scan])):
								'''
								
								# For Debug
								area = 'V1'
								zstat = 'zstat2'
								resultsFile = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019_230321/data/functional/phantom/allRuns/topUp/noB0/psc/timeseries_{area}_{zstat}.csv'
								results = open(resultsFile, 'w')
								results.write('topup,b0,HRFmodel,scan,subject,session,region,condition,run,rep,TR,PSC\n')
								topup = 'topUp'
								b0 = 'noB0'
								HRFmodel = 'doubleGamma'
								fieldStrength = '7T'
								subject = 'sub-F019'
								session = 230321
								sessDir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}'
								scan = 'phantom'
								params = experiment['design'][scan]['params']
								run = 0
								
								'''

								print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | topup : {topup} | b0: {b0} | HRF model: {HRFmodel} | Scan: {scan} | Subject: {subject} | Session: {session} | Run: {run+1} ')
								subject_sess = f'{subject}_{session}'
								runDir = os.path.join('/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal', subject_sess,'data/functional', scan, f'run{run + 1:02}/outputs/topUp/noB0', 'preprocessing.feat')

								print('creating timeseries...')

								masks = sorted(glob.glob(os.path.join(sessDir, f'data/functional/phantom/allRuns/outputs/topUp/noB0/phantomEnhanceMasksNative/{area}_masks/phzstat2_loczstat3_{area}*.nii.gz')))
								os.makedirs(os.path.join(runDir, 'timeseries/enhancement_masks'), exist_ok=True)
								for mask in masks:
									maskName = os.path.basename(mask)[:-7]
									print(maskName)

									# get mean time-series of mask for whole scan
									if not os.path.exists(f'{runDir}/timeseries/enhancement_masks/{maskName}.txt') or overwrite:
										os.system(f'fslmeants -i {runDir}/timeseries/filtered_func_percentSC -o {runDir}/timeseries/enhancement_masks/{maskName}.txt -m {mask}')

									f = open(f'{runDir}/timeseries/enhancement_masks/{maskName}.txt')
									ts_output = f.readlines()
									f.close()

									for c, condition in enumerate(['vertSingleLeft', 'horSingleLeft', 'vertSingleRight', 'horSingleRight', 'horDoubleIndir']):
										# event_file = sorted(glob.glob(os.path.join('/Users/tonglab/Documents/Loic/phantom_mri_data/derivatives/FEAT/events_3column', subject, f'ses-1/task-{scan}', f'*run{run+1:03}_{condition}.txt'))) # make sure conditions are in same order as the event files stored in event_files
										event_file = os.path.join('/Users/tonglab/Documents/Loic/phantom_mri_data/derivatives/FEAT/events_3column',subject_sess, f'ses-1/task-{scan}', f'{subject}_ses-01_task-{scan}_run-{run + 1:03}_{condition}.txt')
										start_times = []
										# open event file and get start times
										f = open(event_file)
										events = f.readlines()
										start_times.append([row.split()[0] for row in events])
										f.close()
										#os.makedirs(os.path.join(runDir, 'timeseries', condition),exist_ok=True)
										start_time = start_times[0]
										for rep in range(len(start_time)):
											for vol in range(int((params['blockDuration']+params['IBI'])/params['TR']+6)): # adding 6 for 3 vols either side of epoch

												thisVol = int(start_time[rep])/params['TR'] + vol - 3 # subtract 3 to begin 3 vols before onset (zeroth vol recorded between 0:1 TRs so start at -3)
												temp = None
												if thisVol < len(ts_output):
													temp = ts_output[int(thisVol)].split(" ")[0]
												results.write(f'{topup},{b0},{HRFmodel},{scan},{subject},{session},{maskName},{condition},{run+1},{rep+1},{vol-3},{temp}\n')

		results.close()
print('FINISHED')
