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
utils = '/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019/code/utils'
sys.path.append(op.expanduser(f'{utils}'))
from experiment import experiment

# set system to process NIFTI_GZ
os.system('export FSLOUTPUTTYPE=NIFTI_GZ')
resultsFile = f'/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019/data/functional/localizer/psc/timeseries_allconditions.csv'
results = open(resultsFile, 'w')
results.write('topup,b0,HRFmodel,scan,subject,session,region,condition,run,rep,TR,PSC\n')
for c, condition in enumerate(['Inducers','Phantom']): #make sure conditions are in same order as the event files stored in event_files
	for topup in ['topUp']:#, 'noTopUp']:
		for b0 in ['noB0']:#,'b0']:
			for HRFmodel in ['doubleGamma']:#, 'singleGamma', 'singleGamma034']:
				fieldStrength = list(experiment['scanInfo'].keys())[0]
				for subject in experiment['scanInfo'][fieldStrength].keys():
					for session in experiment['scanInfo'][fieldStrength][subject].keys():
						sessDir = f'/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}'
						for scan in experiment['scanInfo'][fieldStrength][subject][session]['funcScans'].keys():
							if scan != 'restingState' and scan != 'retinotopy' and scan != 'phantom':
								params = experiment['design'][scan]['params']
								for run in range(len(experiment['scanInfo'][fieldStrength][subject][session]['funcScans'][scan])):

									'''
									# for debugging
									distCor = 'topUp'
									HRFmodel = 'doubleGamma'
									scan = 'figureGround_v3'
									params = experiment['design'][scan]['params']
									subject = 'M015'
									session = '201215'
									run=0
									'''

									print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | topup : {topup} | b0: {b0} | HRF model: {HRFmodel} | Scan: {scan} | Subject: {subject} | Session: {session} | Run: {run+1}')
									runDir = os.path.join('/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal', subject, 'data/functional', scan,
														  f'run{run+1:02}/outputs', 'firstLevel.feat')
									os.makedirs(os.path.join(runDir, 'timeseries',condition), exist_ok=True)

									if not os.path.exists(f'{runDir}/timeseries/{condition}/filtered_func_percentSC.nii.gz') or overwrite:
										print('calculating mean...')
										os.system(f'fslmaths {runDir}/filtered_func_data -Tmean {runDir}/timeseries/{condition}/Tmean')
										print('demeaning...')
										os.system(f'fslmaths {runDir}/filtered_func_data -sub {runDir}/timeseries/{condition}/Tmean {runDir}/timeseries/{condition}/filtered_func_demeaned')
										print('converting to percent MR signal...')
										os.system(f'fslmaths {runDir}/timeseries/{condition}/filtered_func_demeaned -div {runDir}/timeseries/{condition}/Tmean -mul 100 {runDir}/timeseries/{condition}/filtered_func_percentSC')


									# create timeseries
									print('creating timeseries...')
									masks = sorted(glob.glob(os.path.join(sessDir, f'data/functional/{scan}/run{run+1:02}/outputs/masksNative/*.nii.gz')))
									for mask in masks:
										maskName = os.path.basename(mask)[:-7]
										print(maskName)

										# get mean time-series of mask for whole scan
										if not os.path.exists(f'{runDir}/timeseries/{condition}/{maskName}.txt') or overwrite:
											os.system(f'fslmeants -i {runDir}/timeseries/{condition}/filtered_func_percentSC -o {runDir}/timeseries/{condition}/{maskName}.txt -m {mask}')

										f = open(f'{runDir}/timeseries/{condition}/{maskName}.txt')
										ts_output = f.readlines()
										f.close()

										if scan == 'localizer':
											event_files = sorted(glob.glob(os.path.join('/home/tonglab/Documents/Loic/phantom_mri_data/derivatives/FEAT/events_3column', subject, f'ses-1/task-{scan}','*.txt')))
											# conds = []
											start_times = []
											for event_file in event_files:
												# conds.append(os.path.basename(event_file)[15:-4])
												# open event file and get start times
												f = open(event_file)
												events = f.readlines()
												start_times.append([row.split()[0] for row in events])
												f.close()

										# elif scan == 'localizer':
										# 	start_times = [list(range(12,300,48)),list(range(36,300,48))]
										# 	conds = ['center','surround']

										# for c, cond in enumerate(conds):
										start_time = start_times[c]
										for rep in range(len(start_time)):
											for vol in range(int((params['blockDuration']+params['IBI'])/params['TR']+6)): # adding 6 for 3 vols either side of epoch

												thisVol = int(start_time[rep])/params['TR'] + vol - 3 # subtract 3 to begin 3 vols before onset (zeroth vol recorded between 0:1 TRs so start at -3)
												temp = None
												if thisVol < len(ts_output):
													temp = ts_output[int(thisVol)].split(" ")[0]
												results.write(f'{topup},{b0},{HRFmodel},{scan},{subject},{session},{maskName},{condition},{run+1},{rep+1},{vol-3},{temp}\n')

results.close()
print('FINISHED')
