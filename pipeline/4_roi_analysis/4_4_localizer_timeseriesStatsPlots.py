'''
runs featquery analysis
'''

import os
import os.path as op
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# get scan info from experiment file
utils = '/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019/code/utils'
sys.path.append(op.expanduser(f'{utils}'))
from experiment import experiment
#buffer = 3 # 3 volumes either side of epoch

# get conditions
#relOri = ['-90', '-45', '-15', '0', '15', '45', '90']
contrasts = ['Both conditions', 'Inducer patches only', 'Phantom patch only', 'Phantom > Inducers', 'Inducers > Phantom']
figSize = (4, 3)
dataFile = f'/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019/data/functional/localizer/psc/timeseries_allconditions.csv'
data = pd.read_csv(dataFile)
# for condition in ['Inducers','Phantom']:

# add columns for region, hemisphere and voxel counts
regions, hemispheres, nVoxels, masks = [[], [], [], []]
for row in data.index:
	regions.append(data['region'][row].split(sep='_')[1])
	masks.append(data['region'][row].split(sep='_')[0])
	# hemispheres.append(data['region'][row].split(sep='_')[1]) #use these lines if the mask has the name regions_hemisphere_nVoxels
	# nVoxels.append(data['region'][row].split(sep='_')[2])
data['region'] = regions
data['mask'] = masks
#data['hemisphere'] = hemispheres
#data['nVoxels'] = nVoxels

# normed PSC
baseWindow=np.array([-1]).astype(int)


for topup in ['topUp']:#, 'noTopUp']:
	for b0 in ['noB0']:#,'b0']:
		for HRFmodel in ['doubleGamma']:#, 'singleGamma', 'singleGamma034']:
			fieldStrength = list(experiment['scanInfo'].keys())[0]
			for subject in list(experiment['scanInfo'][fieldStrength].keys()):
				for session in experiment['scanInfo'][fieldStrength][subject].keys():


					'''
					# debugging
					distCor = 'topUp'
					HRFmodel = 'doubleGamma'
					subject = 'M015'
					session = 201215
					'''

					# EXPERIMENT
					scan = 'localizer'
					params = experiment['design'][scan]['params']
					epochDur = int((params['blockDuration'] + params['IBI']) / params['TR'])
					buffer = 3  # 3 volumes either side of epoch
					xticks = np.arange(-buffer, epochDur + buffer + 1) * params[
						'TR']  # time points fall in between measurements
					x_pos = xticks[:-1] + 1  # points are shifted to centre on middle of TR


					outDir = os.path.join('/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019/analysis/results', topup, b0, HRFmodel, subject, session, scan)
					os.makedirs(outDir, exist_ok=True)
					dataConf = data[(data['topup'] == topup) &
									(data['b0'] == b0) &
									(data['HRFmodel'] == HRFmodel) &
									(data['subject'] == subject) &
									(data['session'] == int(session)) &
									(data['scan'] == scan)].reset_index()

					# baseline
					lastRow = 0
					PSCbaselined = []
					while lastRow < len(dataConf):
						theseRows = np.arange(lastRow, lastRow + 15)
						baselineVals = dataConf['PSC'][theseRows[tuple([baseWindow + 3])]].astype(float)
						baseline = np.mean(baselineVals)
						for row in theseRows:
							if dataConf['PSC'][row] != 'None':
								PSCbaselined.append(float(dataConf['PSC'][row]) - baseline)
							else:
								PSCbaselined.append('None')
						lastRow = row + 1
					dataConf = dataConf.assign(PSCbaselined=PSCbaselined)


					# make plots for each configuration of parameters
					# for region in dataConf['region'].unique():
					# 	dataRegion = dataConf[dataConf['region'] == region]
						# for size in dataRegion['nVoxels'].unique():
						# 	dataSize = dataRegion[dataRegion['nVoxels'] == size]

							# plt.figure(figsize=figSize)
							#
							# # for ori in relOri:
							# means = np.empty([epochDur + buffer * 2])
							# errorbar = np.empty([epochDur + buffer * 2])
							# for timepoint in range(epochDur + buffer * 2):
							# 	values_exp = np.array(dataSize['PSCbaselined'][(dataSize['condition']==ori) & (dataSize['TR'] == timepoint - buffer)])
							# 	values_exp = np.delete(values_exp, [x == 'None' for x in values_exp])
							# 	values_exp = values_exp.astype(float)
							# 	means[timepoint] = np.mean(values_exp)
							# 	errorbar[timepoint] = stats.sem(values_exp)
							# plt.errorbar(x_pos,
							# 			 means,
							# 			 yerr=errorbar,
							# 			 marker='.',
							# 			 markersize=5,
							# 			 label=ori)
							# plt.legend(bbox_to_anchor=(1.04, .75), borderaxespad=0)  # put legend outside plot
							# plt.xticks(xticks)
							# plt.xlabel('time (s)')
							# plt.ylabel('signal change (%)')
							# plt.title(f'timeseries: {region}, {size} voxels')
							# #plt.tight_layout(rect=[0, 0, 0.75, 1])  # ensure everything is placed in the canvas
							# plt.subplots_adjust(right=1.4)  # allow space for legend
							# plt.savefig(os.path.join(outDir, f'{region}_{size}_timeseries.png'), bbox_inches='tight')
							# plt.show()
							# plt.close()
							#
							# # PARAMETER ESTIMATES FROM TIMESERIES
							# # alternative to featquery analyses
							# timeWindow = (8, 14)  # in seconds, inclusive
							# volWindow = np.array(timeWindow) / params['TR']
							# vols2use = np.arange(volWindow[0], volWindow[1] + 1)
							# dataWindow = dataSize[(dataSize['TR'].isin(vols2use))]
							# dataWindow = dataWindow.astype({'PSCbaselined': 'float'})  # convert PSC column to float
							# #dataWindowMeanAcrossReps = dataWindow.groupby(['orientation', 'subject']).mean()
							# #dataWindowMeanAcrossReps = dataWindowMeanAcrossReps[['PSC']]
							# dataWindowMeans = dataWindow.groupby(['condition']).mean()
							# dataWindowSems = dataWindow.groupby(['condition']).sem()
							#
							# means = []
							# sems = []
							# condsFull = []
							# plt.figure(figsize=(5,3))
							# for o, ori in enumerate(dataWindowMeans.index):
							# 	condsFull.append(f'{ori}')
							# 	means.append(dataWindowMeans['PSCbaselined'][ori])
							# 	sems.append(dataWindowSems['PSCbaselined'][ori])
							# #barColors = ['red', 'darkred', 'blue', 'darkblue']
							# plt.errorbar(relOri,
							# 			 means,
							# 			 yerr=sems,
							# 			 marker='.',
							# 			 markersize=5)
							# plt.xticks(relOri)
							# plt.xlabel('relative orientation')
							# plt.ylabel('signal change (%)')
							# plt.tight_layout()  # ensure everything is placed in the canvas
							# plt.savefig(f'{outDir}/{region}_{size}_PSCts.png')
							# plt.show()
							# plt.close()

					# LOCALIZER
					scan = 'localizer'
					params = experiment['design']['localizer']['params']
					epochDur = int((params['blockDuration'] + params['IBI']) / params['TR'])
					buffer = 3  # 3 volumes either side of epoch
					xticks = np.arange(-buffer, epochDur + buffer + 1) * params[
						'TR']  # time points fall in between measurements
					x_pos = xticks[:-1] + 1  # points are shifted to centre on middle of TR

					outDir = os.path.join('/home/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/sub-F019/analysis/results', topup, b0, HRFmodel, subject, session, scan)
					os.makedirs(outDir, exist_ok=True)
					dataConf = data[(data['topup'] == topup) &
									(data['b0'] == b0) &
									(data['HRFmodel'] == HRFmodel) &
									(data['subject'] == subject) &
									(data['session'] == int(session)) &
									(data['scan'] == scan)].reset_index()

					# baseline
					lastRow = 0
					PSCbaselined = []
					while lastRow < len(dataConf):
						theseRows = np.arange(lastRow, lastRow + 18)
						baselineVals = dataConf['PSC'][theseRows[tuple([baseWindow + 3])]].astype(float)
						baseline = np.mean(baselineVals)
						for row in theseRows:
							if dataConf['PSC'][row] != 'None':
								PSCbaselined.append(float(dataConf['PSC'][row]) - baseline)
							else:
								PSCbaselined.append('None')
						lastRow = row + 1
					dataConf = dataConf.assign(PSCbaselined=PSCbaselined)

					# make plots for each configuration of parameters
					for region in dataConf['region'].unique():
						dataRegion = dataConf[dataConf['region'] == region]
						# for size in dataRegion['nVoxels'].unique():
						# 	dataSize = dataRegion[dataRegion['nVoxels'] == size]
						for m, mask in enumerate(dataConf['mask'].unique()):
							plt.figure(figsize=figSize)
							for stiLoc in ['Inducers','Phantom']:
								means = np.empty([epochDur + buffer * 2])
								errorbar = np.empty([epochDur + buffer * 2])
								for timepoint in range(epochDur + buffer * 2):
									values_loc = np.array(dataRegion['PSCbaselined'][(dataRegion['condition']==stiLoc)&(dataRegion['TR'] == timepoint - buffer)&(dataRegion['mask']==mask)])
									values_loc = np.delete(values_loc, np.isnan(values_loc))
									values_loc = values_loc.astype(float)
									means[timepoint] = np.mean(values_loc)
									errorbar[timepoint] = stats.sem(values_loc)
								plt.errorbar(x_pos,
											 means,
											 yerr=errorbar,
											 marker='.',
											 markersize=5,
											 label = stiLoc)
							plt.legend(bbox_to_anchor=(1.04, .75), borderaxespad=0)  # put legend outside plot
							plt.xticks(xticks)
							plt.xlabel('time (s)')
							plt.ylabel('signal change (%)')
							plt.title(f'timeseries: {region}, {contrasts[m]} ') #{size} voxels
							#plt.tight_layout(rect=[0, 0, 0.75, 1])  # ensure everything is placed in the canvas
							plt.subplots_adjust(right=1.4)  # allow space for legend
							plt.savefig(os.path.join(outDir, f'{region}_{mask}_timeseries.png'), bbox_inches='tight')#{size}
							plt.show()
							plt.close()


print('Done.')
