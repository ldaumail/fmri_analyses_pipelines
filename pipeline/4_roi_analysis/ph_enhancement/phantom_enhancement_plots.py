
import os
import os.path as op
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# get scan info from experiment file
subject = 'sub-F019'
session = '230321'
utils = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/code/utils'
sys.path.append(op.expanduser(f'{utils}'))
from experiment import experiment
#buffer = 3 # 3 volumes either side of epoch

# get conditions
contrasts = ['Inducers - (horizontal single right + single left) > 0']
condition = ['horDoubleIndir', 'horSingleLeft', 'horSingleRight', 'vertSingleLeft','vertSingleRight']

# First plot: baseline corrected timeseries
figSize = (4, 3)

for area in ['V1_lh','V1_rh', 'V2_lh', 'V2_rh', 'V3_lh', 'V3_rh', 'hV4_lh', 'hV4_rh', 'V3a_lh', 'V3a_rh', 'V3b_lh', 'V3b_rh']: #  ['V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b']:
	for z, zstat in enumerate(['zstat2', 'zstat3']):
		dataFile = f'/users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/topUp/noB0/psc/improved_roi/timeseries_{area}_{zstat}.csv'
		data = pd.read_csv(dataFile)
		# for condition in ['Inducers','Phantom']:

		# add columns for region, hemisphere and voxel counts
		regions, hemispheres, nVoxels, masks = [[], [], [], []]
		for row in data.index:
			regions.append("".join([data['region'][row].split(sep='_')[1],'_',data['region'][row].split(sep='_')[2]]))
			cont = data['region'][row].split(sep='_')[0]
			thr = data['region'][row].split(sep='_')[3]
			masks.append(f"{cont}{thr}")
			nVoxels.append(data['region'][row].split(sep='_')[4])
			# hemispheres.append(data['region'][row].split(sep='_')[1]) #use these lines if the mask has the name regions_hemisphere_nVoxels
			# nVoxels.append(data['region'][row].split(sep='_')[2])
		data['region'] = regions
		data['mask'] = masks
		#data['hemisphere'] = hemispheres
		data['nVoxels'] = nVoxels

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
							scan = 'phantom'
							params = experiment['design'][scan]['params']
							epochDur = int((params['blockDuration'] + params['IBI']) / params['TR'])
							buffer = 3  # 3 volumes either side of epoch
							xticks = np.arange(-buffer, epochDur + buffer + 1) * params[
								'TR']  # time points fall in between measurements
							x_pos = xticks[:-1] + 1  # points are shifted to centre on middle of TR

							outDir = os.path.join(f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/analysis/results', topup, b0, HRFmodel, subject, session, scan, 'improved_roi/baselined')
							os.makedirs(outDir, exist_ok=True)
							dataConf = data[(data['topup'] == topup) &
											(data['b0'] == b0) &
											(data['HRFmodel'] == HRFmodel) &
											(data['subject'] == subject) &
											(data['session'] == int(session)) &
											(data['scan'] == scan)].reset_index()

							# baseline
							#lastRow = 0
							PSCbaselined = [np.array([len(dataConf)])]
							bsDat = []
							for r, run in enumerate(dataConf['run'].unique()):
								dataRun = dataConf[dataConf['run'] == run]
								for start in range(0, len(dataRun), 20):
									rawData = dataRun['PSC'][start:start+20]
									bsDat = np.concatenate((bsDat, rawData[1:4].to_numpy()), axis=0)
								baseline = np.mean(bsDat)
								PSCbaselinedRun = dataRun['PSC']-baseline
								PSCbaselined[0+r*len(PSCbaselinedRun):(r+1)*len(PSCbaselinedRun)] = PSCbaselinedRun.to_numpy()
							dataConf = dataConf.assign(PSCbaselined=PSCbaselined)

							# make plots for each configuration of parameters
							for region in dataConf['region'].unique():
								dataRegion = dataConf[dataConf['region'] == region]
								for m, mask in enumerate(dataConf['mask'].unique()):
									numvoxels = dataRegion['nVoxels'][dataRegion['mask'] == mask].unique()
									numvoxels = numvoxels[0][9:]
									plt.figure(figsize=figSize)
									for s, stiLoc in enumerate(condition):
										means = np.empty([epochDur + buffer * 2])
										errorbar = np.empty([epochDur + buffer * 2])
										for timepoint in range(epochDur + buffer * 2):
											values_loc = np.array(dataRegion['PSCbaselined'][(dataRegion['condition']==stiLoc)&(dataRegion['TR'] == timepoint - buffer)&(dataRegion['mask'] == mask)])
											values_loc = np.delete(values_loc, np.isnan(values_loc))
											values_loc = values_loc.astype(float)
											means[timepoint] = np.mean(values_loc)
											errorbar[timepoint] = stats.sem(values_loc)
											# print(f'Time = {timepoint} | Stim: {stiLoc}')
										plt.errorbar(x_pos,
													 means,
													 yerr=errorbar,
													 marker='.',
													 markersize=5,
													 label = condition[s])
									plt.legend(bbox_to_anchor=(1.04, .75), borderaxespad=0)  # put legend outside plot
									plt.xticks(xticks)
									plt.xlabel('time (s)')
									plt.ylabel('signal change (%)')

									plt.title(f'timeseries: {region},Contrast: {contrasts[z]}, Mask threshold: {mask[-3:]}, nVoxels = {numvoxels}') #{size} voxels
									#plt.tight_layout(rect=[0, 0, 0.75, 1])  # ensure everything is placed in the canvas
									plt.subplots_adjust(right=1.4)  # allow space for legend
									plt.savefig(os.path.join(outDir, f'{region}_{mask}_nvox{numvoxels}_timeseries_baslined.png'), bbox_inches='tight')#{size}
									plt.show()
									plt.close()
print('Done.')


# Second plot: difference between response to pair and single
figSize = (4, 3)

for area in ['V1_lh','V1_rh', 'V2_lh', 'V2_rh', 'V3_lh', 'V3_rh', 'hV4_lh', 'hV4_rh', 'V3a_lh', 'V3a_rh', 'V3b_lh', 'V3b_rh']: #  ['V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b']:
	for z, zstat in enumerate(['zstat2']): # 'zstat3'
		dataFile = f'/users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/data/functional/phantom/allRuns/outputs/topUp/noB0/psc/phEnhancement_roi/timeseries_{area}_ph{zstat}.csv'
		data = pd.read_csv(dataFile)
		# for condition in ['Inducers','Phantom']:

		# add columns for region, hemisphere and voxel counts
		regions, hemispheres, nVoxels, masks = [[], [], [], []]
		for row in data.index:
			regions.append("".join([data['region'][row].split(sep='_')[1],'_',data['region'][row].split(sep='_')[2],'',data['region'][row].split(sep='_')[3]]))
			cont = data['region'][row].split(sep='_')[0]
			thr = data['region'][row].split(sep='_')[4]
			masks.append(f"{cont}{thr}")
			nVoxels.append(data['region'][row].split(sep='_')[5])
			# hemispheres.append(data['region'][row].split(sep='_')[1]) #use these lines if the mask has the name regions_hemisphere_nVoxels
			# nVoxels.append(data['region'][row].split(sep='_')[2])
		data['region'] = regions
		data['mask'] = masks
		#data['hemisphere'] = hemispheres
		data['nVoxels'] = nVoxels

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
							scan = 'phantom'
							params = experiment['design'][scan]['params']
							epochDur = int((params['blockDuration'] + params['IBI']) / params['TR'])
							buffer = 3  # 3 volumes either side of epoch
							xticks = np.arange(-buffer, epochDur + buffer + 1) * params[
								'TR']  # time points fall in between measurements
							x_pos = xticks[:-1] + 1  # points are shifted to centre on middle of TR

							outDir = os.path.join(f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/analysis/results', topup, b0, HRFmodel, subject, session, scan, 'phEnhancement_roi/baselined/difference')
							os.makedirs(outDir, exist_ok=True)
							dataConf = data[(data['topup'] == topup) &
											(data['b0'] == b0) &
											(data['HRFmodel'] == HRFmodel) &
											(data['subject'] == subject) &
											(data['session'] == int(session)) &
											(data['scan'] == scan)].reset_index()

							# baseline
							#lastRow = 0
							PSCbaselined = [np.array([len(dataConf)])]
							bsDat = []
							for r, run in enumerate(dataConf['run'].unique()):
								dataRun = dataConf[dataConf['run'] == run]
								for start in range(0, len(dataRun), 20):
									rawData = dataRun['PSC'][start:start+20]
									bsDat = np.concatenate((bsDat, rawData[1:4].to_numpy()), axis=0)
								baseline = np.mean(bsDat)
								PSCbaselinedRun = dataRun['PSC']-baseline
								PSCbaselined[0+r*len(PSCbaselinedRun):(r+1)*len(PSCbaselinedRun)] = PSCbaselinedRun.to_numpy()
							dataConf = dataConf.assign(PSCbaselined=PSCbaselined)

							# make plots for each configuration of parameters
							for region in dataConf['region'].unique():
								dataRegion = dataConf[dataConf['region'] == region]
								for m, mask in enumerate(dataConf['mask'].unique()):
									numvoxels = dataRegion['nVoxels'][dataRegion['mask'] == mask].unique()
									numvoxels = numvoxels[0][9:]
									plt.figure(figsize=figSize)
									# for s, stiLoc in enumerate(condition):
									pairVals = np.empty([epochDur + buffer * 2, len(dataConf['run'].unique())*len(dataConf['rep'].unique())])
									# pairErrorbar = np.empty([epochDur + buffer * 2])
									for timepoint in range(epochDur + buffer * 2):
										values_loc = np.array(dataRegion['PSCbaselined'][(dataRegion['condition'] == 'horDoubleIndir')&(dataRegion['TR'] == timepoint - buffer)&(dataRegion['mask'] == mask)])
										# values_loc = np.delete(values_loc, np.isnan(values_loc))
										values_loc = values_loc.astype(float)
										pairVals[timepoint,:] = values_loc
										# pairMeans[timepoint] = np.mean(values_loc)
										# pairErrorbar[timepoint] = stats.sem(values_loc)
										# print(f'Time = {timepoint} | Stim: {stiLoc}')

									singleVals = np.empty([epochDur + buffer * 2, len(dataConf['run'].unique())*len(dataConf['rep'].unique())])
									# singleErrorbar = np.empty([epochDur + buffer * 2])
									for timepoint in range(epochDur + buffer * 2):
										values_left = np.array(dataRegion['PSCbaselined'][(dataRegion['condition'] == 'horSingleRight') & (dataRegion['TR'] == timepoint - buffer) & (dataRegion['mask'] == mask)])
										values_right = np.array(dataRegion['PSCbaselined'][(dataRegion['condition'] == 'horSingleLeft') & (dataRegion['TR'] == timepoint - buffer) & (dataRegion['mask'] == mask)])
										# values_loc = np.delete(values_loc, np.isnan(values_loc))
										values_left = values_left.astype(float)
										values_right = values_right.astype(float)
										singleVals[timepoint, :] = np.add(values_left, values_right)
										# singleMeans[timepoint] = np.mean(values_loc)
										# singleErrorbar[timepoint] = stats.sem(values_loc)
									# print(f'Time = {timepoint} | Stim: {stiLoc}')
									diff_pair_vs_single = np.subtract(pairVals, singleVals)
									mean_diff = np.mean(diff_pair_vs_single, axis= 1)
									errorbar = stats.sem(diff_pair_vs_single, axis= 1)
									plt.errorbar(x_pos,
												 mean_diff,
												 yerr=errorbar,
												 marker='.',
												 markersize=5,
												 label = 'Mean difference: pair - (single left + single right)')

									plt.legend(bbox_to_anchor=(1.04, .75), borderaxespad=0)  # put legend outside plot
									plt.xticks(xticks)
									plt.xlabel('time (s)')
									plt.ylabel('signal change (%)')

									plt.title(f'timeseries: {region},Contrast: {contrasts[z]}, Mask threshold: {mask[-3:]}, nVoxels = {numvoxels}') #{size} voxels
									#plt.tight_layout(rect=[0, 0, 0.75, 1])  # ensure everything is placed in the canvas
									plt.subplots_adjust(right=1.4)  # allow space for legend
									plt.savefig(os.path.join(outDir, f'{region}_{mask}_nvox{numvoxels}_timeseries_baselined.png'), bbox_inches='tight')#{size}
									plt.show()
									plt.close()
print('Done.')
