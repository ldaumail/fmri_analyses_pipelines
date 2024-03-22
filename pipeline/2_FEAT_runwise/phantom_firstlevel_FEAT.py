#!/usr/bin/python
'''
takes first level design files (each run) made prior to this script using feat_gui, edits and submits to feat
Requirement: create masks using the localizer task
'''

import os
import os.path as op
import glob
import pickle
import datetime
import shutil
import time
#time.sleep(10000)
# get scan info from experiment file
# get scan info from experiment file
subject = 'sub-F019'
session = '230321'
utils = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}/code/utils'
sys.path.append(op.expanduser(f'{utils}'))
from experiment import experiment
# Preprocessing with FEAT
subjanaldir = f'/Users/tonglab/Documents/Loic/phantom_mri_data/analysis/loic_anal/{subject}_{session}'
dataDir = f'{subjanaldir}/data/sourceData'

scan = 'phantom'
fieldStrength = list(experiment['scanInfo'].keys())[0]
timeseries = []
session = list(experiment['scanInfo'][fieldStrength][subject].keys())[0]
nRuns = len(experiment['scanInfo'][fieldStrength][subject][session]['funcScans'][scan])

for r in range(0,nRuns):

    outDir = f"{subjanaldir}/data/functional/{scan}/run{r+1:02}/outputs/topUp/noB0/preprocessing.feat"
    os.makedirs(f"{subjanaldir}/data/functional/{scan}/run{r+1:02}/outputs/topUp/noB0", exist_ok=True)
    #if os.path.isdir(os.path.dirname(outDir)): # dont bother checking for eg with topup dirs for 3T scans

    print(f'\n{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Subject: {subject} | fieldStrength: {fieldStrength} | Session: {session} | Scan: {scan} | Run: {r + 1} | TopUp')

    #if not os.path.isdir(outDir):
    print('Analysis not found / just deleted, analysing...')

    # replace relevant settings in design file. The original should be set up for the first run of the first subject,
    # and the files should be set up so that all that needs changing is the subject ID and the run number
    designFile = f'{subjanaldir}/designs/task-{scan}/first_level/design.fsf'
    with open(designFile, 'r') as file:
        fileData = file.read()

    # replace topup and b0 types for input field and output field
    # fileData = fileData.replace('rawData', f'rawData{topupString}{b0String}')
    # fileData = fileData.replace('noTopUp', topup)
    # fileData = fileData.replace('noB0', b0)

    # replace field strength
    # fileData = fileData.replace('3T', fieldStrength)

    # replace subject ID
    # fileData = fileData.replace('M015', subject)

    # replace session
    # fileData = fileData.replace('190612', session)

    # replace run number
    fileData = fileData.replace('run01', f'run{r+1:02}')
    fileData = fileData.replace('run-001', f'run-{r + 1:03}')

    # write the file out again
    designFileTemp = f'{designFile[:-4]}_temp.fsf'
    with open(designFileTemp, 'w') as file:
        file.write(fileData)
    file.close()

    # check out dir in design file is correct
    lines = fileData.splitlines()
    for line in lines:
        if line.startswith(f'set fmri(outputdir)'):
            actualOutDir = line.split()[2][1:-1]
            targetedOutDir = os.path.join(os.getcwd(), outDir)
            if not actualOutDir == targetedOutDir:
                print(targetedOutDir)
                print(actualOutDir)
                raise Exception(
                    f'Destination directory does not match, check feat analysis at {outDir}')

    # run analysis
    os.system(f'feat {designFileTemp}')

    # else:
    #     print('Analysis found and overwrite = False, skipping...')
