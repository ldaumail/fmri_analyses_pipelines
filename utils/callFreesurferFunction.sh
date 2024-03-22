#!/bin/bash

while getopts s: option
do
case "${option}"
in
s) STRING=${OPTARG};;
esac
done

FSLDIR=/Users/tonglab/fsl
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/share/fsl/bin:${PATH}
export FSLDIR PATH
export FREESURFER_HOME=/Applications/freesurfer/7.4.1
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=$HOME/Documents/Loic/all_mri_data/freesurfer

$STRING

#export FREESURFER_HOME=/Applications/freesurfer/7.4.1
#export FSFAST_HOME=/Applications/freesurfer/7.4.1/fsfast
#export FSF_OUTPUT_FORMAT=nii.gz
#export SUBJECTS_DIR=$HOME/Documents/Loic/all_mri_data/freesurfer
#export MNI_DIR=/Applications/freesurfer/7.4.1/mni
#export FSL_DIR=/Applications/fsl


