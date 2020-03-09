#!/bin/bash
# Loads FSL and creates fieldmap. 
# Input: (1) Subject number

# 0. Load FSL and set paths
module load fsl
source $FSLDIR/etc/fslconf/fsl.sh

cd ..
cd data
data_in=$PWD/topup_datain.txt

cd $1
dir_subj=$PWD

# 1. Use fslmerge to combine EPIs
cd realign
fslmerge -t merge.nii.gz AP_00004.nii PA_00004.nii
merged_img=$PWD/merge.nii.gz

# 2. Run topup
topup --imain=$merged_img --datain=$data_in --config=b02b0.cnf --fout=fpm_fieldmap.nii.gz --iout=merge_unwarp.nii.gz 

# 3. Convert from nii.gz to nii
gunzip fpm_fieldmap.nii.gz

