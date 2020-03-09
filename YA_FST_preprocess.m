%% YA_FST_preprocess.m
% Preprocesses images for analysis. 
% 
% MM/DD/YY: Changelog
% 02/10/20: Forked from isss_multi. 
%   + Fieldmaps are created using FSL (batch/unwarp.sh) and can be run from
%     here or from bash. 
% 03/09/20: Subjects 5976YL and 5977YL ran without issues. 

clc; clearvars

%% Pathing
dir_batch = pwd;
dir_preprocess = fullfile(dir_batch, 'preprocess'); 
dir_snr = fullfile(dir_batch, 'Noise_script'); 
YA_FST_params

cd ..
dir_data = fullfile(pwd, 'data'); 

tic
for ii = 2:length(subj)
    %% Setup for this subject
    thissubj = subj(ii);
    dir_subj = fullfile(dir_data, thissubj.name); 
    disp(['Preprocessing subj ' thissubj.name '...'])
    
    %% Unwarp and realignment
    % Inspect FM and convert to VDM manually!
    disp('Unwarping and realigning data')
    cd(dir_preprocess); realign_unwarp_v3(thissubj, study)
    disp('Done!')

    %% Coregistration
    disp('Coregistration begins')
    cd(dir_preprocess); coregister(thissubj, study)
    disp('Done coregistering!')
    
    %% Normalization
    disp('Normalizing to MNI-space')
    cd(dir_preprocess); normalize(thissubj, study)
    disp('Done normalizing!')

    %% Smoothing
    disp('Smoothing data now')
    cd(dir_preprocess); smooth(thissubj, study)
    disp('Done smoothing!')

    %% SNR
    disp('Running SNR scripts')
    cd(dir_snr); snr_sd_v4(thissubj, study)
    disp('Done with SNR!')

end

disp('Batch complete!')
toc
