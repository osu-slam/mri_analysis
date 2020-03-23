%% YA_FST_MVPA.m
% Batch script to run all MVPA analysis on supercomputer. 

% MM/DD/YY -- CHANGELOG
% 
% 06/13/18 -- File initialized. Cloned to supercomputer. 
% 02/14/20 -- Cloned for YA_FST. Making scripts universal. MJH
% 02/18/20 -- Dropping physio and rerunning
% 02/20/20 -- New strategy to drop physio, rerunning. 
% 02/25/20 -- Model with all betas, wu- data, but only investigating TRs 2
%   through 4. 
% 03/05/20 -- New experiment structure. 

clearvars; clc; 

%% Pathing and parameters
dir_batch = pwd;
dir_MVPA  = fullfile(dir_batch, 'MVPA'); 
dir_process = fullfile(dir_batch, 'process'); 
dir_preprocess = fullfile(dir_batch, 'preprocess'); 
YA_FST_params

tic
for ii = 1:length(subj)
    %% Subject-specific parameters
    thissubj = subj(ii);
    disp(['MVPA on subject ' thissubj.name '...'])

    %% Specify and estimate unsmoothed GLM using FIR
%     disp(['Specifying 1st level GLM for subject ' thissubj.name '.'])
%     cd(dir_process); clear_spm_mat(thissubj, study, 5)
%     cd(dir_MVPA); spec_est_GLM_MVPA_v2(thissubj, study, 5, 0, 1)
%     disp('Done!')
    
    %% Extract time course information
%     disp(['Extracting time course information for ' thissubj.name '...'])
%     cd(dir_MVPA); extract_betas_v2(thissubj, study, 5, 'plot')
%     disp('Done!')

    %% Create coordinates spheres for all subjects
%     disp(['Creating spheres for ' thissubj.name '...'])
%     try parpool(28); catch; delete(gcp('nocreate')); parpool(28); end
%     cd(dir_MVPA); coord_sphere(thissubj, study, 5)
%     disp('Done!')
    
    %% Run the searchlight analysis
    disp(['Running searchlights on ' thissubj.name '...'])
    try parpool(28); catch; delete(gcp('nocreate')); parpool(28); end
% 
%     cd(dir_MVPA); searchlight_lang_clear_noise(thissubj, study, 5)
    cd(dir_MVPA); searchlight_lang_clear_noise_spatiotemporal(thissubj, study, 5)
%     cd(dir_MVPA); searchlight_lang_clear_noise_SVM(thissubj, study, 5)
%     cd(dir_MVPA); searchlight_or_sr_clear(thissubj, study, 5)
%     cd(dir_MVPA); searchlight_babble(thissubj, study, 5)
%     cd(dir_MVPA); searchlight_rate(thissubj, study, 5)
    disp('Done!')
    
    %% Make histograms for each subject!
    disp(['Making histograms for ' thissubj.name '...'])
    cd(dir_MVPA); plot_acc(thissubj, study, 5, 'GNB_lng_clear_noi_spatiotemporal') 
    disp('Done!')
    
end
toc

%% Second-level stats
% tic
% % cd(dir_batch_MVPA); second_level

% toc
% Subject3 OR > SR 

disp('Completed batch!')