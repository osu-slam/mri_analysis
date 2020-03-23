%% YA_FST_process
% Builds models for each subject
% CHANGELOG (MM/DD/YY)
% 11/13/17  Initialized file -- MH
% 02/28/17  Running an analysis to see which acquisition window is best
% 05/30/18  Supercomputer analysis. Forking from original. 
% 05/15/19  Modeling each event as its own beta for MVPA. 
% 12/23/19  Revisiting this monster over the holidays so we can avoid the
% sinister timing error. The idea is as follows:
%   1. Model each event (including wrong ones)
%   2. Create AUE using only the correct events
%   3. Average first within runs (correct trials) and then across runs. 
%   4. Second level SPA analysis as usual. 
% 02/12/20  Forked for YA_FST. Another timing error found, should be
%   resolved with how I calculated onsets. Adding universal functions.
% 02/17/20  Modifying 1st level design so it specifies design with 1 run
%   for visual inspection. Re-preprocessed with 3 x 3 x 3.5mm voxels. 
% 02/18/20  Dropping physio.
% 02/21/20  New timing scheme, dropped first 5 images
% 03/10/20  New design with CONN regressors

clc; clearvars 

%% Pathing
dir_batch = pwd;
dir_process = fullfile(dir_batch, 'process');
YA_FST_params

for ii = 2:length(subj)
    thissubj = subj(ii);
    disp(['Starting to process ' thissubj.name '.'])
    
    %% Extract timing
%     disp('Extracting timing from behav for lang...')
% %     cd(dir_process); extract_timing_nowrong_v3(subj(ii), study);
%     cd(dir_process); extract_timing_all(thissubj, study); 
%     disp('Done!')
        
    %% Specify and estimate GLM using FIR
    disp(['Specifying 1st level GLM for subject ' thissubj.name '.'])
    cd(dir_process); clear_spm_mat(thissubj, study, 6)
    cd(dir_process); spec_est_GLM(thissubj, study, 6, 0, 0)
    % Input (thissubj, study, design#, do first run, do interactive)
    disp('Done!')
    
    %% Build contrasts
    disp(['Building contrasts for subject ' thissubj.name '.'])
    cd(dir_process); build_contrasts(thissubj, study, 6)
    disp('Done!')
    
    %% Calculate AUE (stopped here!)
%     cd(dir_process)
%     disp('Calculating AUE...')
%     SPA_calculate
%     cd(dir_process)
%     SPA_manipulate
%     disp('Done!')  

    %% Calculate AUE for windows
%     cd(dir_process)
%     SPA_calculate_window
%     disp('Done!')

    %% Window analysis
%     cd(dir_process)
%     disp(['Building window contrasts for subject ' thissubj '.'])
%     build_contrasts_window
%     disp('Done!')  

    %% Check results
%     disp('Checking results...')
%     cd(dir_process)
%     results_report_swau
%     disp('Done!')    
    
end

%% Build second-level GLMs
% disp('Building, evaluating, and generating results for second-level GLMs...')
% cd(dir_process)
% second_level_AUE_and_results
% disp('Done!')

%% Run ROI analysis
% disp('ROI analysis...')
% cd(dir_roi)
% beta_plot_ROI_sphere
% disp('Done!')






