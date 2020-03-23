%% YA_FST_conn_batch
% Uses CONN to perform aCompCorr, creates regressors for WM and CSF signal
% fluctuations. These regressors are used to remove noise from GLMs. 
% Author: Matthew Heard

% MM/DD/YY -- CHANGELOG
% 03/16/20 -- Changelog started. 

% Sample CONN Script for the OSU Workshop
% Created by Andrew Jahn, University of Michigan, 02.27.2020
% Adapted from Alfonso Nieto-Castanon's script, www.alfnie.com 

%% Parameters and path
YA_FST_params

NSUBJECTS = length(subj);

%% Set up batch
clear batch
batch.filename=fullfile(pwd,'CONN_physio_v3.mat'); 

batch.Setup.isnew = 1;
batch.Setup.nsubjects = NSUBJECTS;
batch.Setup.RT = study.scan.TR; 
batch.Setup.functionals = repmat({{}}, [NSUBJECTS,1]);  

for nsub = 1:length(subj)
    thissubj = subj(nsub); 
    %% Paths for each subject
    dir_subj = fullfile(study.path, 'data', thissubj.name); 
    
    dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 
    dir_anat = fullfile(dir_subj, 'ANATOMICAL'); 
    
    %% Parameters
    nsessions = thissubj.runs - length(thissubj.drop); 
    
    %% Functional and structural data for each subject
    for nses = 1:nsessions
        target = fullfile(dir_func, ... 
            [study.runname num2str(nses) '*.nii']); % un-preprocessed!
        files_func = dir(target); 
        files_func = fullfile(dir_func, {files_func(:).name})';
    
        batch.Setup.functionals{nsub}{nses}{1} = files_func; 
    end 
    
    target = fullfile(dir_anat, study.anat); % we may need preprocessed one?
    % May need to preprocess anatomical images in this script...
    files_anat = dir(target); 
    files_anat = fullfile(dir_anat, files_anat(:).name); 
    
    batch.Setup.structurals{nsub, 1} = files_anat;
    
    %% Load condition/session info
%     nconditions = nsessions; 
%     batch.Setup.conditions.names = [{'task'}, ... 
%         arrayfun(@(n) sprintf('Session%d', n), 1:nconditions, 'uni', 0)];

    nconditions = 1; 
    batch.Setup.conditions.names = {'task'}; 
    for ncond = 1
        for nses = 1:nsessions
            batch.Setup.conditions.onsets{ncond}{nsub}{nses} = 0;
            batch.Setup.conditions.durations{ncond}{nsub}{nses} = inf;
        end

    end
    
%     for ncond = 1:nconditions
%         for nses = 1:nsessions
%             batch.Setup.conditions.onsets{1+ncond}{nsub}{nses} = [];
%             batch.Setup.conditions.durations{1+ncond}{nsub}{nses} = []; 
%         end
% 
%     end
    
end

% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
% Prepares batch structure

     % Point to functional volumes for each subject/session

batch.Setup.preprocessing.steps='default_mni';
batch.Setup.preprocessing.sliceorder='interleaved (Siemens)';
batch.Setup.done=1;
batch.Setup.overwrite='Yes';

%Uncomment the following 2 lines if you want to use Andy's custom atlas
%batch.Setup.rois.files{1}='ROIs/AndyROIs.nii';
%batch.Setup.rois.multiplelabels = 1;

% uncomment the following 3 lines if you prefer to run one step at a time:
% conn_batch(batch); % runs Preprocessing and Setup steps only
% clear batch;
% batch.filename=fullfile(cwd,'Arithmetic_Scripted.mat');            % Existing conn_*.mat experiment name

% DENOISING step
% CONN Denoising                                    % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options 
batch.Denoising.filter=[0.01, 0.1];                 % frequency filter (band-pass values, in Hz)
batch.Denoising.done=1;
batch.Denoising.overwrite='Yes';

% uncomment the following 3 lines if you prefer to run one step at a time:
% conn_batch(batch); % runs Denoising step only
% clear batch;
% batch.filename=fullfile(cwd,'Arithmetic_Scripted.mat');            % Existing conn_*.mat experiment name

% FIRST-LEVEL ANALYSIS step
% CONN Analysis                                     % Default options (uses all ROIs in conn/rois/ as connectivity sources); see conn_batch for additional options 
batch.Analysis.done=1;
batch.Analysis.overwrite='Yes';

% Run all analyses
conn_batch(batch);

% CONN Display
% launches conn gui to explore results
conn
conn('load',fullfile(pwd,'Arithmetic_Scripted.mat'));
conn gui_results