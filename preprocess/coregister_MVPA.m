%% coregister
% Coregister structural and functional images

% MM/DD/YY: CHANGELOG
% 02/11/20: Forked from isss_multi, turned into universal function. 
% 02/21/20: New dataset with fewer realigned images. Updated!
% 02/24/20: Forked for MVPA

function coregister_MVPA(subj, study)
%% Check inputs
if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_figure('GetWin','Graphics'); % Thanks Guillaume
disp('Done!')

%% Pathing
cd ..
batch = fullfile(pwd, 'matlabbatch', 'SPM_coregister.mat'); 
load(batch)

dir_subj = fullfile(study.path, 'data', subj.name); 

dir_anat = fullfile(dir_subj, 'ANATOMICAL'); 
dir_func = fullfile(dir_subj, 'FUNCTIONAL_v2_MVPA'); 
dir_thissubj_batch = fullfile(dir_subj, 'batch'); 

%% Poke into matlabbatch
anat = fullfile(dir_anat, study.anat); 
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(anat);

meanimg = ['meanu' study.runname '1_00007.nii']; 
reference = fullfile(dir_func, meanimg); 
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(reference);

batchname = fullfile(dir_thissubj_batch, 'coregister_meanu_MVPA.mat');
save(batchname, 'matlabbatch')

spm_jobman('run',matlabbatch);    
    
end