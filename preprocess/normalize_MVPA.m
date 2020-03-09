%% normalize
% Normalize functional images to MNI space. 

% MM/DD/YY: CHANGELOG
% 02/11/20: Changelog started, file forked from isss_multi
% 02/17/20: Estimate, no reslice, and switching to 3x3x3.5mm voxel
% 02/17/20: Last change is an oxymoron? Est + Res, 3x3x3.5mm voxel. 
% 02/24/20: Forked for MVPA

function normalize_MVPA(subj, study)
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
dir_batch = pwd; 
batch = fullfile(dir_batch, 'matlabbatch', 'SPM_normalize.mat'); 
load(batch)

cd ..
cd ..
dir_spm = fullfile(pwd, 'spm12'); 
cd(dir_batch)

dir_subj = fullfile(study.path, 'data', subj.name); 

dir_anat = fullfile(dir_subj, 'ANATOMICAL'); 
dir_func = fullfile(dir_subj, 'FUNCTIONAL_v2_MVPA'); 
dir_thissubj_batch = fullfile(dir_subj, 'batch'); %% Pathing

%% Poke into matlabbatch
tpm = fullfile(dir_spm, 'tpm', 'TPM.nii'); 
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {tpm}; 
anat = fullfile(dir_anat, study.anat); 
matlabbatch{1}.spm.spatial.normalise.estwrite.subj(1).vol = {anat};

func_data = {}; 
for rr = 1:subj.runs
    [boldFiles, ~] = spm_select('List', dir_func, ['^u' study.runname num2str(rr) '_00*.*\.nii$']);
    boldFiles = [repmat([dir_func filesep], length(boldFiles),1), boldFiles];
    boldFiles = cellstr(boldFiles);
    func_data = vertcat(func_data, boldFiles); 
end

matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = func_data;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [3 3 3.5]; 

filename = fullfile(dir_thissubj_batch, 'normalize_z_3_5_MVPA.mat'); 
save(filename, 'matlabbatch')

spm_jobman('run',matlabbatch);

end