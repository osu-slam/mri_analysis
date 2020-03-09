%% convert_dicm
% Converts all image dicm to nii, and put them in their proper places
%
% MM/DD/YY: CHANGELOG
% 02/07/20: Changelog started, forked from isss_multi, now universal!
% 03/09/20: Updated to use subj/study input. 

function convert_dicm(subj, study)
%% Check input
if ~isstruct(subj)
    error('Input (subj) where subj is a SINGLE structure')
end

if ~isstruct(study)
    error('Input (study), which is a structure of study info!')
end

%% Path
thissubj = subj.name; 
dir_subj = fullfile(study.path, 'data', thissubj); 
dir_nii  = fullfile(dir_subj, 'nii');
dir_dicm = fullfile(dir_subj, 'dicm');
dir_anat = fullfile(dir_subj, 'ANATOMICAL');
dir_func = fullfile(dir_subj, 'FUNCTIONAL');
dir_real = fullfile(dir_subj, 'realign');

%% Check if dicm2nii is on path
if ~contains(path, 'dicm2nii')
    disp('Adding dicm2nii to path!')
    cd ../../..; addpath dicm2nii
    if ~contains(path, 'dicm2nii')
        error('Adding dicm2nii to path failed.')
    end
    
end

%% Convert dicms to nii for scan
disp(['Converting scan dicms for ' dir_subj '...'])
dicm2nii(dir_dicm, dir_nii, '.nii 3D')
disp('Done!')

%% Move niis to proper locations
cd(dir_nii)
mprage = dir(study.anat);
AP     = dir('AP_*.nii'); 
PA     = dir('PA_*.nii'); 
APPA   = [AP; PA]; 
epi    = dir([study.runname '*.nii']);

disp(['Moving scans into proper folders for ' dir_subj])
movefile(mprage.name, dir_anat);

for ii = 1:length(epi)
    movefile(epi(ii).name, dir_func);
end

for ii = 1:length(APPA)
    movefile(APPA(ii).name, dir_real);
end

disp('Done!')

end