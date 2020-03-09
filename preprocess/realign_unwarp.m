%% realign_unwarp
% Realigns and unwarps functional images using SPM. 
% 
% MM/DD/YY: Changelog
% 02/10/20: Forked from isss_multi. 

function realign_unwarp(dir_subj, subj, study)
%% Check input
if ~ischar(dir_subj)
    error('Input ("dir_subj") where dir_subj is a string')
end

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
batch = fullfile(pwd, 'matlabbatch', 'SPM_realign_unwarp.mat'); 
load(batch)

dir_realign = fullfile(dir_subj, 'realign');
dir_reg = fullfile(dir_subj, 'reg');
dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 
dir_thissubj_batch = fullfile(dir_subj, 'batch'); 
dir_motionps = fullfile(dir_subj, 'ps', 'preprocessing');

regNames = cell(1, subj.runs);

for rr = 1:subj.runs
    %% Load names and put into mattlabbatch
    [boldFiles, ~] = spm_select('List', dir_func, ['^' study.runname num2str(rr) '_00*.*\.nii$']);
    boldFiles = [repmat([dir_func filesep], length(boldFiles),1), boldFiles];
    boldFiles = cellstr(boldFiles);
    matlabbatch{1}.spm.spatial.realignunwarp.data(rr).scans = boldFiles;
    
    [vdmFiles, ~] = spm_select('List', dir_realign, '^vdm5_*.*\.nii$');
    vdmFiles = cellstr([dir_realign filesep vdmFiles]);
    matlabbatch{1}.spm.spatial.realignunwarp.data(rr).pmscan = vdmFiles;
    
    regNames{rr} = fullfile(dir_func, ['rp_' study.runname num2str(rr) '_00001.txt']);
    
    if strcmp(matlabbatch{1}.spm.spatial.realignunwarp.data(rr).scans, '<UNDEFINED>')
        matlabbatch{1}.spm.spatial.realignunwarp.data = [];
    end
    
end

filename = fullfile(dir_thissubj_batch, 'realign_unwarp.mat');
save(filename, 'matlabbatch') 

disp('Starting realign_unwarp...')
spm_jobman('run',matlabbatch);

%% Move regressors and motion graph to their proper directory
for ii = 1:length(regNames)
    movefile(regNames{ii}, dir_reg);
end

psfiles = dir('*.pdf'); 
for ii = 1:length(psfiles)
    target = fullfile(dir_motionps, ['realign_unwarp_' num2str(ii) '.pdf']); 
    movefile(psfiles(ii).name, target)
end

end