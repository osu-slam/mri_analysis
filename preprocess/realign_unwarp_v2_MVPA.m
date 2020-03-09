%% realign_unwarp_v2
% Realigns and unwarps functional images using SPM. 
% 
% MM/DD/YY: Changelog
% 02/10/20: Forked from isss_multi. 
% 02/21/20: Version 2, drops the first five images
% 02/24/20: Forked for MVPA, where we drop the first TR of each event. 

function realign_unwarp_v2_MVPA(subj, study)
%% Check input
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

dir_subj = fullfile(study.path, 'data', subj.name);

dir_reg  = fullfile(dir_subj, 'reg');
dir_func = fullfile(dir_subj, 'FUNCTIONAL_v2_MVPA'); 
dir_realign  = fullfile(dir_subj, 'realign');
dir_motionps = fullfile(dir_subj, 'ps', 'preprocessing');
dir_thissubj_batch = fullfile(dir_subj, 'batch'); 

regNames = cell(1, subj.runs);

%% Drop the first five images
% The first image is a dummy scan. The next four are not associated with
% any stimuli. 
target = fullfile(dir_func, '*.nii'); 
funcfiles = dir(target); funcfiles = {funcfiles(:).name}';
deletethese = {'00001', '00002', '00003', '00004', '00005'}; 

these = false(length(funcfiles), 1);
for ii = 1:length(deletethese)
    these = these | contains(funcfiles, deletethese{ii}); 
end

if ~any(these) 
    warning('Not deleting any extra EPIs!')
else
    files = funcfiles(these)';
    files = fullfile(dir_func, files);
    for ii = 1:length(files)
        delete(files{ii})
    end
    
end

%% Drop the first TR of each event
target = fullfile(dir_func, '*.nii'); 
funcfiles = dir(target); funcfiles = {funcfiles(:).name}';

numscans = strsplit(funcfiles{end}, '.'); numscans = numscans{1}; 
numscans = str2double(numscans(end-1:end));
targets = 6:study.scan.epis:numscans; 
deletethese = cell(length(targets), 1); 
for ii = 1:length(targets)
     if targets(ii) < 10
         deletethese{ii} = ['0000' num2str(targets(ii))];
     elseif targets(ii) < 100
         deletethese{ii} = ['000' num2str(targets(ii))];
     elseif targets(ii) < 1000
         deletethese{ii} = ['00' num2str(targets(ii))];
     else
         error('Too large a number!')
     end
     
end

these = false(length(funcfiles), 1);
for ii = 1:length(deletethese)
    these = these | contains(funcfiles, deletethese{ii}); 
end

if ~any(these) 
    warning('Not deleting any extra EPIs!')
else
    files = funcfiles(these)';
    files = fullfile(dir_func, files);
    for ii = 1:length(files)
        delete(files{ii})
    end
    
end

for rr = 1:subj.runs
    %% Load names and put into mattlabbatch
    [boldFiles, ~] = spm_select('List', dir_func, ['^' study.runname num2str(rr) '_00*.*\.nii$']);
    boldFiles = [repmat([dir_func filesep], length(boldFiles),1), boldFiles];
    boldFiles = cellstr(boldFiles);
    matlabbatch{1}.spm.spatial.realignunwarp.data(rr).scans = boldFiles;
    
    [vdmFiles, ~] = spm_select('List', dir_realign, '^vdm5_*.*\.nii$');
    vdmFiles = cellstr([dir_realign filesep vdmFiles]);
    matlabbatch{1}.spm.spatial.realignunwarp.data(rr).pmscan = vdmFiles;
    
    if strcmp(matlabbatch{1}.spm.spatial.realignunwarp.data(rr).scans, '<UNDEFINED>')
        matlabbatch{1}.spm.spatial.realignunwarp.data = [];
    end
    
end

filename = fullfile(dir_thissubj_batch, 'realign_unwarp_v2_MVPA.mat');
save(filename, 'matlabbatch') 

disp('Starting realign_unwarp...')
spm_jobman('run',matlabbatch);

%% Move regressors and motion graph to their proper directory
target = fullfile(dir_func, '*.txt'); 
regs = dir(target);
regNames = fullfile(dir_func, {regs(:).name}');
for ii = 1:length(regNames)
    movefile(regNames{ii}, dir_reg);
end

target = fullfile(pwd, '*.png');
psfiles = dir(target); 
psNames = fullfile(pwd, {psfiles(:).name}'); 
for ii = 1:length(psfiles)
    target = fullfile(dir_motionps, ['realign_unwarp_v2_' num2str(ii) '_MVPA.png']); 
    movefile(psNames{ii}, target)
end

end