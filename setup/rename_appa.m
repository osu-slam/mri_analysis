%% rename_appa
% Renames AP-PA files to reflect resting state and task-based maps. Lifted
% mostly from rename_files. 

% MM/DD/YY -- CHANGELOG
% 03/09/20 -- File initialized. 

function rename_appa(subj, study)
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
dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 

dir_nii  = fullfile(dir_subj, 'nii');
dir_dicm = fullfile(dir_subj, 'dicm');
dir_anat = fullfile(dir_subj, 'ANATOMICAL');
dir_real = fullfile(dir_subj, 'realign');

target = subj.appa; 

%% Rename files
for tt = 1:length(target)
    % Get sessions
    temp = fullfile(dir_real, [target{tt} '*']); 
    files = dir(temp); files = fullfile(dir_real, {files(:).name}); 
    split = cellfun((@(x) strsplit(x, '_')), files, 'UniformOutput', false); 
    imgnum = cellfun((@(x) x{end}), split, 'UniformOutput', false); 
    
    % Build new name
    type = target{tt}(1:2); 
    newname = cell(length(imgnum), 1);
    for ii = 1:length(newname)
        newname{ii} = fullfile(dir_real, [type '_' imgnum{ii}]); 
    end
    
    for ff = 1:length(files)
        movefile(files{ff}, newname{ff}); 
    end
    
end

end
