%% unzip_data
% Unzips all scans and physio regressors
% 
% MM/DD/YY: CHANGELOG
% 02/07/20: Changelog started. Forked from isss_multi for YA_FST. MJH
% 03/09/20: Updated to use subj/study input. 

function unzip_data(subj, study)
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

%% Do the thing
try
    cd(dir_subj)
catch err
    disp(['Could not change to subj dir ', dir_subj, '. Does the folder exist?'])
    rethrow(err)
end

cd('zip')
files = dir('*.zip');
if length(files) == 1
    disp('Found one file, no physio present?')
elseif length(files) == 2
    disp('Found two files, images and physio I hope!')
else
    error('Strange number of files!')
end

filenames = {files(:).name};
physio = contains(filenames, 'physio'); 
dicom = ~physio; 

dir_dicm = fullfile(dir_subj, 'dicm');
disp('Unzipping dicm zip now...')
unzip(filenames{dicom}, dir_dicm);
disp('Done!')

if any(physio)
    dir_physio = fullfile(dir_subj, 'PHYSIO');
    disp('Unzipping physio zip now...')
    unzip(filenames{physio}, dir_physio);
    disp('Done!')
end

end