%% convert_dicm_physio
% Converts all physio dicm to nii, and put them in their proper places
%
% MM/DD/YY: CHANGELOG
% 02/07/20: File initialized, new physio practices at CCBBI
% 03/09/20: Updated to use subj/study input. Fixed errors regarding missing
%   folders. 

function convert_dicm_physio(subj, study)
%% Check input
if ~isstruct(subj)
    error('Input (subj) where subj is a SINGLE structure')
end

if ~isstruct(study)
    error('Input (study), which is a structure of study info!')
end

%% Path
thissubj  = subj.name; 
dir_subj  = fullfile(study.path, 'data', thissubj); 
dir_phys  = fullfile(dir_subj, 'PHYSIO');
dir_reg   = fullfile(dir_phys, 'reg');
dir_trace = fullfile(dir_phys, 'trace'); 

%% Convert physio dicms to nii
mkdir('dicoms')
mkdir('trace')
mkdir('reg')

target = fullfile(dir_phys, '*.dcm'); 
dcms = dir(target);

if isempty(dcms)
    error('There are no .dcm to extract in PHYSIO. Have you unzipped?')
end

if ~exist(dir_reg, 'file'); mkdir(dir_reg); end
if ~exist(dir_trace, 'file'); mkdir(dir_trace); end

disp('Converting physio dicms...')
for ii = 1:length(dcms)
    dcm = fullfile(dcms(ii).folder, dcms(ii).name);
    dicm2nii(dcm, dir_phys, '.nii 3D')
    
    movefile(dcm, 'dicoms')

    temp = fullfile(dir_phys, 'PulseRespiratoryRegressors.nii');
    if exist(temp, 'file')
        type = 'reg';
    else
        temp = fullfile(dir_phys, 'PulseRespiratoryTraces.nii');
        type = 'trace'; 
    end
    
    if ~exist(temp, 'file')
        error('No file found!')
    end
    
    
    if ii < 10
        final = fullfile(dir_phys, type, ['physio_' type '_0' num2str(ii) '.nii']);
    else
        final = fullfile(dir_phys, type, ['physio_' type '_' num2str(ii) '.nii']);
    end
    movefile(temp, final)

    temp = fullfile(dir_phys, 'dcmHeaders.mat');
    if ii < 10
        final = fullfile(dir_phys, type, ['header_' type '_0' num2str(ii) '.mat']);
    else
        final = fullfile(dir_phys, type, ['header_' type '_' num2str(ii) '.mat']);
    end
    movefile(temp, final)

end

end