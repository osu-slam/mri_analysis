%% BETA_extract_isss_multi
% Extracts betas for each subject, updated to work under MVPA_batch script. 
% 
% MM/DD/YY -- CHANGELOG
% 02/14/20 -- Cloned for YA_FST, made universal. MJH

function extract_betas_nophysio(subj, study, dd)
%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input (subj, study, dd) where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input (subj, study, dd) where study is struct')
end

if ~isnumeric(dd)
    error('Input (subj, study, dd) where dd specifies which design')
end

%% Pathing and parameters
dir_subj = fullfile(study.path, 'data', subj.name); 
dir_data_MVPA = fullfile(dir_subj, 'MVPA'); 
if ~isdir(dir_data_MVPA)
    mkdir(dir_data_MVPA) 
end

design = study.design(dd); 

dir_design = fullfile(dir_subj, 'design', [design.name '_nophysio']); 
SPM_file = fullfile(dir_design, 'SPM.mat');
mask_file = fullfile(dir_design, 'mask.nii'); 

%% Determine which voxels are inside the mask
mask_hdr = spm_vol(mask_file); % mask is ROI when applicable! 
mask_mat = spm_read_vols(mask_hdr);
mask_idx = find(mask_mat);     % The indices of where mask = 1
% mask_vox = length(mask_idx);
[mask_x, mask_y, mask_z] = ind2sub(size(mask_mat), mask_idx);
XYZ = [mask_x, mask_y, mask_z]';
% XYZ has three rows, and one col for every voxel in the mask

numCons = length(design.cond); 
numBetas = study.scan.epis * numCons * subj.runs; 
betas_all = nan(numBetas, length(XYZ));

%% Load design
load(SPM_file) 
labels = SPM.xX.name; 
path = [dir_design filesep];

%% Start calculating TC
beta_idx = 1;
idx_con = cell(numBetas, 1); 

for rr = 1:subj.runs 
    disp(['Run: ' num2str(rr) ' of ' num2str(subj.runs)]);

    for cc = 1:numCons
    % Select betas
        this_con_betas = find(contains(labels, ['Sn(' num2str(rr) ') ' design.cond{cc}]))';

        for tr = 1:study.scan.epis
            if this_con_betas(tr) < 10 % 1-9
                idx_con{beta_idx} = [path 'beta_000' num2str(this_con_betas(tr)) '.nii'];
            elseif this_con_betas(tr) < 100 % 10-99
                idx_con{beta_idx} = [path 'beta_00' num2str(this_con_betas(tr)) '.nii'];
            elseif this_con_betas(tr) < 1000 % 100-999
                idx_con{beta_idx} = [path 'beta_0' num2str(this_con_betas(tr)) '.nii'];
            else % 1000-9999
                idx_con{beta_idx} = [path 'beta_' num2str(this_con_betas(tr)) '.nii'];
            end

            betas_all(beta_idx, :) = spm_get_data(idx_con{beta_idx}, XYZ);
            beta_idx = beta_idx + 1;
            % betas_all(idx, 1) therefore is organized by the 
            % following logic, which is preserved through the code:
            % NOI1 -> SIL1 -> Each event run1 -> NOI2 -> SIL2 -> Each event -> etc...
        end

    end

end

% Manipulate entire matrix
disp('Zero meaning betas...')
betas_mean = mean(betas_all, 1);
betas_zero_mean = bsxfun(@minus, betas_all, betas_mean); %#ok<NASGU>

%% Sanity check
figure % a great sanity check!
histogram(betas_zero_mean)

%% Save data
disp('Saving data...')
filename = fullfile(dir_data_MVPA, [subj.name '_zero_meaned_betas_nophysio.mat']);
save(filename, 'betas_zero_mean', 'idx_con'); 

end