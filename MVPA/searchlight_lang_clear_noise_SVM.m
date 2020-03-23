%% searchlight_or_sr
% Does a binary classification using searchlight-based GNB.
% Adapted from Yune's code on 12/14/2017
% Updated for SLAM_isss_multi_MVPA_batch on 6/13/18
% Made parallel on 6/22/18
% Split for noise/silence and lng/noi on 6/26
% 
% MM/DD/YY -- CHANGELOG
% 02/14/20 -- Cloned for YA_FST. Made into function. 
% 02/17/20 -- Forked for language (clear) vs noise classification. 
% 02/18/20 -- Found error in HIT/FA logic! Fixed. 
% 02/25/20 -- New design. Does it work?
% 03/10/20 -- Cloned for SVM

function searchlight_lang_clear_noise_SVM(subj, study, dd)
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

%% Pathing and some parameters
if ~contains(path, 'libsvm-3.24')
    disp('Adding libSVM to path!')
    home = pwd; 
    cd ../../../libsvm-3.24/matlab
    addpath(pwd); 
    cd(home)
end

radius = study.mvpa.radius;

dir_subj = fullfile(study.path, 'data', subj.name); 
dir_data_MVPA = fullfile(dir_subj, 'MVPA'); 

design = study.design(dd); 
dir_design = fullfile(dir_subj, 'design', design.name); 
mask_file = fullfile(dir_design, 'mask.nii'); 
SPM_file = fullfile(dir_design, 'SPM.mat'); 

%% Load mask and SPM.mat
Vmask = spm_vol(mask_file);
mask_mat = spm_read_vols(Vmask);
mask_idx = find(mask_mat); 
num_mask_voxels = length(mask_idx);

load(SPM_file)

%% Load TC and spheres
filename = fullfile(dir_data_MVPA, ...
    [subj.name '_' design.name '_betas.mat']);
load(filename)

filename = fullfile(dir_data_MVPA, ...
    [subj.name '_' design.name '_spheres_radius' num2str(radius) '.mat']);  
load(filename);  

%% More parameters and preallocation
numvox = length(sphere_XYZ_indices_cell);  %#ok<USENS>

crossval_acc_mat_folds = zeros(subj.runs, numvox);

blocks = 1:subj.runs-length(subj.drop); 

%% Read design matrix
% Read the column-names from the design matrix, to figure out which
% regressors to use for defining our time-points
theselabels = SPM.xX.name;
regressors = ~cellfun(@isempty, regexp(theselabels, 'R\d'));
constant = contains(theselabels, 'constant');
removethese = regressors | constant;
theselabels(removethese) = [];

for fold = blocks % for each fold... 
    disp(['Cross-validation fold ' num2str(fold)]);
    
    %% Preallocate 
    crossval_acc_vector = zeros(1, numvox);
    
    %% making training & testing sets 
    test = fold; 
    tsttag = {['Sn(' num2str(test) ')']}; 
    
    block_tst_set = contains(theselabels, tsttag{1});
    
    lang_set = contains(theselabels, 'clear'); 
    lang_tst_set = block_tst_set & lang_set; 
    
    noi_set = contains(theselabels, 'NOI');
    noi_tst_set = block_tst_set & noi_set; 
    
    train = blocks; train(fold) = []; 
    trntag = cell(length(train), 1); 
    block_trn_set = false(size(block_tst_set)); 
    for ii = 1:length(train)
        trntag{ii} = ['Sn(' num2str(train(ii)) ')']; 
        block_trn_set = block_trn_set | contains(theselabels, trntag{ii});
    end
    
    lang_trn_set = block_trn_set & lang_set; 
    noi_trn_set  = block_trn_set & noi_set; 
    
    repfactor = sum(lang_trn_set)/sum(noi_trn_set);

%%% PARFOR START %%%
    parfor voxel_num = 1:numvox
        sphere_inds_for_this_voxel = sphere_XYZ_indices_cell{voxel_num};
        these_betas = betas_zero_mean(:,sphere_inds_for_this_voxel); %#ok<NODEF>

        %% Load data into machine
        % Data
        lang_trn_betas = these_betas(lang_trn_set, :); 
        noi_trn_betas  = these_betas(noi_trn_set, :); 
        noi_trn_betas  = repmat(noi_trn_betas, [repfactor 1]); 
        Instances_trn  = vertcat(lang_trn_betas, noi_trn_betas); 
        
        lang_tst_betas = these_betas(lang_tst_set, :); 
        noi_tst_betas  = these_betas(noi_tst_set, :); 
        noi_tst_betas  = repmat(noi_tst_betas, [repfactor 1]); 
        Instances_tst  = vertcat(lang_tst_betas, noi_tst_betas); 
        
        % Labels
        Labels_trn_lang =    ones(size(lang_trn_betas, 1), 1);
        Labels_trn_noi  = -1*ones(size(noi_trn_betas, 1),  1);
        Labels_tst_lang =    ones(size(lang_tst_betas, 1), 1);
        Labels_tst_noi  = -1*ones(size(noi_tst_betas, 1),  1);

        Labels_trn = [Labels_trn_lang; Labels_trn_noi];
        Labels_tst = [Labels_tst_lang; Labels_tst_noi];

        %%%%%%%%LibSVM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        model = svmtrain(Labels_trn, Instances_trn, ['-q', 'libsvm_options']);
        Outputs_testing = svmpredict(Labels_tst, Instances_tst, model, ...
            ['-q','libsvm_options']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%GNB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nb = fitcnb(Instances_trn, Labels_trn);
% %         nb = NaiveBayes.fit(Instances_trn,Labels_trn);
%         Outputs_testing = predict(nb, Instances_tst);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Accuracy
        num_testing_points = size(Instances_tst,1);
        testing_Prop_correct = sum(Outputs_testing == Labels_tst) / num_testing_points;
        crossval_acc_vector(voxel_num) = testing_Prop_correct-0.5;

        if rem(voxel_num,10000)==0
            disp([  'Test-set accuracy ' num2str(voxel_num) ...
                ' out of ' num2str(num_mask_voxels) ...
                ' = ' num2str(testing_Prop_correct)  ]);
        end

    end  %%% End of parfor loop through voxel spheres

    %% Merge vectors from each fold into matrices
    crossval_acc_mat_folds(fold, :) = crossval_acc_vector;

end  %%% End of loop through folds

%% Take the mean across folds
acc_vector = mean(crossval_acc_mat_folds, 1); 

%% Convert to brain-space matrix
acc = nan(Vmask.dim);
acc(mask_idx) = acc_vector; 

%% Save the result
file = fullfile(dir_design, 'beta_0001.nii'); 
V = spm_vol(file);

V.fname = fullfile(dir_data_MVPA, ... 
    [subj.name '_' design.name '_beta_SVM_lng_clear_noi_rad' num2str(radius) '.nii']);  
spm_write_vol(V, acc);

end
