%% searchlight_or_sr
% Does a binary classification using searchlight-based SVM (I think).
% Adapted from Yune's code on 12/14/2017
% Updated for SLAM_isss_multi_MVPA_batch on 6/13/18
% Made parallel on 6/22/18
% Split for noise/silence and lng/noi on 6/26
% 
% MM/DD/YY -- CHANGELOG
% 02/14/20 -- Cloned for YA_FST. Made into function. 
% 02/14/20 -- Cloned for YA_FST babble, a 3-way classification

function searchlight_babble(subj, study, dd)
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
filename = fullfile(dir_data_MVPA, [subj.name '_zero_meaned_betas.mat']);
load(filename)

filename = fullfile(dir_data_MVPA, [subj.name '_spheres_radius' num2str(radius) '.mat']);  
load(filename);  

%% More parameters and preallocation
numvox = length(sphere_XYZ_indices_cell);  %#ok<USENS>

crossval_ora_HIT_mat_folds = zeros(subj.runs, numvox);
crossval_sra_HIT_mat_folds = zeros(subj.runs, numvox);
crossval_ora_FA_mat_folds  = zeros(subj.runs, numvox);
crossval_sra_FA_mat_folds  = zeros(subj.runs, numvox);
crossval_acc_mat_folds     = zeros(subj.runs, numvox);

blocks = 1:subj.runs; 

%% Read design matrix
%%%% Read the column-names from the design matrix, to figure out which
%%%% regressors to use for defining our time-points
% all_column_cond_names = char(SPM.xX.name);
% num_cols = size(SPM.xX.X,2);
theselabels = SPM.xX.name;
regressors = ~cellfun(@isempty, regexp(theselabels, 'R\d'));
constant = contains(theselabels, 'constant');
removethese = regressors | constant;
theselabels(removethese) = [];
numbetas = length(theselabels); 

for fold = 1:subj.runs % for each fold... 
    disp(['Cross-validation fold ' num2str(fold)]);
    
    %% Preallocate 
    crossval_clear_HIT_vector = zeros(1, numvox);
    crossval_snr2_HIT_vector = zeros(1, numvox);
    crossval_clear_FA_vector  = zeros(1, numvox);
    crossval_snr2_FA_vector  = zeros(1, numvox);
    crossval_acc_vector     = zeros(1, numvox);
    
    %% making training & testing sets 
    test = fold; 
    tsttag = {['Sn(' num2str(test) ')']}; 
    tst_block = contains(theselabels, tsttag{1});
    
    clear_tst_set = tst_block & contains(theselabels, 'clear');
    snr2_tst_set = tst_block & contains(theselabels, 'snr2'); 
    snr_2_tst_set = tst_block & contains(theselabels, 'snr_2'); 
    
    tst_all = clear_tst_set | snr2_tst_set | snr_2_tst_set;
    
    train = blocks; train(fold) = []; 
    trntag = cell(length(train), 1); 
    trn_block = false(1, numbetas); 
    for ii = 1:length(train)
        trntag{ii} = ['Sn(' num2str(train(ii)) ')']; 
        trn_block = trn_block | contains(theselabels, trntag{ii});
    end
    
    clear_trn_set = false(1, numbetas); 
    snr2_trn_set  = false(1, numbetas); 
    snr_2_trn_set = false(1, numbetas); 
    
    for ii = 1:length(train)
        clear_trn_set = clear_trn_set | (trn_block & contains(theselabels, 'clear'));
        snr2_trn_set  = snr2_trn_set  | (trn_block & contains(theselabels, 'snr2'));
        snr_2_trn_set = snr_2_trn_set | (trn_block & contains(theselabels, 'snr_2'));
    end
    
    trn_all = clear_trn_set | snr2_trn_set | snr_2_trn_set; 

%%% PARFOR START %%%
    for voxel_num = 1:numvox
        sphere_inds_for_this_voxel = sphere_XYZ_indices_cell{voxel_num};
        these_betas = betas_zero_mean(:,sphere_inds_for_this_voxel); %#ok<NODEF>

        %% Load data into machine

        %%% Each data-point is one time-slice, i.e. one row of this
        Instances_trn = these_betas(trn_all, :);
        Instances_tst = these_betas(tst_all, :); 

        Labels_trn_clear =    ones(length(find(clear_trn_set)), 1);
        Labels_trn_snr2  =  2*ones(length(find(snr2_trn_set)), 1);
        Labels_trn_snr_2 =  3*ones(length(find(snr_2_trn_set)), 1);
        Labels_tst_clear = 1*ones(length(find(clear_tst_set)), 1);
        Labels_tst_snr2  = 2*ones(length(find(snr2_tst_set)), 1);
        Labels_tst_snr_2 = 3*ones(length(find(snr_2_tst_set)), 1);

        Labels_trn = [Labels_trn_clear; Labels_trn_snr2; Labels_trn_snr_2];
        Labels_tst = [Labels_tst_clear; Labels_tst_snr2; Labels_tst_snr_2];

%         %%%%%%%%LibSVM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         model = fitcsvm(Instances_trn,Labels_trn);
%         [Outputs_testing, ~] = predict(model, Instances_tst);
%         testing_Prop_correct=accuracy(1,1);
%         cross_validation_mats(center_x,center_y,center_z,fold) = ... 
%             (testing_Prop_correct -12.5)/100;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%GNB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nb = fitcnb(Instances_trn, Labels_trn);
%         nb = NaiveBayes.fit(Instances_trn,Labels_trn);
        Outputs_testing = predict(nb, Instances_tst);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Accuracy
        Predicting_clear = Outputs_testing(Labels_tst==1); % All data that should be clear
        Predicting_snr2  = Outputs_testing(Labels_tst==2); % All data that should be snr2
        Predicting_snr_2 = Outputs_testing(Labels_tst==3); % All data that should be snr-2

        % Accuracies of prediction
        Hit_clear = sum(Predicting_clear==Labels_tst(Labels_tst== 1))/length(Predicting_clear);
        Hit_snr2  = sum(Predicting_snr2==Labels_tst(Labels_tst==2))/length(Predicting_snr2);
        Hit_snr_2  = sum(Predicting_snr_2==Labels_tst(Labels_tst==3))/length(Predicting_snr_2);

        % False alarm rate
        Predicting_non_clear = Outputs_testing(~(Labels_tst== 1)); % All data that should NOT be clear
        Predicting_non_snr2  = Outputs_testing(~(Labels_tst== 2)); % All data that should NOT be snr2
        Predicting_non_snr_2 = Outputs_testing(~(Labels_tst== 3)); % All data that should NOT be snr-2

        % False alarm
        FA_clear = sum(Predicting_non_clear == 1*ones(length(Predicting_non_clear),1)) / length(Predicting_non_clear);
        FA_snr2  = sum(Predicting_non_snr2  == 2*ones(length(Predicting_non_snr2),1)) / length(Predicting_non_snr2);
        FA_snr_2 = sum(Predicting_non_snr_2 == 3*ones(length(Predicting_non_snr_2),1)) / length(Predicting_non_snr_2);

        num_testing_points   = size(Instances_tst,1);
        testing_Prop_correct = sum(Outputs_testing == Labels_tst) / num_testing_points;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% store then in the matrix

        % Results are logit compatible
        if Hit_clear ~=0
            crossval_clear_HIT_vector(voxel_num) = Hit_clear; 
        elseif Hit_clear == 0
            crossval_clear_HIT_vector(voxel_num) = 1/length(Predicting_clear);
        elseif Hit_clear == 1
            crossval_clear_HIT_vector(voxel_num) = (length(Predicting_clear)-1)/length(Predicting_clear);
        end

        if Hit_snr2 ~=0
            crossval_snr2_HIT_vector(voxel_num) = Hit_snr2; 
        elseif Hit_snr2 == 0
            crossval_snr2_HIT_vector(voxel_num) = 1/length(Predicting_snr2);
        elseif Hit_snr2 == 1
            crossval_snr2_HIT_vector(voxel_num) = (length(Predicting_snr2)-1)/length(Predicting_snr2);
        end

        if Hit_snr_2 ~=0
            crossval_snr_2_HIT_vector(voxel_num) = Hit_snr_2; 
        elseif Hit_snr_2 == 0
            crossval_snr_2_HIT_vector(voxel_num) = 1/length(Predicting_snr_2);
        elseif Hit_snr_2 == 1
            crossval_snr_2_HIT_vector(voxel_num) = (length(Predicting_snr_2)-1)/length(Predicting_snr_2);
        end
        
        if FA_clear ~=0
            crossval_clear_FA_vector(voxel_num) = FA_clear;
        elseif FA_clear ==0
            crossval_clear_FA_vector(voxel_num) = 1/length(Predicting_non_clear);
        elseif FA_clear ==1
            crossval_clear_FA_vector(voxel_num) = (length(Predicting_non_clear)-1)/length(Predicting_non_clear);
        end

        if FA_snr2 ~=0
            crossval_snr2_FA_vector(voxel_num) = FA_snr2;
        elseif FA_snr2 ==0
            crossval_snr2_FA_vector(voxel_num) = 1/length(Predicting_non_snr2);
        elseif FA_snr2 ==1
            crossval_snr2_FA_vector(voxel_num) = (length(Predicting_non_snr2)-1)/length(Predicting_non_snr2);
        end
        
        if FA_snr_2 ~=0
            crossval_snr_2_FA_vector(voxel_num) = FA_snr_2;
        elseif FA_snr_2 ==0
            crossval_snr_2_FA_vector(voxel_num) = 1/length(Predicting_non_snr_2);
        elseif FA_snr_2 ==1
            crossval_snr_2_FA_vector(voxel_num) = (length(Predicting_non_snr_2)-1)/length(Predicting_non_snr_2);
        end

        % Accuracy is not logit-compatible
        crossval_acc_vector(voxel_num) = testing_Prop_correct-(1/3); 
        
        if rem(voxel_num,10000)==0
            disp([  'Test-set accuracy ' num2str(voxel_num) ...
                ' out of ' num2str(num_mask_voxels) ...
                ' = ' num2str(testing_Prop_correct)  ]);
        end

    end  %%% End of parfor loop through voxel spheres

    %% Merge vectors from each fold into matrices: PICK UP HERE!
    crossval_ora_HIT_mat_folds(fold, :) = crossval_clear_HIT_vector;
    crossval_sra_HIT_mat_folds(fold, :) = crossval_snr2_HIT_vector;
    crossval_ora_FA_mat_folds(fold, :)  = crossval_clear_FA_vector;
    crossval_sra_FA_mat_folds(fold, :)  = crossval_snr2_FA_vector;
    crossval_acc_mat_folds(fold, :)     = crossval_acc_vector;

end  %%% End of loop through folds

%% Take the mean across folds
ora_HIT_vector = mean(crossval_ora_HIT_mat_folds, 1); 
sra_HIT_vector = mean(crossval_sra_HIT_mat_folds, 1); 
ora_FA_vector  = mean(crossval_ora_FA_mat_folds, 1); 
sra_FA_vector  = mean(crossval_sra_FA_mat_folds, 1); 
acc_vector     = mean(crossval_acc_mat_folds, 1); 

%% Convert to brain-space matrix
ora_HIT = zeros(Vmask.dim);
sra_HIT = zeros(Vmask.dim);
ora_FA  = zeros(Vmask.dim);
sra_FA  = zeros(Vmask.dim);
acc     = zeros(Vmask.dim);

ora_HIT(mask_idx) = ora_HIT_vector; 
sra_HIT(mask_idx) = sra_HIT_vector; 
ora_FA(mask_idx)  = ora_FA_vector; 
sra_FA(mask_idx)  = sra_FA_vector; 
acc(mask_idx)     = acc_vector; 

%% Save the result
file = fullfile(dir_design, 'beta_0001.nii'); 
V = spm_vol(file);

dir_tp = fullfile(dir_data_MVPA, 'BETAS');
mkdir(dir_tp)

V.fname = fullfile(dir_tp, [subj.name '_BETA_GNB_ora_against_sra_HIT_rad' num2str(radius) '.nii']);  
spm_write_vol(V, ora_HIT);

V.fname = fullfile(dir_tp, [subj.name '_BETA_GNB_sra_against_ora_HIT_rad' num2str(radius) '.nii']);  
spm_write_vol(V,mean(sra_HIT,4));

V.fname = fullfile(dir_tp, [subj.name '_BETA_GNB_ora_for_sra_FA_rad' num2str(radius) '.nii']);  
spm_write_vol(V,mean(ora_FA,4));

V.fname = fullfile(dir_tp, [subj.name '_BETA_GNB_sra_for_ora_FA_rad' num2str(radius) '.nii']);  
spm_write_vol(V,mean(sra_FA,4));

V.fname = fullfile(dir_tp, [subj.name '_BETA_GNB_ora_sra_rad' num2str(radius) '.nii']);  
spm_write_vol(V,mean(acc,4));

end