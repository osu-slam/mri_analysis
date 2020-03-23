%% YA_FST_unpack
% Unpacks data for all subjects. 

% MM/DD/YY -- CHANGELOG
% 11/08/17 -- File initialized
% 02/05/20 -- Forked for YA_FST
% 02/10/20 -- V1 done and run for first 3 subjects. MJH
% 02/21/20 -- New design with shifted onsets. 
% 03/09/20 -- Updated (most) code to follow subj/study convention!

clc; clearvars

%% Pathing
dir_setup = fullfile(pwd, 'setup');
YA_FST_params

for ii = 1:length(subj)
    %% Subject-specific parameters
    thissubj = subj(ii);
    disp(['Unpacking ' thissubj.name '...'])
    
    %% Make directories
%     disp('Making directories')
%     cd(dir_setup); make_dir(thissubj, study)
%     disp('Done!')
    
    %% Unzip files
%     disp('Unzipping all files now...')
%     cd(dir_setup); unzip_data(thissubj, study)
%     disp('Done!')

     %% Convert dicms
%     disp('Converting all dicms...')
%     cd(dir_setup); convert_dicm(thissubj, study)
%     cd(dir_setup); convert_dicm_physio(thissubj, study)
%     disp('Done!')
    
    %% Rename files
%     if ~isempty(subj(ii).rename)
%         disp('Renaming necessary runs...')
%         cd(dir_setup); rename_files % functional data
%         disp('Done!')
%     end
    
%     if ~isempty(subj(ii).appa)
%         disp('Renaming AP-PA images...')
%         cd(dir_setup); rename_appa(thissubj, study)
%         disp('Done!')
%     end
    
    %% Make physio regressors
    disp('Making physio regressors...')
    cd(dir_setup); make_phys_reg(thissubj, study)
    disp('Done!')
    
    %% Make fieldmap
%     % DONE EXTERNALLY WITH unwarp.sh
%     path = fullfile(dir_batch, 'unwarp.sh'); 
%     cmd  = [path ' ' thissubj.name]; 
%     disp('Creating fieldmap with FSL...')
%     cd(dir_batch); system(cmd); 
%     disp('Done!')
% 
%     disp(['Finished unpacking ' thissubj.name '.'])
end
