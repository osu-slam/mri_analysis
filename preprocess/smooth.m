%% smooth
% Adds Gaussian smoothing to functional data for GLM processing. 
% 
% MM/DD/YY: CHANGELOG
% 02/11/20: Forked from isss_multi, made universal

function smooth(subj, study)
%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
disp('Done!')

%% Pathing
cd ..
dir_batch = pwd; 
batch = fullfile(dir_batch, 'matlabbatch', 'SPM_smooth_GLM.mat'); 
load(batch)

cd ..
cd ..
dir_spm = fullfile(pwd, 'spm12'); 
cd(dir_batch)

dir_subj = fullfile(study.path, 'data', subj.name); 

dir_anat = fullfile(dir_subj, 'ANATOMICAL'); 
dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 
dir_thissubj_batch = fullfile(dir_subj, 'batch'); %% Pathing
dir_glm  = fullfile(dir_subj, 'FUNC_GLM');

%% Poke values into batch
func_data = {}; 
for rr = 1:subj.runs
    [boldFiles, ~] = spm_select('List', dir_func, ['^wu' study.runname num2str(rr) '_00*.*\.nii$']);
    boldFiles = [repmat([dir_func filesep], length(boldFiles),1), boldFiles];
    boldFiles = cellstr(boldFiles);
    func_data = vertcat(func_data, boldFiles); 
end
matlabbatch{1}.spm.spatial.smooth.data = func_data;

filename = fullfile(dir_thissubj_batch, 'smooth_GLM.mat'); 
save(filename, 'matlabbatch')

spm_jobman('run',matlabbatch);
disp('Done with GLM smoothing!')

%% Move files
disp('Moving files to GLM directory')
target = fullfile(dir_func, ['swu' study.runname '*']);
movefile(target, dir_glm);

end