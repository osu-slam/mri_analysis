%% spec_est_GLM
% Specifies and estimates first level GLM. 

% CHANGELOG (MM/DD/YY)
% 11/13/17  Changelog initialized -- MH
% 11/13/17  Updated to match new file conventions, added exception for
%   physio-less subjects, removed ART from code, added code to save 
%   functions -- MH
% 02/28/18  Running new sub-analysis to look at which acquisition time
%   window gives best SNR. Splitting into new file. -- MH
% 05/30/18  Updated for 14 subjects. Also updated to use filesep instead of
%   back or forward slash, should work on all OS. -- MH
% 02/12/20  Forked for YA_FST
% 02/17/20  Now specifies model with 1 run for visual inspection. Updated
%   image selection to use study.prefix
% 02/20/20  New way to drop physio
% 02/24/20  Silent trial is now implicitly modeled. 
% 02/26/20  Tried to make an HRF version. This function does not work. 

function spec_est_GLM_HRF(subj, study, design)
%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

if ~isnumeric(design)
    error('Input ("dd") which specifies which designs (#) should be used')
end

%% Prepare SPM
disp('Loading SPM...')
close all; 
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_figure('GetWin','Graphics'); % Thanks Guillaume
disp('Done!')

%% Parameters (CHECK THESE BEFORE RUNNING)
% pmod = struct('name', [], 'param', [], 'poly', []); 
% pmod(1) = []; 
% regress = struct('name', [], 'val', []); 
% regress(1) = []; 

%% Pathing
cd ..
dir_batch = pwd; 
dir_mlb   = fullfile(dir_batch, 'matlabbatch'); 
dir_subj  = fullfile(study.path, 'data', subj.name); 
dir_reg   = fullfile(dir_subj, 'reg');
dir_func_GLM    = fullfile(dir_subj, 'FUNC_GLM'); 
dir_ps_design   = fullfile(dir_subj, 'ps', 'designs');
dir_design_root = fullfile(dir_subj, 'design');

%% Specify GLM
data = fullfile(dir_subj, 'onsets_nowrong_v2.mat'); % just for now i guess
load(data)

for dd = design
    %% GLM with 1 run
    batch = fullfile(dir_mlb, 'SPM_specify_GLM_HRF.mat'); 
    load(batch)
    thisdesign = study.design(dd); 
    dir_thisdesign = fullfile(dir_design_root, [thisdesign.name '_1run']); 
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {dir_thisdesign}; 
    
    % Extract names for all bold files and put into job struct
    [boldFiles, ~] = spm_select('List', dir_func_GLM, ...
        ['^' study.prefix study.runname '1_00*.*\.nii$']);
    
    boldFiles = [repmat([dir_func_GLM filesep], length(boldFiles), 1), boldFiles]; %#ok<AGROW> 
    boldFiles = cellstr(boldFiles); matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = boldFiles;

    % Condition names and onsets. 
    numCons = length(thisdesign.cond);
    flags = false(1, numCons);

    for cc=1:numCons
        % Preallocate structure
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cc).pmod = ...
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod; 

        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cc).name = ...
            thisdesign.cond{cc}; 

        theseevents = onsets.(thisdesign.cond{cc})(:, 1); 
        theseevents = theseevents(~isnan(theseevents)); 
        if isempty(theseevents)
            flags(cc) = true; 
        end

        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cc).onset = ...
            theseevents; 

        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cc).duration = ...
            zeros(1, length(theseevents)); 

        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cc).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cc).orth = 1; 
        
        flags_run = any(flags); 
    end
    
    % Drop silence to prevent non-unique betas
%     matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2) = []; 
    
    % Remove second run and clean up info
    matlabbatch{1}.spm.stats.fmri_spec.sess(2) = []; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''}; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; 
    
    % Add regressors
    clear multiRegs
    multiRegs(1) = string([dir_reg filesep 'rp_' study.runname '1_00006.txt']);
    if thisdesign.physio
        multiRegs(2) = string([dir_reg filesep 'physio_1st_' study.runname '1_00006.txt']);
    end

    multiRegs = cellstr(multiRegs');
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = multiRegs;
    
    % Run the GLM-specify job
    filename = fullfile(dir_subj, 'batch', ['specify_GLM_' thisdesign.name '_1run_HRF.mat']); 
    save(filename, 'matlabbatch')
    
    mkdir(dir_thisdesign)
    spm_jobman('run',matlabbatch);
    
    %% Modify the regressors before estimating
    file_spm = fullfile(dir_thisdesign, 'SPM.mat');
    copy_spm = fullfile(dir_thisdesign, 'SPM_original.mat');
    copyfile(file_spm, copy_spm); 
    load(file_spm)
    
    % Lift and modify the HRF from other trials
    fmri_t = matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t; 
    template = logical([ones(study.scan.epis, 1); zeros(study.scan.silence, 1)]); 
    idx = 1; hrf_mod = false(fmri_t*3, 1); 
    for xx = 1:fmri_t*3
        hrf_mod(xx) = template(idx); 
        if idx == length(template)
            idx = 0; 
        end
        
        idx = idx + 1; 
    end
    
    matrix = SPM.xX.X; 
    for cc = 1:numCons
        key = find(SPM.xX.X(:, cc)); 
        diffkey = diff(key); diffkey = [1; diffkey];
        if any(diffkey > 1)
            gap = find(diffkey > 1); 
            
            key_parts = cell(length(gap) + 1, 1); 
            key_parts{1} = key(1:gap(1)-1); 
            key_parts{end} = key(gap(end):end); 
            if length(gap) > 1
                for xx = 2:length(gap)-1
                    key_parts{xx} = key(gap(xx):gap(xx+1)-1); 
                end
                
            end
            
        else
            key_parts = cell(1, 1); 
            key_parts{1} = key;
        end
        
        for xx = 1:length(key_parts)
            this_key = key_parts{xx};
            hrf = SPM.xX.X(this_key, cc);  
            new_hrf = hrf(hrf_mod(1:length(this_key))); 
            new_key = this_key(hrf_mod(1:length(this_key))); 

            first = this_key(1); last = this_key(1)+length(new_key)-1;  
            matrix(:, cc) = zeros(size(matrix(:, cc))); % clear old value;
            matrix(first:last, cc) = new_hrf; 
        end
        
    end
    
    SPM.xX.X = matrix; 
    
    % Save the new design
    save(file_spm, 'SPM')
    
    %% GLM with all runs
    batch = fullfile(dir_mlb, 'SPM_specify_GLM_HRF.mat'); 
    load(batch)
    thisdesign = study.design(dd); 
    dir_thisdesign = fullfile(dir_design_root, thisdesign.name); 
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {dir_thisdesign}; 
    
    for rr = 1:subj.runs
        %% Preallocate structure
        matlabbatch{1}.spm.stats.fmri_spec.sess(rr) = ...
        matlabbatch{1}.spm.stats.fmri_spec.sess(1); 
        
        %% Extract names for all bold files and put into job struct
        [boldFiles, ~] = spm_select('List', dir_func_GLM, ... 
            ['^' study.prefix study.runname num2str(rr) '_00*.*\.nii$']);
        boldFiles = [repmat([dir_func_GLM filesep], length(boldFiles), 1), boldFiles]; %#ok<AGROW>
        boldFiles = cellstr(boldFiles);
        matlabbatch{1}.spm.stats.fmri_spec.sess(rr).scans = boldFiles;

        %% Condition names and onsets. 
        numCons = length(thisdesign.cond);
        flags = false(1, numCons);
        
        for cc=1:numCons
            % Preallocate structure
            matlabbatch{1}.spm.stats.fmri_spec.sess(rr).cond(cc).pmod = ...
                matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod; 
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(rr).cond(cc).name = ...
                thisdesign.cond{cc}; 
            
            theseevents = onsets.(thisdesign.cond{cc})(:, rr); 
            theseevents = theseevents(~isnan(theseevents)); 
            if isempty(theseevents)
                flags(cc) = true; 
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(rr).cond(cc).onset = ...
                theseevents; 
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(rr).cond(cc).duration = ...
                zeros(1, length(theseevents)); 
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(rr).cond(cc).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(rr).cond(cc).orth = 1; 
        end
        
        %% Finish cleaning session info
        matlabbatch{1}.spm.stats.fmri_spec.sess(rr).multi = {''}; 
        
        %% Add regressors
        clear multiRegs
        multiRegs(1) = string([dir_reg filesep 'rp_' study.runname '1_00006.txt']);
        if thisdesign.physio
            multiRegs(2) = string([dir_reg filesep 'physio_1st_' study.runname '1_00001.txt']);
        end
        
        multiRegs = cellstr(multiRegs');

        matlabbatch{1}.spm.stats.fmri_spec.sess(rr).multi_reg = multiRegs;
        
        %% Add highpass filter
        matlabbatch{1}.spm.stats.fmri_spec.sess(rr).hpf = 128; 
        
        flags_run(rr) = any(flags); 
    end
    
    %% Drop any runs that have a gap. 
    if any(flags_run)
        warning(['Run ' find(flags_run) ' has an empty regressor! Dropping run...'])
        matlabbatch{1}.spm.stats.fmri_spec.sess(flags_run) = []; 
    end
    
    %% Run the GLM-specify job
    filename = fullfile(dir_subj, 'batch', ['specify_GLM_' thisdesign.name '.mat']); 
    save(filename, 'matlabbatch')
    
    mkdir(dir_thisdesign)
    spm_jobman('run',matlabbatch);
    fg = spm_figure('FindWin','Graphics');
    figname = ['design_' thisdesign.name '.png']; 
    saveas(fg, figname)
    movefile(fullfile(pwd, figname), fullfile(dir_ps_design));
    
    %% Modify the regressors before estimating
    file_spm = fullfile(dir_thisdesign, 'SPM.mat');
    copy_spm = fullfile(dir_thisdesign, 'SPM_original.mat');
    copyfile(file_spm, copy_spm); 
    load(file_spm)
    
    fmri_t = matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t; 
    names  = SPM.xX.name; 
    
    conds = false(size(names)); 
    for cc = 1:numCons
        conds = conds | contains(names, thisdesign.cond{cc});
    end
    
    % Build HRF template
    template = logical([ones(study.scan.epis, 1); zeros(study.scan.silence, 1)]); 
    idx = 1; hrf_mod = false(fmri_t*3); 
    for xx = 1:fmri_t*3 % max 3 consecutive events I guess...
        hrf_mod(xx) = template(idx); 
        if idx == length(template)
            idx = 0; 
        end
        
        idx = idx + 1; 
    end
    
    % Make new HRF
    matrix_conds = SPM.xX.X(:, conds);
    for cc = 1:size(matrix_conds, 2)
        key = find(matrix_conds(:, cc)); 
        hrf = matrix(key, cc); 
        
        new_hrf = hrf(hrf_mod(1:length(key))); 
        new_key = key(hrf_mod(1:length(key))); 
        
        first = key(1); last = key(1)+length(new_key)-1;  
        matrix_conds(:, cc) = zeros(size(matrix_conds(:, cc))); % clear old value;
        matrix_conds(first:last, cc) = new_hrf; 
    end
    
    SPM.xX.X(:, conds) = matrix_conds; 
    
    % Poke back into SPM.xX.X (NEEDS WORK!)
    matrix_new = matrix_conds; 
    for cc = 1:size(matrix_conds, 2)
        key = find(matrix_conds(:, cc)); 
        
        new_hrf = hrf(hrf_mod(1:length(key))); 
        new_t = length(new_hrf);
        
        first = key(1); 
        matrix_new(:, cc) = zeros(size(matrix_conds(:, cc))); % clear old value;
        matrix_new(first:first+new_t-1, cc) = new_hrf; 
    end
    
    SPM.xX.X(:, conds) = matrix_new; 
    
    % Save the new design
    save(file_spm, 'SPM')
    
    %% Run the GLM-estimate job
    batch = fullfile(dir_mlb, 'SPM_estimate_GLM.mat'); 
    load(batch)
    
    matlabbatch{1}.spm.stats(1).fmri_est.spmmat = {[dir_thisdesign filesep 'SPM.mat']};
    
    filename = fullfile(dir_subj, 'batch', ['estimate_GLM_' thisdesign.name '.mat']); 
    save(filename, 'matlabbatch')
    spm_jobman('run',matlabbatch);
end

end