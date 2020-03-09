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
% 02/24/20  Dropping first TR of each event due to artifact. Also not
%   modeling silent event. 

function spec_est_GLM_MVPA(subj, study, design)
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
dir_func  = fullfile(dir_subj, 'FUNCTIONAL_v2'); % grab unsmoothed data
dir_ps_design   = fullfile(dir_subj, 'ps', 'designs');
dir_design_root = fullfile(dir_subj, 'design');

%% Specify GLM
data = fullfile(dir_subj, 'onsets_MVPA_v2.mat'); 
load(data)

for dd = design
    %% Evaluate design for 1 run
    batch = fullfile(dir_mlb, 'SPM_specify_GLM.mat'); 
    load(batch)
    thisdesign = study.design(dd); 
    dir_thisdesign = fullfile(dir_design_root, [thisdesign.name '_1run']); 
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {dir_thisdesign}; 
    
    % Extract names for all bold files and put into job struct
    [boldFiles, ~] = spm_select('List', dir_func, ...
        ['^' study.mvpa.prefix study.runname '1_00*.*\.nii$']);
    boldFiles = [repmat([dir_func filesep], length(boldFiles), 1), boldFiles]; %#ok<AGROW>
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
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2) = []; 
    
    % Merge noise and silence
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).name = ...
        'NOI_SIL'; 

    theseevents = sort([onsets.(thisdesign.cond{1})(:, 1); onsets.(thisdesign.cond{2})(:, 1)]);
    theseevents = theseevents(~isnan(theseevents)); 
    if isempty(theseevents)
        flags(1) = true; 
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).onset = ...
        theseevents; 

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).duration = ...
        zeros(1, length(theseevents)); 
    
    % Clean up session info
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(2) = []; % drop other runs
    
    % Add regressors
    clear multiRegs
    multiRegs(1) = string([dir_reg filesep 'rp_' study.runname '1_00006.txt']);
    if thisdesign.physio
        multiRegs(2) = string([dir_reg filesep 'physio_1st_' study.runname '1_00006.txt']);
    end

    multiRegs = cellstr(multiRegs');

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = multiRegs;
    
    % Run the job
    filename = fullfile(dir_subj, 'batch', ['specify_GLM_' thisdesign.name '_1run.mat']); 
    save(filename, 'matlabbatch')
    
    mkdir(dir_thisdesign)
    spm_jobman('run',matlabbatch);
    fg = spm_figure('FindWin','Graphics');
    figname = ['design_' thisdesign.name '_1run.png']; 
    saveas(fg, figname)
    movefile(fullfile(pwd, figname), fullfile(dir_ps_design));
    
    %% Evaluate all runs
    load(batch)
    thisdesign = study.design(dd); 
    dir_thisdesign = fullfile(dir_design_root, thisdesign.name); 
    matlabbatch{1}.spm.stats.fmri_spec.dir = {dir_thisdesign}; 
    flags_run = false(1, subj.runs); 
    
    for rr = 1:subj.runs
        %% Preallocate structure
        matlabbatch{1}.spm.stats.fmri_spec.sess(rr) = ...
        matlabbatch{1}.spm.stats.fmri_spec.sess(1); 
        
        %% Extract names for all bold files and put into job struct
        [boldFiles, ~] = spm_select('List', dir_func, ['^' study.mvpa.prefix study.runname num2str(rr) '_00*.*\.nii$']);
        boldFiles = [repmat([dir_func filesep], length(boldFiles), 1), boldFiles]; %#ok<AGROW>
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
%         matlabbatch{1}.spm.stats.fmri_spec.sess(rr).regress = regress; 
        
        %% Add regressors
        clear multiRegs
        multiRegs(1) = string([dir_reg filesep 'rp_' study.runname num2str(rr) '_00006.txt']);
        if thisdesign.physio
            multiRegs(2) = string([dir_reg filesep 'physio_1st_' study.runname num2str(rr) '_00006.txt']);
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
    
    %% Run the GLM-estimate job
    batch = fullfile(dir_mlb, 'SPM_estimate_GLM.mat'); 
    load(batch)
    
    matlabbatch{1}.spm.stats(1).fmri_est.spmmat = {[dir_thisdesign filesep 'SPM.mat']};
    
    filename = fullfile(dir_subj, 'batch', ['estimate_GLM_' thisdesign.name '.mat']); 
    save(filename, 'matlabbatch')
    spm_jobman('run',matlabbatch);
end

end