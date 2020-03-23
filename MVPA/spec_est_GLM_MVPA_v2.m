%% spec_est_GLM
% Specifies and estimates first level GLM. 

% CHANGELOG (MM/DD/YY)
% 11/13/17  Changelog initialized (hybrid_isss) -- MH
% 02/12/20  Forked for YA_FST
% 02/17/20  Now specifies model with 1 run for visual inspection. Updated
%   image selection to use study.prefix
% 02/20/20  New way to drop physio
% 02/24/20  Silent trial is now implicitly modeled. 
% 03/02/20  Dropped first 5 images and first TR of each. New matlabbatch
%   which has autoregression turned OFF and updated FIR. There was an error
%   loading regressors--we selected the first run's regressor for every
%   run...
% 03/05/20  Added "incorrect" regressor when we model all events. Also
%   added flag for "do first run". 
% 03/05/20  Cloned for MVPA. 
% 03/10/20  Added interactive, first run flags. Now skips incorrect
%   regressor for perfect runs. 

function spec_est_GLM_MVPA_v2(varargin)
%% Check input
subj = varargin{1}; 
study = varargin{2}; 
design = varargin{3}; 
if length(varargin) > 3 % do first run
    do_first_run = varargin{4}; 
else
    do_first_run = 1; 
end

if length(varargin) > 4 % interactive
    do_interactive = varargin{5};
else
    do_interactive = 0; 
end

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

%% Pathing
cd ..
dir_batch = pwd; 
dir_mlb   = fullfile(dir_batch, 'matlabbatch'); 
dir_subj  = fullfile(study.path, 'data', subj.name); 
dir_reg   = fullfile(dir_subj, 'reg');
dir_func    = fullfile(dir_subj, 'FUNCTIONAL'); 
dir_ps_design   = fullfile(dir_subj, 'ps', 'designs');
dir_design_root = fullfile(dir_subj, 'design');

for dd = design
    thisdesign = study.design(dd); 
    if thisdesign.correct
        data = fullfile(dir_subj, 'onsets_nowrong_v3.mat'); 
    else
        data = fullfile(dir_subj, 'onsets_all.mat'); 
    end
    load(data)
        
    %% GLM with 1 run
    if do_first_run
        batch = fullfile(dir_mlb, 'SPM_specify_GLM_v2.mat'); 
        load(batch)
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

            try 
                theseevents = onsets.(thisdesign.cond{cc})(:, 1); 
            catch
                disp('Building onset vector...')
                eventtypes = fields(onsets); 
                tag = strsplit(thisdesign.cond{cc}, '_'); 

                idx = contains(eventtypes, tag{1}); 
                for tt = 2:length(tag)
                    idx = idx & contains(eventtypes, tag{tt}); 
                end

                events = eventtypes(idx); 
                theseevents = []; 
                for tt = 1:length(events)
                    theseevents = [theseevents; onsets.(events{tt})(:, 1)];
                end

            end

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

        % Remove second run and clean up info
        matlabbatch{1}.spm.stats.fmri_spec.sess(2) = []; 
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''}; 
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; 

        % Add regressors
        clear multiRegs
        multiRegs(1) = string([dir_reg filesep 'rp_' study.runname '1_00007.txt']);
        idx = 2; 

        if thisdesign.physio
            multiRegs(idx) = string([dir_reg filesep 'physio_1st_' study.runname '1_00001.txt']);
            idx = idx + 1; 
        end

        if ~thisdesign.correct 
            multiRegs(idx) = string(fullfile(dir_reg, 'incorrect_run1.txt'));
            idx = idx + 1; 
        end

        multiRegs = cellstr(multiRegs');
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = multiRegs;

        % Run the GLM-specify job
        filename = fullfile(dir_subj, 'batch', ['specify_GLM_' thisdesign.name '_1run.mat']); 
        save(filename, 'matlabbatch')

        mkdir(dir_thisdesign)
        spm_jobman('run',matlabbatch);
        if do_interactive
            disp('Waiting for OK to continue!'); pause
        end
        
        fg = spm_figure('FindWin','Graphics');
        figname = ['design_' thisdesign.name '_1run.png']; 
        saveas(fg, figname)
        movefile(fullfile(pwd, figname), fullfile(dir_ps_design));
    end
    
    %% GLM with all runs
    batch = fullfile(dir_mlb, 'SPM_specify_GLM_v2.mat'); 
    load(batch)
    thisdesign = study.design(dd); 
    dir_thisdesign = fullfile(dir_design_root, thisdesign.name); 
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {dir_thisdesign}; 
    
    runs = 1:subj.runs; runs(subj.drop) = []; 
    for rr = 1:length(runs)
        %% Preallocate structure
        matlabbatch{1}.spm.stats.fmri_spec.sess(rr) = ...
        matlabbatch{1}.spm.stats.fmri_spec.sess(1); 
        
        %% Extract names for all bold files and put into job struct
        [boldFiles, ~] = spm_select('List', dir_func, ... 
            ['^' study.mvpa.prefix study.runname num2str(runs(rr)) '_00*.*\.nii$']);
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
            
            try 
                theseevents = onsets.(thisdesign.cond{cc})(:, runs(rr)); 
            catch
                disp('Building onset vector...')
                eventtypes = fields(onsets); 
                tag = strsplit(thisdesign.cond{cc}, '_'); 

                idx = contains(eventtypes, tag{1}); 
                for tt = 2:length(tag)
                    idx = idx & contains(eventtypes, tag{tt}); 
                end

                events = eventtypes(idx); 
                theseevents = []; 
                for tt = 1:length(events)
                    theseevents = [theseevents; onsets.(events{tt})(:, runs(rr))];
                end

            end
            
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
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(rr).multi = {''}; 
        
        %% Add regressors
        clear multiRegs
        multiRegs(1) = string([dir_reg filesep 'rp_' study.runname num2str(runs(rr)) '_00007.txt']);
        idx = 2; 
        if thisdesign.physio
            multiRegs(idx) = string([dir_reg filesep 'physio_1st_' study.runname num2str(runs(rr)) '_00001.txt']);
            idx = idx + 1; 
        end
        
        if ~thisdesign.correct 
            if ~perfect(runs(rr)) % from extract_timing_all!
                multiRegs(idx) = string(fullfile(dir_reg, ['incorrect_run' num2str(runs(rr)) '.txt']));
                idx = idx + 1; 
            end
            
        end
        
        multiRegs = cellstr(multiRegs');

        matlabbatch{1}.spm.stats.fmri_spec.sess(rr).multi_reg = multiRegs;
        
        %% Add highpass filter
        matlabbatch{1}.spm.stats.fmri_spec.sess(rr).hpf = 256; 
        
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
    if do_interactive
        disp('Waiting for OK to continue!'); pause
    end
        
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