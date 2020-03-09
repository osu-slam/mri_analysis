%% build_contrasts
% Builds contrasts for 1st level analysis for one subject

% CHANGELOG 
% 11/15/17  Initialized file (hybrid_isss)
% 02/13/20  Forked for YA_FST, made universal
% 03/02/20  Dropping first TR from each event

function build_contrasts(subj, study, design)
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
disp('Done!')

%% Pathing
cd ..
dir_batch = pwd; 
dir_mlb = fullfile(dir_batch, 'matlabbatch'); 
dir_subj = fullfile(study.path, 'data', subj.name); 
dir_design_root = fullfile(dir_subj, 'design');

for dd = design
    thisdesign = study.design(dd); 
    
    %% Load matlabbatch, add SPM.mat
    batch = fullfile(dir_mlb, 'SPM_build_contrasts.mat'); 
    load(batch)
    
    spmmat = fullfile(dir_design_root, thisdesign.name, 'SPM.mat');
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(spmmat);
    
    %% For each contrast...
    for cc = 1:length(thisdesign.con.name)
        if contains(thisdesign.name, 'HRF')
            thisvec = thisdesign.con.vec(cc, :); 
        else
            thisvec = repelem(thisdesign.con.vec(cc, :), study.scan.epis-1); 
        end
        
        matlabbatch{1}.spm.stats.con.consess{cc} = ... 
            matlabbatch{1}.spm.stats.con.consess{1}; % Preallocate
        
        matlabbatch{1}.spm.stats.con.consess{cc}.tcon.name = ...
            thisdesign.con.name{cc};
        matlabbatch{1}.spm.stats.con.consess{cc}.tcon.weights = ...
            thisvec;     
        matlabbatch{1}.spm.stats.con.consess{cc}.tcon.sessrep = ...
            'repl'; 
    end
    
    if size(matlabbatch{1}.spm.stats.con.consess, 2) > cc
        matlabbatch{1}.spm.stats.con.consess(cc+1:end) = [];
    end
    
    %% Save files
    filename = fullfile(dir_subj, 'batch', ['build_contrasts_' thisdesign.name '.mat']);
    save(filename, 'matlabbatch');
    
    % Run job
    spm_jobman('run', matlabbatch);
end
