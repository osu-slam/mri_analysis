%% make_dir
% Creates directories for a subject
% MM/DD/YY -- CHANGELOG
% 02/05/20 -- Changelog started, updated code to run for any experiment.
% 03/09/20 -- Updated to use subj/study input. 

function make_dir(subj, study)
%% Check input
if ~isstruct(subj)
    err('Input (subj) where subj is a SINGLE structure')
end

if ~isstruct(study)
    err('Input (study), which is a structure of study info!')
end

%% Path
thissubj = subj.name; 
dir_subj = fullfile(study.path, 'data', thissubj); 

%% Do the things
try
    cd(dir_subj)
    disp(['Directories already exist for ' thissubj])
    cd ..
catch
    disp(['Creating directory for ' thissubj])
    mkdir(dir_subj)
    cd(dir_subj)

    mkdir('ANATOMICAL')
%     mkdir('art')
    mkdir('batch')

    mkdir('behav')
    cd('behav')
    mkdir('analysis')
    mkdir('prescan')
    mkdir('scan')
    cd ..

    mkdir('design')
    mkdir('dicm')
    mkdir('FUNC_GLM')
    mkdir('FUNC_MVPA')
    mkdir('FUNCTIONAL')
    mkdir('nii')

    mkdir('PHYSIO')

    mkdir('ps')
    cd('ps')
    mkdir('contrast')
    mkdir('designs')
    mkdir('preprocessing')
    cd ..

    mkdir('realign')
    mkdir('reg')
    mkdir('SNR')
    mkdir('zip')

end

end