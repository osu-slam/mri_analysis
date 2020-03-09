%% extract_timing_MVPA
% Loads and processes timing information. 
% Input: subj

% CHANGELOG
% 09/01/17  File inception
% 09/06/17  Worked on timing and duration extraction
% 09/07/17  Added button press condition, accuracy regressor
% 09/08/17  Changed accuracy regressor to a condition, started function
% 09/15/17  Updated to use cells to store onsets, events, and durations
% 09/26/17  Changing to use cells instead of structs. Way easier to handle
%   in SPM batch
% 09/29/17  Combined with convert to scans. 
% 11/08/17  Updated with new naming conventions
% 11/21/17  Went through to try and catch errors. Turned into a new,
%   simpler script. 
% 02/07/20  Forked for YA_FST. Loads in xlsx files. 
% 02/10/20  v1 finished, timing extracted for 3 subjects. 
% 02/13/20  Errors found, trying again. 
% 02/14/20  New version for MVPA, does not remove incorrect events. 

function extract_timing_MVPA(subj, study)
%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input (subj, study) where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input (subj, study) where study is struct')
end

%% Params and preallocation
events.NOI = 1; 
events.SIL = 1; 
events.lang = 3; % Needs to be updated for each experiment. 
events.synt = 9; 
events.all = 20; 
% Can I eliminate this? Not sure, requires lots of effort right now...
 
scan = study.scan; 

onsets = struct(... % Needs to be updated for each experiment. 
    'NOI', [], ...
    'SIL', [],  ...
    'OR_075_clear', [], ...
    'OR_075_snr2', [], ...
    'OR_075_snr_2', [], ...
    'OR_1_clear', [], ...
    'OR_1_snr2', [], ...
    'OR_1_snr_2', [], ...
    'OR_125_clear', [], ...
    'OR_125_snr2', [], ...
    'OR_125_snr_2', [], ...
    'SR_075_clear', [], ...
    'SR_075_snr2', [], ...
    'SR_075_snr_2', [], ...
    'SR_1_clear', [], ...
    'SR_1_snr2', [], ...
    'SR_1_snr_2', [], ...
    'SR_125_clear', [], ...
    'SR_125_snr2', [], ...
    'SR_125_snr_2', [] ...
    ); 

conNames    = fields(onsets);
numCons     = length(conNames); 
scan.order  = scan.first:scan.epis:scan.epis*(events.all+1); 

%% Load data for subject
disp('Now loading behavioral data...')
dir_subj = fullfile(study.path, 'data', subj.name); 
dir_behav = fullfile(dir_subj, 'behav', 'scan'); 
cd(dir_behav)
var = dir('*.xlsx');
disp(['Found ' num2str(length(var)) ' variable files.']) % There might be many from aborted runs

if isempty(var)
    error('No files found!')
elseif length(var) ~= 1
    %% Combining aborted runs (NEEDS WORK!)
    disp('Let us begin the splicing!')
    
    for vv = 1:length(var)
        % Load all tables
        thisT = readtable(var(vv).name); 
        if vv == 1
            T = thisT; % preallocate table
        end
        
        % Check which blocks we need from each
        blocks = unique(thisT.BLOCK); 
        for bb = 1:length(blocks)
            t = thisT(thisT.BLOCK == blocks(bb), :); 
            key = isnan(t.ActualEventDuration); 
            
            if ~any(key(1:end-1)) % old version has 1 nan at end
                T(T.BLOCK == blocks(bb), :) = t; 
            end
            
        end
        
    end
    
elseif length(var) == 1
    T = readtable(var.name); 
else
    error('Unknown error when loading files')
end

disp('Done!')

%% Populate onsets structure with preallocated matrixes
% Needs to be updated per experiment. 
for cc = 1:numCons
    onsets.(conNames{cc}) = nan(1, subj.runs); 
end

%% Get onset info for each event (correct and incorrect)
blocks = 1:subj.runs; 

for bb = blocks
    t = T(T.BLOCK == bb, :); 
    
    babble = t.Babble; 
    rate = t.Rate; 
    syntax = t.Syntax; 
    sentence = t.Sentence; 
    
    idx = 0; 
    idx_cond = zeros(numCons, 1); 
    for ev = 1:height(t)
        idx = idx + 1; 
        
        if isstrprop(sentence{ev}, 'alpha') % noise or silence
            if strcmp(sentence{ev}, 'N') % noise
                onsets.NOI(bb) = scan.order(idx); 
            elseif strcmp(sentence{ev}, 'S') % silence
                onsets.SIL(bb) = scan.order(idx); 
            else
                error('unknown stim')
            end

        elseif isstrprop(sentence{ev}, 'digit') % sentence presentation
            % Construct what event it was
            thissent = [syntax{ev} 'R_']; % start with syntax 

            if strcmp(rate{ev}, '1')
                tag_r = '1'; 
            elseif strcmp(rate{ev}, '1.25')
                tag_r = '125'; 
            elseif strcmp(rate{ev}, '0.75')
                tag_r = '075'; 
            else
                error('unknown stim')
            end

            if isempty(babble{ev})
                tag_b = 'clear'; 
            elseif strcmp(babble{ev}, 'SNR-2')
                tag_b = 'snr_2'; 
            elseif strcmp(babble{ev}, 'SNR2')
                tag_b = 'snr2'; 
            else
                error('unknown stim')
            end

            thissent_all = [thissent tag_r '_' tag_b]; % for all types
            thiscond = strcmp(conNames, thissent_all); 
            idx_cond(thiscond) = idx_cond(thiscond) + 1; 
            onsets.(thissent_all)(idx_cond(thiscond), bb) = scan.order(idx); 

        else
            error('unknown stim')
        end
        
    end
    
end

%% Save
filename = fullfile(dir_subj, 'onsets_MVPA.mat'); 
save(filename, 'onsets')

end