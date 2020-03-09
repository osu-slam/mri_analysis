%% extract_timing_nowrong
% Loads and processes timing information. Spits out a regressor of
% incorrect scans. 
% Input: dir_subj

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
% 02/21/20  New timing scheme. 
% 02/29/20  Updated for new version of experiment. New strategy for
%   encoding onsets. 
% 03/02/20  We are dropping the first TR of each event. Updating to reflect
%   this. 
% 03/09/20  Added 'perfect' output to represent any blocks with perfect
%   accuracy. This will be used to select regressors. 

function onsets = extract_timing_all(subj, study)
%% Parameters and Preallocation
if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

%% Parameters and path
dir_subj  = fullfile(study.path, 'data', subj.name);
dir_reg   = fullfile(dir_subj, 'reg'); 
dir_behav = fullfile(dir_subj, 'behav', 'scan'); 
events_each = 2; 
events_all  = 28; 
 
scan = study.scan; 
dir_behav = fullfile(dir_subj, 'behav', 'scan'); 
onsets = struct(... % Needs to be updated for each experiment. 
    'NOI', [], ... 
    'SIL', [], ... 
    'OR_rate075_clear', [], ... 
    'OR_rate075_snr2',  [], ...
    'OR_rate075_snrn2', [], ...
    'OR_rate125_clear', [], ...
    'OR_rate125_snr2',  [], ...
    'OR_rate125_snrn2', [], ...
    'SR_rate075_clear', [], ...
    'SR_rate075_snr2',  [], ...
    'SR_rate075_snrn2', [], ...
    'SR_rate125_clear', [], ...
    'SR_rate125_snr2',  [], ...
    'SR_rate125_snrn2', []  ...
    ); 

conNames   = fields(onsets);
numCons    = length(conNames); 
scan.order = 0:(scan.epis - 1):(scan.epis - 1)*(events_all + 1) - study.scan.first; 
% dropped first 5 and first TR of each run. 

%% Load data for subject
disp('Now loading behavioral data...')
var = dir(fullfile(dir_behav, '*.xlsx'));
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
    T = readtable(fullfile(var.folder, var.name)); 
else
    error('Unknown error when loading files')
end

disp('Done!')

%% Populate onsets structure with preallocated matrixes
% Needs to be updated per experiment. 
for cc = 1:numCons
    thiscon = conNames{cc}; 
    onsets.(thiscon) = nan(events_each, subj.runs); 
end

%% Grab the correct events
stim = T.Stim;
answerKey = nan(size(stim)); 

male    = cellfun((@(x) contains(x, 'M')), stim); 
female  = cellfun((@(x) contains(x, 'F')), stim); 
noise   = cellfun((@(x) contains(x, 'noise')), stim); 
silence = cellfun((@(x) contains(x, 'silence')), stim); 

answerKey(male)    = 2; 
answerKey(female)  = 1; 
answerKey(noise)   = 3; 
answerKey(silence) = 3;

blocks = ismember(T.BLOCK, 1:subj.runs); 

resp = T.SubjResponse; 
correct = resp == answerKey; 

accuracy = correct; 
dropthese = ~blocks | noise | silence; 
accuracy(dropthese) = []; 
accuracy = mean(accuracy)*100; 
disp(['Subject answered ' num2str(accuracy) '% sentences correctly!'])

correct(noise)   = true; 
correct(silence) = true; 
correct(~blocks) = false; 

%% Get onset info for each event
blocks = 1:subj.runs; 
perfect = false(size(blocks)); 

for bb = blocks
    %% Create regressor
    c = correct(T.BLOCK == bb); % for creating a regressor
    if all(c)
        disp(['Perfect accuracy in block ' num2str(bb) '!!!'])
        perfect(bb) = true; 
    else
        c = repelem(c, study.scan.epis-1); 

        fname = fullfile(dir_reg, ['incorrect_run' num2str(bb) '.txt']); 
        fid = fopen(fname, 'w'); 
        fprintf(fid, '%d\n', ~c); 
        fclose(fid); 
    end
    
    %% Create onsets
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
                idx_cond(1) = idx_cond(1) + 1; % not future proof...
                onsets.NOI(idx_cond(1), bb) = scan.order(idx); 
            elseif strcmp(sentence{ev}, 'S') % silence
                idx_cond(2) = idx_cond(2) + 1;
                onsets.SIL(idx_cond(2), bb) = scan.order(idx); 
            else
                error('unknown stim')
            end

        elseif isstrprop(sentence{ev}, 'digit') % sentence presentation
            % Construct what event it was
            thissent = [syntax{ev} 'R_']; % start with syntax 

            if strcmp(rate{ev}, '1.25')
                tag_r = 'rate125'; 
            elseif strcmp(rate{ev}, '0.75')
                tag_r = 'rate075'; 
            else
                error('unknown stim')
            end

            if isempty(babble{ev})
                tag_b = 'clear'; 
            elseif strcmp(babble{ev}, 'SNR-2')
                tag_b = 'snrn2'; 
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
fname = fullfile(dir_subj, 'onsets_all.mat'); 
if exist(fname, 'file')
    delete(fname)
end

save(fname, 'onsets', 'accuracy', 'perfect')

end
