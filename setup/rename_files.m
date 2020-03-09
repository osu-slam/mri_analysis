%% rename_files
% Run this if there was a terminated run that left unseemly, strange names.

%% Actual code
cd(fullfile(dir_subj, 'FUNCTIONAL'))
target = subj(ii).rename; 

for tt = 1:length(target)
    % Get sessions
    files = dir(['*' target{tt} '*']); files = {files(:).name}; 
    split = cellfun((@(x) strsplit(x, '_')), files, 'UniformOutput', false); 
    sessnum = cellfun((@(x) x{3}), split, 'UniformOutput', false); 
    sessnum = unique(sessnum); 
    
    num_sessnum = zeros(size(sessnum)); 
    for ss = 1:length(sessnum)
        num_sessnum(ss) = sum(contains(files, sessnum{ss})); 
    end
    
    [~, idx] = max(num_sessnum); 
    disp(['Keeping ' sessnum{idx} ' because it has the most images!'])
    
    keep = dir(['*' sessnum{idx} '*']); keep = {keep(:).name}; 
    dest = cellfun((@(x) strsplit(x, '_')), keep, 'UniformOutput', false); 
    dest = cellfun((@(x) strjoin(x([1, 2, 4]), '_')), dest, 'UniformOutput', false); 
    for ff = 1:length(keep)
        movefile(keep{ff}, dest{ff}); 
    end
    
    other = size(sessnum); 
    other(idx) = []; 
    drop = dir(['*' sessnum{other} '*']); drop = {drop(:).name}; 
    for ff = 1:length(drop)
        delete(drop{ff}); 
    end
    
end
