%% conn_get_regressors
% Run this code while CONN toolbox is open. 
%
% MM/DD/YY - CHANGELOG
% 03/19/20 - Changelog started. Wrote conn_module that gets regressors. Now
%   customizing regressor (dropping 2 columns). 
%
% TODO:
% 1. Write batch that loads all data into CONN (new subjects)
clearvars; 

%% Parameters
dir_batch = pwd;
YA_FST_params

cd ..
dir_data = fullfile(pwd, 'data'); 

%% This code creates the regressors, saves them as dp_swu*.txt
conn_module( 'preprocessing', ...
   'steps', {'functional_regression'}, ...
   'reg_names', {'White Matter','CSF'}, ...
   'reg_dimensions',[2, 2], ...
   'reg_skip', true, ...
   'reg_deriv', [0, 0]);

for ii = 1:length(subj)
    %% Setup for this subject
    thissubj = subj(ii);
    dir_subj = fullfile(dir_data, thissubj.name); 
    dir_func = fullfile(dir_subj, 'FUNCTIONAL'); % location may change
    dir_reg  = fullfile(dir_subj, 'reg'); 
    disp(['Writing regressor for subj ' thissubj.name '...'])
    blocks = 1:thissubj.runs; 
    blocks(thissubj.drop) = []; 
    
    for rr = blocks
        %% Move regressors to regressor directory
        disp(['Block ' num2str(rr)])
        
        fname = fullfile(dir_func, ...
            ['dp_' study.prefix study.runname num2str(rr) '_00007.txt']); 
        if exist(fname, 'file')
            movefile(fname, dir_reg)
        end
        
        %% Load our new regressors
        fname = fullfile(dir_reg, ...
            ['dp_' study.prefix study.runname num2str(rr) '_00007.txt']); 
        T = readtable(fname); 
        
        %% Trim off the first two columns of the data 
        % First two columns appear to be constant and linear terms
        newname = fullfile(dir_reg, ...
            ['dp_trimmed_pc2_' study.prefix study.runname num2str(rr) '_00007.txt']); 
        
        T = T(:, 3:end); 
        writetable(T, newname, 'Delimiter', '\t', 'WriteVariableNames', false)
    end
    
end
