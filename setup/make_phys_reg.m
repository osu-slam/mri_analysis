%% make_physio_regressor(run)
% Does exactly what it says on the tin. Save your zip file holding the
% regressor in the \PHYSIO\run folder, and sit back to watch the magic
% happen. Make sure to update the number of images!!!

% CHANGELOG
% 09/11/17  File inception
% 09/15/17  Updated file to follow current conventions of isss_multi_params
% 02/08/20  Forked for YA_FST
% 02/29/20  Adding input for number of scans as part of subj struct, 
%   updated to use subj input. Revamped for second analysis. 
% 03/20/20  Updated to match new (drop first, etc.) design

function make_phys_reg(subj, study)
%% Check input
if ~isstruct(subj) || (length(subj) ~= 1)
    error('Input ("subj") where subj is ONE structure')
end

if ~isstruct(study)
    error('Input ("study") where study is structure')
end

%% Pathing and params
dir_subj  = fullfile(study.path, 'data', subj.name); 

% dir_nii  = fullfile(dir_subj, 'nii'); 
% dir_func = fullfile(dir_subj, 'FUNCTIONAL'); 
dir_physio = fullfile(dir_subj, 'PHYSIO'); 
dir_physio_reg = fullfile(dir_physio, 'reg');
dir_reg = fullfile(dir_subj, 'reg'); 
scan = study.scan; 
nscans = scan.first + (scan.epis + scan.silence) * scan.events; 
nscans_drop = (scan.epis-1) * scan.events; 

%% Read in raw data
niis = dir(fullfile(dir_physio_reg, '*.nii'));
hdrs = dir(fullfile(dir_physio_reg, '*.mat'));

V = cell(1, length(niis));
P = cell(1, length(niis));
for n = 1:length(niis)
    file = fullfile(hdrs(n).folder, hdrs(n).name); 
    load(file)
    P{n} = h.PulseRespiratoryRegressors.ProtocolName; 
    
    file = fullfile(niis(n).folder, niis(n).name); 
    hdr = spm_vol(file); 
    V{n} = spm_read_vols(hdr); 
    V{n} = reshape(V{n}, size(V{n}, 3), size(V{n}, 2));
end

skipRuns = []; 
for ii = 1:length(V)
    if length(V{ii}) ~= nscans
        skipRuns = [skipRuns ii];
        disp('Found a terminated run. Skipping...')
    end
end

if ~isempty(skipRuns)
    for rr = skipRuns
        P(rr) = [];
        V(rr) = [];
    end
end

disp(['Found ' num2str(length(P)) ' files total.'])

lengths = cellfun(@length, V);     
if length(unique(lengths)) ~= 1
    error('Inconsistent number of elements!')
end

for rr = 1:length(P)
    %% Prepare to make regressors
    skipToEnd = 0;
    firstScans = false(scan.first, 1); 
    % Not modeling first four scans. And remember the sinister TR bug...
    template = vertcat(false(scan.silence+1, 1), true(scan.epis-1, 1));
    % All models skip the first TR because of artifact
    extract = vertcat(firstScans, repmat(template, [scan.events, 1])); 
    
%     if lengths(1) ~= 125 % if we have one too few images...
%         extract(end) = []; 
%     end

    try
        reg = []; 
        for zz = 1:subj.runs
            temp = V{rr}; 
            
            for cc = 1:size(temp, 2)
                reg = [reg, temp(extract, cc)];                
            end
            
        end
        
        assert(size(reg, 1) == nscans_drop)
    catch
        disp(['Something happened during ' P{rr} ' and the regressor was not created.'])
        skipToEnd = 1;
    end

    if ~skipToEnd
        disp(['Saving ' P{rr}])

        filename1 = fullfile(dir_reg, ['physio_full_' P{rr} '_00007.txt']);
        fid = fopen(filename1, 'w');
        for line = 1:length(reg)
            eline = '\t';
            for ii = 1:8
                if reg(line, ii) < 0
                    eline = [eline, '%e   '];
                elseif reg(line, ii) >= 0
                    eline = [eline, ' %e   '];
                end
            end
            fprintf(fid, eline, reg(line, :));
            fprintf(fid, '\n');
        end
        fclose(fid); 

        filename2 = fullfile(dir_reg, ['physio_1st_' P{rr} '_00007.txt']);
        fid = fopen(filename2, 'w');
        for line = 1:length(reg)
            eline = '\t';
            for ii = [1 2 5 6]
                if reg(line, ii) < 0
                    eline = [eline, '%e   '];
                elseif reg(line, ii) >= 0
                    eline = [eline, ' %e   '];
                end
            end
            fprintf(fid, eline, reg(line, [1 2 5 6]));
            fprintf(fid, '\n');
        end
        fclose(fid); 
    end

end

end
