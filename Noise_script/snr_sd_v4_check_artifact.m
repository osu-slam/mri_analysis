function snr_sd_v4_check_artifact(subj, study)
% Runs a signal-to-noise ratio analysis on fMRI data. Check that the
% pathing works for your specific study, and adjust it in the "paths" 
% section of the code. 

% CHANGELOG
% Code runs to check SNR of 06Jul17 test with Lee. Works with 3D nii files -- Matt
% 20 July 17 v5 -- adding a couple of new inputs for use with isss_multiband
% 29 Sep 17 -- Updated input. -- MH
% 20 Sep 19 -- SNR analysis: this version drops out the "silent" period of
% the multi scans. Figures removed. -- MH
% 02/11/20 -- Made universal. 
% 02/13/20 -- Added PDF back in
% 02/20/20 -- Checking artifact that Xiangrui found. MJ

%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

f = filesep;

%% Root directory (where you save all your fMRI data)
dir_root = study.path;
dir_subj = fullfile(dir_root, 'data', subj.name); % Specific dir for study subjects

%% Paths (edit for your lab's conventions)
disp(' ')
disp('--------------------------')
disp(['Starting analysis for subject ', subj.name])
disp('--------------------------')

%% The Gnarly Stuff 
for rr = 1:subj.runs % For each run... 
    thisrun = [lower(study.prefix), study.runname num2str(rr)]; 
    data_subj = fullfile(dir_subj, 'FUNC_GLM', [thisrun '_0*.nii']); 

    niis = dir(data_subj); 

    if isempty(niis)
        warning(['Files ' thisrun ' not found!'])
    else
        disp(' ')
        disp(['### Found data for ' thisrun ' and will start calculating snr'])

        % Select niis for the current run
        niis_filename = cell(length(niis), 1); 
        for jj = 1:length(niis)
            niis_filename{jj} = fullfile(niis(jj).folder, niis(jj).name);
        end

        % Preallocate avg and st_tmp using the first scan's dimensions
        V_files = spm_vol(niis_filename);
        % ERROR HERE, usually happens if path is wrong, or no files selected
        avg    = zeros([V_files{1}.dim, 4]);
        sd_tmp = zeros([V_files{1}.dim, 4]);
        num_scans = ceil((length(V_files)-1*ones(1, 4))./4); 
        if mod(length(V_files)-1, 4) ~= 0
            num_scans(mod(length(V_files) - 1, 4)+1) = num_scans(1) - 1; 
        end
        
        % Calculate average signal, sd, and SNR 
        disp(['Calculating average signal for ', thisrun, ':'])
        idx = 1; 
        for jj = 1:length(V_files)
            data = spm_read_vols(V_files{jj});
            fprintf('.')
            if jj ~= 1 % drop first file
                avg(:, :, :, idx) = avg(:, :, :, idx) + (data / num_scans(idx));
                
                if idx == 4
                    idx = 0; 
                end
            
                idx = idx + 1; 
            end
            
        end
        disp(' ')

        disp(['Calculating standard deviation for ', thisrun, ':'])
        idx = 1; 
        for jj = 1:length(V_files)
            data = spm_read_vols(V_files{jj});
            fprintf('.')
            if jj ~= 1 % drop first file
                sd_tmp(:, :, :, idx) = ...
                    sd_tmp(:, :, :, idx) + ... 
                    (avg(:, :, :, idx) - data).^2; 
                if idx == 4
                    idx = 0; 
                end
            
                idx = idx + 1; 
            end
            
        end
        disp(' ')

        sd = sqrt(sd_tmp); 
        for ss = 1:4
            sd(:, :, :, ss) = sd(:, :, :, ss) / (num_scans(ss) - 1); 
        end
        
        snr = avg ./ sd;

        % Clean up snr varaible
        snr(isnan(snr)) = 0;
        snr(snr>5000)   = 0; 
        % trim absurdly high values that can occur outside the brain

        % Make sure SNR directories exists
        dir_snr = [dir_subj f 'SNR']; 
        unix(['mkdir ' dir_snr]);
        dir_snr_bold = [dir_snr f thisrun f 'artifact'];
        unix(['mkdir ' dir_snr_bold]);
        unix(['mkdir ' dir_snr_bold f 'AVERAGE']); 
        unix(['mkdir ' dir_snr_bold f 'SD']);
        unix(['mkdir ' dir_snr_bold f 'SNR']);

        % Output files to the SNR dir
        avg_output = V_files{1};
        sd_output  = V_files{1};
        snr_output = V_files{1};
        for ss = 1:4
            avg_output.fname = [dir_snr_bold f 'AVERAGE' f thisrun '_average_TR' num2str(ss) '.nii'];
            sd_output.fname  = [dir_snr_bold f 'SD'      f thisrun '_sd_TR' num2str(ss) '.nii'];
            snr_output.fname = [dir_snr_bold f 'SNR'     f thisrun '_snr_TR' num2str(ss) '.nii'];

            % The following method keeps the scaling factor set to 1
            avg_output = spm_create_vol(avg_output);
            sd_output  = spm_create_vol(sd_output);
            snr_output = spm_create_vol(snr_output);
            disp('### Writing out volumes:')
            for jj=1:avg_output.dim(3)
                avg_output = spm_write_plane(avg_output, avg(:,:,jj, ss), jj);
                sd_output  = spm_write_plane(sd_output,  sd(:,:,jj, ss),  jj); %#ok<NASGU>
                snr_output = spm_write_plane(snr_output, snr(:,:,jj, ss), jj); 
                fprintf('.')
            end
            
            disp(' ')
        end
        
        cd(dir_root);   
    end

end
end
