function snr_sd_v2_hybrid_multi_evens(subjs, runs, prefix)
% Runs a signal-to-noise ratio analysis on fMRI data. Check that the
% pathing works for your specific study, and adjust it in the "paths" 
% section of the code. 

% CHANGELOG
% Code runs to check SNR of 06Jul17 test with Lee. Works with 3D nii files -- Matt
% 20 July 17 v5 -- adding a couple of new inputs for use with isss_multiband
% 29 Sep 17 -- Updated input. -- MH

%% Checks inputs
if (~iscell(subjs) || ~iscell(runs) || ~ischar(prefix))
    msg1  = 'Input must be snr_sd_v2({subjects}, {runs}, "prefix"), where subjs and runs are cells, prefix is str'; 
    error(msg1);    
end

numsubjs = length(subjs); 
numruns  = length(runs); 

f = filesep;
warning off

%% Root directory (where you save all your fMRI data)
cd ..
% root_dir = [pwd f study];
root_dir = pwd; % Changed 05/29/18 for supercomputer -- MH

for s = 1:numsubjs % For each subject...
    %% Paths (edit for your lab's conventions)
    disp(' ')
    disp('--------------------------')
    disp(['Starting analysis for subject ', subjs{s}])
    disp('--------------------------')
    
    subj_dir       = [root_dir f 'data_14subjanalysis' f subjs{s}]; % Specific dir for study subjects
    functional_dir = [subj_dir f 'FUNC_GLM'];    % Specific dir for functional data per subject
    
    firstrun = '_00001.nii'; % How do you number your first acquisition? 
 
    %% The Gnarly Stuff (change at your own risk)
    for r = 1:numruns % For each run... 
        
        % Checks pathing
        bold = [prefix, runs{r}];

        if ~isdir(functional_dir)    
            error('no downsample')
        end

        cd(functional_dir);
        files = dir([bold '*.nii']); 

        if isempty(files)
            error('Files not found!')
        end
        
        % bold = dir(['.', f, runs{r}, '_0*.nii']);

    % Determine current preprocessed bold.
    % This will prioritize working with preprocessed files, and defaults to
    % unprocessed files if none are found.
    
%         if   exist([functional_dir f 'swau', runs{r}, firstrun], 'file')
%             bold = ['swau', runs{r}];  
%         elseif exist([functional_dir f 'swr', runs{r}, firstrun], 'file')
%             bold = ['swr', runs{r}];
%         elseif exist([functional_dir f 'swu', runs{r}, firstrun], 'file')
%             bold = ['swu', runs{r}];
%         elseif exist([functional_dir f 'sw',  runs{r}, firstrun], 'file')
%             bold = ['sw',  runs{r}];   
%         elseif exist([functional_dir f 'sra', runs{r}, firstrun], 'file')
%             bold = ['sra', runs{r}];   
%         elseif exist([functional_dir f 'sr',  runs{r}, firstrun], 'file')
%             bold = ['sr',  runs{r}];    
%         elseif exist([functional_dir f 'su',  runs{r}, firstrun], 'file')
%             bold = ['su',  runs{r}];     
%         elseif exist([functional_dir f 'r',  runs{r}, firstrun], 'file')
%             bold = ['r',  runs{r}];     
%         elseif exist([functional_dir f 'u',  runs{r}, firstrun], 'file')
%             bold = ['u',  runs{r}];     
%         else
%             bold = runs{r};   
%         end

        disp(' ')
        disp(['### Found data for ', bold, ' and will start calculating snr'])

        % Select niis for the current run
        niis = dir([functional_dir f bold '_0*.nii']); 
        
        if contains(runs{r}, 'hybrid')
            niis = niis(2:2:end);
        elseif contains(runs{r}, 'multiband')
            out = sort([11:14:248, 12:14:248, 13:14:248, 14:14:248]);
            niis(out) = []; niis = niis(2:2:end);
        end
        
        
        niis_filename = cell(1, length(niis)); 
        for j = 1:length(niis)
            niis_filename{j} = [niis(j).folder f niis(j).name];
        end

        % Preallocate avg and st_tmp using the first scan's dimensions
        files = spm_vol(niis_filename);
        tmp_data = spm_read_vols(files{1}); % ERROR HERE, usually happens if path is wrong, or no files selected
        avg      = zeros(size(tmp_data));
        sd_tmp   = zeros(size(tmp_data));

        % Calculate average signal, sd, and SNR
        % original way in the code
        n = length(files);
        disp(['Calculating average signal for ', bold, ':'])
        for j = 1:n
            data = spm_read_vols(files{j});
            avg = avg + data / n;
            fprintf('.')
        end
        disp(' ')
        
        disp(['Calculating standard deviation for ', bold, ':'])
        for j = 1:n
            data = spm_read_vols(files{j});
            sd_tmp = sd_tmp + (avg - data).^2;  
            fprintf('.')
        end
        disp(' ')

        sd = sqrt(sd_tmp / (n - 1));
        snr = avg ./ sd;

        % Clean up snr varaible
        snr(isnan(snr)) = 0;
        snr(snr>5000)   = 0; % eliminates the absurdly high values that can occur outside the brain

        % Make sure SNR directories exists
        snr_dir = [subj_dir f 'SNR']; 
        unix(['mkdir ' snr_dir]);
        snr_bold_dir = [snr_dir f bold '_EVENS'];
        unix(['mkdir ' snr_bold_dir]);
        unix(['mkdir ' snr_bold_dir f 'AVERAGE']); 
        unix(['mkdir ' snr_bold_dir f 'SD']);
        unix(['mkdir ' snr_bold_dir f 'SNR']);


        % Output files to the SNR dir
        avg_output = files{1};
        sd_output  = files{1};
        snr_output = files{1};
        avg_output.fname = [snr_bold_dir f 'AVERAGE' f bold '_average.nii'];
        sd_output.fname  = [snr_bold_dir f 'SD'      f bold '_sd.nii'];
        snr_output.fname = [snr_bold_dir f 'SNR'     f bold '_snr.nii'];

        % % This way outputs as float, but with a scaling factor
        % spm_write_vol(avg_output, avg); 
        % spm_write_vol(sd_output,  sd);
        % spm_write_vol(snr_output, snr);

        % The following method keeps the scaling factor set to 1
        avg_output = spm_create_vol(avg_output);
        sd_output  = spm_create_vol(sd_output);
        snr_output = spm_create_vol(snr_output);
        disp('### Writing out volumes:')
        for j=1:avg_output.dim(3)
            avg_output = spm_write_plane(avg_output, avg(:,:,j), j);
            sdoutput   = spm_write_plane(sd_output,  sd(:,:,j),  j); %#ok<NASGU>
            snr_output = spm_write_plane(snr_output, snr(:,:,j), j); 
            fprintf('.')
        end
        disp(' ')
        %avg_output=spm_close_vol(avg_output);
        %sd_output=spm_close_vol(sd_output);
        %snr_output=spm_close_vol(snr_output);
        %end snr_sd on preprocessed data for all runs
    %else
    %    msg4=['snr_sd will only operate on the RAW bold images'];
    %    disp(msg4)
    %end

    %Return to root
        cd(root_dir);

    %Removes previous .ps files so as not to over-append
        if exist([snr_dir f subjs{s} '_' runs{r} '_snr.ps'], 'file')
            delete([snr_dir f subjs{s} '_' runs{r} '_snr.ps']);
        end

    % if exist([functional_dir,bold,num2str(numruns),'_0003.nii'])
        subj2 = regexprep(subjs{s},'_','\_'); %small fix for figure printing
        
        %for j=1:numruns
        %load files
        files_avg = spm_vol([snr_bold_dir f 'AVERAGE' f bold '_average.nii']);
        files_sd  = spm_vol([snr_bold_dir f 'SD'      f bold '_sd.nii']);
        files_snr = spm_vol([snr_bold_dir f 'SNR'     f bold '_snr.nii']);
        
        data_avg  = spm_read_vols(files_avg);
        data_sd   = spm_read_vols(files_sd);
        data_snr  = spm_read_vols(files_snr);
        
        %Slices
        slice1_avg  = squeeze(data_avg(:,:,20));
        slice2_avg  = squeeze(data_avg(:,:,24));
        
        slice1_sd   = squeeze(data_sd(:,:,20));
        slice2_sd   = squeeze(data_sd(:,:,24));
        
        axial1_snr  = squeeze(data_snr(:,:,20));
        axial2_snr  = squeeze(data_snr(:,:,24));
        
        sagital_snr = squeeze(data_snr(26,:,:));
        
        coronal_snr = squeeze(data_snr(:,32,:));
        if r == 1
            %Do Figure
            width = 8.5;
            height = 11;
            % Get the screen size in inches
            set(0, 'units', 'inches')
            scrsz = get(0, 'screensize');
            % Calculate the position of the figure
            position = [scrsz(3)/2 - width/2, scrsz(4)/2 - height/2, width, height];
            figure(1), clf;
            h=figure(1);
            set(h,'units','inches')
            % Place the figure
            set(h,'position',position)
            % Do not allow Matlab to resize the figure while printing
            set(h,'paperpositionmode','auto')
            % Set screen and figure units back to pixels
            set(0,'units','pixel')
            set(h,'units','pixel')
            %Set colors
            set(gcf,'color',[1 1 1])
            colormap(hot)
        end
        %Start Plotting
        if r == 1
            titlepos = 0.985;
            subpos = 1;
            reduction = 0;
        else
            titlepos = 0.5;
            subpos = 9;
            reduction = 0.475;
        end
        %Plots
        subplot(4,4,subpos)
            set(gca,'position',[0.05,(0.75-reduction),0.2,0.2])
            imagesc(flipud(slice1_avg'),[10,900])
            hold on
            axis equal
            axis off
            title('Avg1','fontweight','bold','position',[27,0.5])
       subplot(4,4,(subpos+1))
            set(gca,'position',[0.26,(0.75-reduction),0.2,0.2])
            imagesc(flipud(slice2_avg'),[10,900])
            hold on
            axis equal
            axis off
            title('Avg2','fontweight','bold','position',[27,0.5])
       subplot(4,4,(subpos+2))
            set(gca,'position',[0.54,(0.75-reduction),0.2,0.2])
            imagesc(flipud(slice1_sd'),[2,20])
            hold on
            axis equal
            axis off
            title('SD1','fontweight','bold','position',[27,0.5])
       subplot(4,4,(subpos+3))
            set(gca,'position',[0.75,(0.75-reduction),0.2,0.2])
            imagesc(flipud(slice2_sd'),[2,20])
            hold on
            axis equal
            axis off
            title('SD2','fontweight','bold','position',[27,0.5])
       subplot(4,4,(subpos+4))
            set(gca,'position',[0.05,(0.53-reduction),0.2,0.2])
            imagesc(flipud(axial1_snr'),[10,350])
            hold on
            axis equal
            axis off
            title('SnR1','fontweight','bold','position',[27,0.5])
       subplot(4,4,(subpos+5))
            set(gca,'position',[0.26,(0.53-reduction),0.2,0.2])
            imagesc(flipud(axial2_snr'),[10,350])
            hold on
            axis equal
            axis off
            title('SnR2','fontweight','bold','position',[27,0.5])
       subplot(4,4,(subpos+6))
            set(gca,'position',[0.54,(0.53-reduction),0.2,0.2])
            imagesc(flipud(sagital_snr'),[10,350])
            hold on
            axis equal
            axis off
            title('SnR Sagital','fontweight','bold','position',[33,0.5])
       subplot(4,4,(subpos+7))
            set(gca,'position',[0.75,(0.53-reduction),0.2,0.2])
            imagesc(flipud(coronal_snr'),[10,350])
            hold on
            axis equal
            axis off
            title('SnR Coronal','fontweight','bold','position',[27,0.5])
       %Title
       ttl = ['SNR\_SD: ' subj2 ' ' runs{r}];
       tax = axes('Position',[0.01,titlepos,1,1]);
       tmp= text(0,0,ttl);
       set(tax,'xlim',[0,1],'ylim',[0,1])
       set(tmp,'FontSize',16,'HorizontalAlignment','left','FontWeight','bold')
       axis off
       %Plot checks
       if r == numruns % if last run
               %Print, close and return
               prnstr = ['print -dpsc2 -painters -append ',[subj_dir f 'SNR' f subjs{s} '_snr.ps']];
               eval(prnstr);
               disp(['SNR_SD output printed to ',f subjs{s} f 'SNR' f subjs{s},'_snr.ps'])
               close (1);
       else % If not the last figure
           %print and clear figure
           prnstr = ['print -dpsc2 -painters -append ',[subj_dir f 'SNR' f subjs{s} '_snr.ps']];
           eval(prnstr);
           pause(0.5);
           clf;
       end     
    end
end
end
