function snr_sd_v3_hybrid_multi_downsample(subjs, runs, prefix)
% Runs a signal-to-noise ratio analysis on fMRI data. Check that the
% pathing works for your specific study, and adjust it in the "paths" 
% section of the code. 

% CHANGELOG
% Code runs to check SNR of 06Jul17 test with Lee. Works with 3D nii files -- Matt
% 20 July 17 v5 -- adding a couple of new inputs for use with isss_multiband
% 29 Sep 17 -- Updated input. -- MH
% 20 Sep 19 -- SNR analysis: this version drops out the "silent" period of
% the multi scans. Figures removed. -- MH

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
dir_root = pwd; % Changed 05/29/18 for supercomputer -- MH
dir_subj  = fullfile(dir_root, 'data_14subjanalysis', subjs); % Specific dir for study subjects

for ss = 1:numsubjs % For each subject...
    %% Paths (edit for your lab's conventions)
    disp(' ')
    disp('--------------------------')
    disp(['Starting analysis for subject ', subjs{ss}])
    disp('--------------------------')
    
%     firstrun = '_00001.nii'; % How do you number your first acquisition? 
 
    %% The Gnarly Stuff 
    for rr = 1:numruns % For each run... 
        % Checks pathing
        data_subj = fullfile(dir_subj{ss}, 'FUNC_GLM', ...
            [lower(prefix), runs{rr} '_0*.nii']); 
            % Specific dir for functional data per subject
            
        niis = dir(data_subj); 

        if isempty(niis)
            warning(['Files ' lower(prefix), runs{rr}, ' not found!'])
        else
            bold = [lower(prefix), runs{rr}]; 
            disp(' ')
            disp(['### Found data for ', lower(prefix), runs{rr}, ' and will start calculating snr'])

            % Select niis for the current run
            niis_filename = cell(length(niis), 1); 
            for jj = 1:length(niis)
                niis_filename{jj} = fullfile(niis(jj).folder, niis(jj).name);
            end
            
            %% Remove silent scans from files list
            if contains(runs{rr}, 'multiband')
                scans_silent = sort([11:14:248, 12:14:248, 13:14:248, 14:14:248], 'descend');
                for ii = scans_silent % remove 68 "silent" slices
                    niis_filename(ii) = [];
                end

                disp('removed silent events!')
            end
            
            %% downsample (just take evens for now)
            if contains(runs{rr}, {'hybrid', 'multiband'})
                scans_silent = sort(1:2:length(niis_filename), 'descend'); 
                % drop odd
                for ii = scans_silent % downsample
                    niis_filename(ii) = [];
                end

                disp('downsampled!')
            end
            
            %% Carry on!
            % Preallocate avg and st_tmp using the first scan's dimensions
            V_files = spm_vol(niis_filename);
            % ERROR HERE, usually happens if path is wrong, or no files selected
            avg      = zeros(V_files{1}.dim);
            sd_tmp   = zeros(V_files{1}.dim);

            % Calculate average signal, sd, and SNR 
            disp(['Calculating average signal for ', bold, ':'])
            for jj = 1:length(V_files)
                data = spm_read_vols(V_files{jj});
                avg = avg + (data / length(V_files));
                fprintf('.')
            end
            disp(' ')

            disp(['Calculating standard deviation for ', bold, ':'])
            for jj = 1:length(V_files)
                data = spm_read_vols(V_files{jj});
                sd_tmp = sd_tmp + (avg - data).^2;  
                fprintf('.')
            end
            disp(' ')

            sd = sqrt(sd_tmp / (length(V_files) - 1));
            snr = avg ./ sd;

            % Clean up snr varaible
            snr(isnan(snr)) = 0;
            snr(snr>5000)   = 0; 
            % trim absurdly high values that can occur outside the brain

            % Make sure SNR directories exists
            dir_snr = [dir_subj{ss} f 'SNR']; 
            unix(['mkdir ' dir_snr]);
            dir_snr_bold = [dir_snr f bold '_nosilent_downsample'];
            unix(['mkdir ' dir_snr_bold]);
            unix(['mkdir ' dir_snr_bold f 'AVERAGE']); 
            unix(['mkdir ' dir_snr_bold f 'SD']);
            unix(['mkdir ' dir_snr_bold f 'SNR']);

            % Output files to the SNR dir
            avg_output = V_files{1};
            sd_output  = V_files{1};
            snr_output = V_files{1};
            avg_output.fname = [dir_snr_bold f 'AVERAGE' f bold '_average.nii'];
            sd_output.fname  = [dir_snr_bold f 'SD'      f bold '_sd.nii'];
            snr_output.fname = [dir_snr_bold f 'SNR'     f bold '_snr.nii'];
            
%             avg_output.pinfo = ;
%             sd_output.pinfo  = [];
%             snr_output.pinfo = []; 
%             % % This way outputs as float, but with a scaling factor
%             disp('### Writing out volumes:')
%             spm_write_vol(avg_output, avg); 
%             spm_write_vol(sd_output,  sd);
%             spm_write_vol(snr_output, snr);

            % The following method keeps the scaling factor set to 1
            avg_output = spm_create_vol(avg_output);
            sd_output  = spm_create_vol(sd_output);
            snr_output = spm_create_vol(snr_output);
            disp('### Writing out volumes:')
            for jj=1:avg_output.dim(3)
                avg_output = spm_write_plane(avg_output, avg(:,:,jj), jj);
                sdoutput   = spm_write_plane(sd_output,  sd(:,:,jj),  jj); %#ok<NASGU>
                snr_output = spm_write_plane(snr_output, snr(:,:,jj), jj); 
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
            cd(dir_root);

%         %Removes previous .ps files so as not to over-append
%             if exist([dir_snr f subjs{ss} '_' runs{rr} '_snr.ps'], 'file')
%                 delete([dir_snr f subjs{ss} '_' runs{rr} '_snr.ps']);
%             end
% 
%         % if exist([functional_dir,bold,num2str(numruns),'_0003.nii'])
%             subj2 = regexprep(subjs{ss},'_','\_'); %small fix for figure printing
% 
%             %for j=1:numruns
%             %load files
%             files_avg = spm_vol([dir_snr_bold f 'AVERAGE' f bold '_average.nii']);
%             files_sd  = spm_vol([dir_snr_bold f 'SD'      f bold '_sd.nii']);
%             files_snr = spm_vol([dir_snr_bold f 'SNR'     f bold '_snr.nii']);
% 
%             data_avg  = spm_read_vols(files_avg);
%             data_sd   = spm_read_vols(files_sd);
%             data_snr  = spm_read_vols(files_snr);
% 
%             %Slices
%             slice1_avg  = squeeze(data_avg(:,:,20));
%             slice2_avg  = squeeze(data_avg(:,:,24));
% 
%             slice1_sd   = squeeze(data_sd(:,:,20));
%             slice2_sd   = squeeze(data_sd(:,:,24));
% 
%             axial1_snr  = squeeze(data_snr(:,:,20));
%             axial2_snr  = squeeze(data_snr(:,:,24));
% 
%             sagital_snr = squeeze(data_snr(26,:,:));
% 
%             coronal_snr = squeeze(data_snr(:,32,:));
%             if rr == 1
%                 %Do Figure
%                 width = 8.5;
%                 height = 11;
%                 % Get the screen size in inches
%                 set(0, 'units', 'inches')
%                 scrsz = get(0, 'screensize');
%                 % Calculate the position of the figure
%                 position = [scrsz(3)/2 - width/2, scrsz(4)/2 - height/2, width, height];
%                 figure(1), clf;
%                 h=figure(1);
%                 set(h,'units','inches')
%                 % Place the figure
%                 set(h,'position',position)
%                 % Do not allow Matlab to resize the figure while printing
%                 set(h,'paperpositionmode','auto')
%                 % Set screen and figure units back to pixels
%                 set(0,'units','pixel')
%                 set(h,'units','pixel')
%                 %Set colors
%                 set(gcf,'color',[1 1 1])
%                 colormap(hot)
%             end
%             %Start Plotting
%             if rr == 1
%                 titlepos = 0.985;
%                 subpos = 1;
%                 reduction = 0;
%             else
%                 titlepos = 0.5;
%                 subpos = 9;
%                 reduction = 0.475;
%             end
%             %Plots
%             subplot(4,4,subpos)
%                 set(gca,'position',[0.05,(0.75-reduction),0.2,0.2])
%                 imagesc(flipud(slice1_avg'),[10,900])
%                 hold on
%                 axis equal
%                 axis off
%                 title('Avg1','fontweight','bold','position',[27,0.5])
%            subplot(4,4,(subpos+1))
%                 set(gca,'position',[0.26,(0.75-reduction),0.2,0.2])
%                 imagesc(flipud(slice2_avg'),[10,900])
%                 hold on
%                 axis equal
%                 axis off
%                 title('Avg2','fontweight','bold','position',[27,0.5])
%            subplot(4,4,(subpos+2))
%                 set(gca,'position',[0.54,(0.75-reduction),0.2,0.2])
%                 imagesc(flipud(slice1_sd'),[2,20])
%                 hold on
%                 axis equal
%                 axis off
%                 title('SD1','fontweight','bold','position',[27,0.5])
%            subplot(4,4,(subpos+3))
%                 set(gca,'position',[0.75,(0.75-reduction),0.2,0.2])
%                 imagesc(flipud(slice2_sd'),[2,20])
%                 hold on
%                 axis equal
%                 axis off
%                 title('SD2','fontweight','bold','position',[27,0.5])
%            subplot(4,4,(subpos+4))
%                 set(gca,'position',[0.05,(0.53-reduction),0.2,0.2])
%                 imagesc(flipud(axial1_snr'),[10,350])
%                 hold on
%                 axis equal
%                 axis off
%                 title('SnR1','fontweight','bold','position',[27,0.5])
%            subplot(4,4,(subpos+5))
%                 set(gca,'position',[0.26,(0.53-reduction),0.2,0.2])
%                 imagesc(flipud(axial2_snr'),[10,350])
%                 hold on
%                 axis equal
%                 axis off
%                 title('SnR2','fontweight','bold','position',[27,0.5])
%            subplot(4,4,(subpos+6))
%                 set(gca,'position',[0.54,(0.53-reduction),0.2,0.2])
%                 imagesc(flipud(sagital_snr'),[10,350])
%                 hold on
%                 axis equal
%                 axis off
%                 title('SnR Sagital','fontweight','bold','position',[33,0.5])
%            subplot(4,4,(subpos+7))
%                 set(gca,'position',[0.75,(0.53-reduction),0.2,0.2])
%                 imagesc(flipud(coronal_snr'),[10,350])
%                 hold on
%                 axis equal
%                 axis off
%                 title('SnR Coronal','fontweight','bold','position',[27,0.5])
%            %Title
%            ttl = ['SNR\_SD: ' subj2 ' ' runs{rr}];
%            tax = axes('Position',[0.01,titlepos,1,1]);
%            tmp= text(0,0,ttl);
%            set(tax,'xlim',[0,1],'ylim',[0,1])
%            set(tmp,'FontSize',16,'HorizontalAlignment','left','FontWeight','bold')
%            axis off
%            %Plot checks
%            if rr == numruns % if last run
%                    %Print, close and return
%                    prnstr = ['print -dpsc2 -painters -append ',[dir_subj{ss} f 'SNR' f subjs{ss} '_snr.ps']];
%                    eval(prnstr);
%                    disp(['SNR_SD output printed to ',f subjs{ss} f 'SNR' f subjs{ss},'_snr.ps'])
%                    close (1);
%            else % If not the last figure
%                %print and clear figure
%                prnstr = ['print -dpsc2 -painters -append ',[dir_subj{ss} f 'SNR' f subjs{ss} '_snr.ps']];
%                eval(prnstr);
%                pause(0.5);
%                clf;
%            end     
        end

    end
end
end
