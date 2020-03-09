%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this program sets up the defaults for the Slices scripts
%% you will need to edit the file as marked below for your setup

%% written by Antonia Hamilton, Dartmouth College, 2007


clear all, close all
disp(' ')
disp('This is the slices defaults script')
disp('You will want to edit this in your favourite text editor')
disp('to set the right paths and thresholds')

%%% Edit the next 2 lines to change the EPI template location
sdef.epi_template = [pwd, '\noise_templates\xEPI.img'];
disp(' ')
disp('The current EPI template location is:')
disp(sdef.epi_template)
if(exist(sdef.epi_template)==2)
  disp('This file is OK')
else
  disp('Cannot find the file, please edit slices_defaults.m')
end

%% edit the next 2 lines to change the brain template location
sdef.brain_template = [pwd, '\noise_templates\xwhole_brain.img'];
disp(' ')
disp('The whole brain template location is')
disp(sdef.brain_template)
if(exist(sdef.brain_template)==2)
  disp('This file is OK')
else
  disp('Cannot find the file, please edit slices_defaults.m')
end

sdef.th = 25;   %% edit this number to change your threshold
disp(' ')
disp(['The current noise threshold is:  ',num2str(sdef.th)])


sdef.wth = 30;   %% edit this number to change your weaker threshold
disp(['The current weaker threshold (not really used) is:  ',num2str(sdef.wth)])
disp(' ')


save slices_def sdef
disp('Slices defaults have been saved in:')
disp([pwd,filesep,'slices_def.mat'])
disp(' ')
disp('Now you are ready to use the slices scripts')





