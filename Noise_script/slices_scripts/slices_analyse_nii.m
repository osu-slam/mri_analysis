
% ANALYSE YOUR FMRI DATA FOR SLICE NOISE
% this it the first script which must be run on all data sets
%
%  written by Antonia Hamilton, Dartmouth College, 2007

clear all, close all, clc;
spm_defaults  %% assume spm2

disp('Welcome to The Slices Noise Analysis Progam')
disp('  ')

try
  load slices_def
catch
  disp('Could not load slices_def.mat')
  disp('Please run the slices_defaults program')
  disp('and make sure you are in the right directory')
  return
end

%templatedir = [pwd,filesep,'noise_templates'];
%if(~exist([templatedir,filesep,'xEPI.img']))
%  disp(['Cannot find EPI template image in ',templatedir])
%  return
%end
%if(~exist([templatedir,filesep,'xwhole_brain.img']))
%  disp(['Cannot find mask template image in ',templatedir])
%  return
%end

disp('Please select the first image of the first block')

% Here is the first line of code I edited -- Matt
tmp = spm_select(1,'.nii','Select the first image of the first block')
[pth,nam,ext] = fileparts(tmp);
path=[pth '\'];

% Here is the second line of code I edited -- Matt
disp('Assuming all images are called bold*.nii')
nses = input('How many sessions? ')

outname = input('Enter the filename for saving stuff (no extension) ','s')

%Turn off annoying warnings
warning off MATLAB:divideByZero;
warning off

for ses = 1:nses
  %P = spm_select('List',pth,['bold',num2str(ses),'*.img']);
  % Here is the third line -- Matt
  image_files_without_path = spm_select('List',pth,['^bold' num2str(ses) '.*\.nii$']); 
  copies_of_directory_path = repmat(path,240,1);
  image_files_with_path=[copies_of_directory_path image_files_without_path];
  V = spm_vol(image_files_with_path);
  nslice = V(1).dim(3);
  nscan = length(V);
  
  %% realign data
  if (1)
    disp('realigning the files now')
    spm_realign(V);
  else
    disp('skipping realignment')
  end
  [pth,nam,ext] = fileparts(V(1).fname);
  
  Rnames = fullfile(pth,['rp_',nam,'.txt'])
  ra = load(Rnames);
  
  disp(['Realignment done for session ',num2str(ses),' now reslicing'])
  
  %% reslice data so we can get the mean
  spm_reslice(V);
  %disp('skipping reslice')
  
  %% load in realigned data
  % Fourth line -- Matt
  Pr = spm_select('Files',pth,['rbold',num2str(ses),'*.nii']);
  Vr = spm_vol(Pr);
  % Fifth line -- Matt
  Vo = spm_vol(spm_select('Files',pth,['meanbold',num2str(ses),'*.nii']));

  %% get characteristics of mean image
  vox_size = diag(Vo.mat); 
  vox_size = vox_size(1:3);

  d = Vo.dim(1:3); 
  % corners in voxel-space
  c = [ 1    1    1    1 
        1    1    d(3) 1 
        1    d(2) 1    1 
        1    d(2) d(3) 1
        d(1) 1    1    1 
        d(1) 1    d(3) 1
        d(1) d(2) 1    1 
        d(1) d(2) d(3) 1 ]'; 
  % corners in world-space
  tc = Vo.mat(1:3,1:4)*c; 
  % reflect in x if required
  if spm_flip_analyze_images; tc(1,:) = -tc(1,:); end; 
  
  % bounding box (world) min and max
  mn = min(tc,[],2)'; 
  mx = max(tc,[],2)'; 
  bb = [mn; mx];
  bb = bb-0.1;

  %% copy mask locally
  % Sixth line -- Matt
  eval(['!cp ',sdef.brain_template,' ',...
        pth,filesep,'mask',num2str(ses),'.nii'])
  eval(['!cp ',sdef.brain_template(1:end-4),'.hdr ',...
        pth,filesep,'mask',num2str(ses),'.hdr'])
%  eval(['!cp ',templatedir,filesep,'xEPI.img ',... 
%        pth,filesep,'EPI',num2str(ses),'.img'])
%  eval(['!cp ',templatedir,filesep,'xEPI.hdr ',... 
%        pth,filesep,'EPI',num2str(ses),'.hdr']) 
   
  %% normalise EPI template to mean BOLD image
  
  %% load up template
  Ve = spm_vol(sdef.epi_template);
  
  defaults.normalise.write.vox = vox_size';   %% voxel size 
  defaults.normalise.write.bb =  bb;  %% my bounding box
  defaults.normalise.write.interp = 7; 
  defaults.normalise.write.wrap = [0 0 0];
  
  % Normalise template to mean EPI image
  matname  = [pth,filesep,'norm_bold',num2str(ses),'.mat'];
  params   = spm_normalise(Vo,Ve,matname,'','',defaults.normalise.estimate);  
 
  %% apply normalisation to mask image
  % Sixth line -- Matt
  Vm = spm_vol(spm_select('Files',pth,['mask',num2str(ses),'.nii'])) 
  Vn = spm_write_sn(Vm,params,defaults.normalise.write);
  spm_write_vol(Vn,Vn.dat);
  
  disp(['Mask created for sesssion ',num2str(ses)])

  %-------------------------------------------------------------

 %% work out good neighbours
 neigh = repmat(-5:2:5,nscan,1)+repmat(1:nscan,6,1)';
 neigh(neigh<1) = NaN;
 neigh(neigh>nscan) = NaN;
 
 for kk=1:nscan
        
   %% calculate good neighbours
   df = abs(ra-repmat(ra(kk,:),nscan,1));
   trans = sum(df(:,1:3),2);
   rot = sum(df(:,4:6),2);
   
   good = find(trans<1 & (rot*180/pi)<1);
 
   overlap = intersect(good,neigh(kk,:));
   
   if(length(overlap)<3)  %% really bad scans
     neigh(neigh==kk) = NaN;     
     neigh(kk,:) = NaN;
   elseif(length(overlap<6))  %% moderate bad scans
     missing = setdiff(neigh(kk,:),good);
     for i=1:length(missing)
       mind = find(neigh(kk,:)==missing(i));
       neigh(kk,mind) = NaN;
     end
   end
 end
 
  %% now look for noise within mask
  % Seventh line -- Matt
 M = spm_vol(spm_select('Files',pth,['wmask',num2str(ses),'.nii']));
 M.dat = spm_read_vols(M);
 
 disp('Now searching every file for noise')
 
 for kk=1:nscan  %% for every input file
   
   gv(kk) = spm_global(V(kk));
   
   [dat1,loc] = spm_read_vols(V(kk),1);  %% read with zero masking   
   nn = neigh(kk,:);
   
   for jj=1:length(nn)  %% read some neighbours
     
     if(isfinite(nn(jj)))  %% for good neighbours
       jind = nn(jj);
       
       [dat2,loc] = spm_read_vols(V(jind),1);  %% read with zero masking
       
       for i=1:nslice
         slice1 = squeeze(dat1(:,:,i));
         slice2 = squeeze(dat2(:,:,i));
         msk = squeeze(M.dat(:,:,i));
         
         df = (slice1-slice2).*(msk>0.5);        
         scan_noise(jj,i) = nanmean(abs(df(msk>0.5)));
       end
     else  %% for bad neighbours
       scan_noise(jj,1:nslice) = NaN;
     end
   end
   noise(:,kk) = nanmean(scan_noise)';
   
   kk
 end
 
 %% save all important info
 a(ses).P = P;
 a(ses).noise = noise;
 a(ses).mask = M;
 a(ses).ra = ra;
 a(ses).neigh = neigh;
 
 slices_summary(a(ses),[outname,'.ps'])
 
 pause(0.1)
 
end

save(outname,'a')
% eval(['print -dpsc2 -painters ',outname,'.ps'])
disp(['Noise analysis saved to ',outname,'.mat and printed to ',outname,'.ps'])

return
