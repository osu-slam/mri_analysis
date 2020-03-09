%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% examine samples of good and bad slices from a dataset
%% you must run slices_analyse first so you have some data

%% written by Antonia Hamilton, Dartmouth College, 2007


clear all, close all

disp('Welcome the the noise examination script')
disp(' ')
disp(' (for all options, please type in lower case, then hit RETURN)')

spm_defaults  %% assume spm2

try
  load slices_def
catch
  disp('Could not load slices_def.mat')
  disp('Please run the slices_defaults program')
  disp('and make sure you are in the right directory')
  return
end


disp('Pick the saved noise file to analyse')
th = sdef.th;
wth = sdef.wth;
disp(['Current noise threshold is ',num2str(th)])
disp(['Weak noise threshold is ',num2str(wth)])

inname = spm_get(1,'*.mat','Pick the saved noise file to analyse')

%% sort out variables
load(inname)

ses = input(['Which session (1 to ',num2str(length(a)),') ']);
nses = length(a);
nslice = size(a(1).noise,1);
nscan = size(a(1).noise,2);

P = a(ses).P;
noise = a(ses).noise;
mask = a(ses).mask;
ra = a(ses).ra;
neigh = a(ses).neigh;

V = spm_vol(P);

%% plot summary of this session
slices_summary(a(ses),th,wth,'')

%% decide on display mode
disp(' ')
disp('What set of slices do you want to see?')
disp('Random selection of good & bad = R')
disp('Start with the worst = W')
disp('Random selection of slices near threshold = T')

opt = input(' option? ','s')

opt_ok = 0;
while(opt_ok < 1)
    if (opt == 'R' | opt == 'r' | opt == 'W' | opt == 'w' | opt == 'T' | opt == 't')
        opt_ok = 1;
    else
        disp('Incorrect selection, please select a valid option (r,w,t,R,W,T)')
        opt = input(' option? ','s')
    end
end


switch opt
 case {'r','R'}
  picks = [0,20; 20,25; 25,30; 25,30; 30,40; 40,200];  %% range to pick images from  
  disp('Rows will plot imges with noise ranging')
  disp('from - to')
  disp(picks)
  disp(' ')
 case {'w','W'}
  disp('rows will plot images, starting with the worst')
  worstnoise = 1000;
 case {'t','T'}
  picks = [th-2,th; th-2,th; th-2,th; th,th+2; th,th+2; th,th+2];
  disp('Rows will plot imges with noise ranging')
  disp('from - to')
  disp(picks)
  disp(' ')
otherwise
    disp('Incorrect selection, please select a valid option (r,w,t,R,W,T)')
    return
end

%% decide on comparison
disp(' ')
disp('How do you want to see the data?')
disp('In comparison to the best neighbour = B')
disp('In comparison to the average of the neighbours = A')
dopt = input(' option? ','s')

dopt_ok = 0;
while(dopt_ok < 1)
    if (dopt == 'B' | dopt == 'b' | dopt == 'A' | dopt == 'a')
        dopt_ok = 1;
    else
        disp('Incorrect selection, please select a valid option (b,a,B,A)')
        dopt = input(' option? ','s')
    end
end
  
done = 0;
while(~done)
  %% sample bad slices at random

  disp('Plotting sample slices, please wait')
  
  figure(2), clf
  set(2,'Position',[50,50,400,600])
  
  for pp=1:6  %% plot 6 images on each plot
  
    switch opt
      case {'r','R'}
       [ppi,ppj] = find(noise>picks(pp,1) & noise<picks(pp,2));
       rord = randperm(length(ppi));
       
       slice = ppi(rord(1));
       scan = ppj(rord(1));
     case {'w','W'}
      tmp = noise;
      tmp(noise>=worstnoise) = NaN;
      rmax = max(tmp);
      scan = find(rmax==max(rmax));
      
      rmax = max(tmp');
      slice = find(rmax==max(rmax));
    
      worstnoise = noise(slice,scan);
      
     case {'t','T'}
      [ppi,ppj] = find(noise>picks(pp,1) & noise<picks(pp,2));
      rord = randperm(length(ppi));
      
      slice = ppi(rord(1));
      scan = ppj(rord(1));
    end
      
    %% load bad image
    bim = V(scan);
    [dat1,loc] = spm_read_vols(bim,1);  %% read with zero masking
    slice1 = squeeze(dat1(:,:,slice));
  
    %% load mask
    M = spm_vol(mask);
    msk = squeeze(M.dat(:,:,slice));
    msk(isnan(msk)) = 0;
    
    switch dopt
      case {'b','B'}
       %% pick a neighbour with low noise
       possn = neigh(scan,:);
       possn = possn(isfinite(possn));
       nn = noise(slice,possn);
       best = find(nn ==min(nn));
       nn = possn(best);
    
       [dat2,loc] = spm_read_vols(V(nn),1); 
       slice2 = squeeze(dat2(:,:,slice));

       df = (slice1-slice2).*(msk>0.5);
       ndf = nanmean(abs(df(msk>0.5))); 
       
     case {'a','A'}
      %% get the average of the neighbours
      possn = neigh(scan,:);
      possn = possn(isfinite(possn));
      clear ndf sdf;
      for i=1:length(possn)
        dat2 = spm_read_vols(V(possn(i)),1); 
        slice2 = squeeze(dat2(:,:,slice));
        tdf = (slice1-slice2).*(msk>0.5);
        sdf(:,:,i) = tdf;  %% save
        
        ndf(i) = nanmean(abs(tdf(msk>0.5)));
      end
      df = squeeze(nanmean(shiftdim(sdf,2)));  %% this is just for the picture
      nn = possn(1);
     end
    
    scan_noise = nanmean(ndf);
     
    %% find edges of mask
    bw = edge(imdilate(msk>0.5,strel('disk',1)));
    [mx,my] = find(bw);
    
    figure(2)
    %% plot everything
    subplot(6,3,pp*3-2)
    imagesc(slice1',[0,1000])
    hold on
    plot(mx,my,'w.')
    axis equal
    axis off
    title(['scan: ',num2str(scan),' slice: ',num2str(slice)])
    
    subplot(6,3,pp*3-1)
    imagesc(slice2',[0,1000])
    hold on
    plot(mx,my,'w.')
    axis equal
    axis off
    title(['scan: ',num2str(nn),' slice: ',num2str(slice)])
    
    subplot(6,3,pp*3)
    nmsk = 1-(msk>0.5);
    nmsk(msk<0.5) = NaN;
    imagesc((df+nmsk)',[-250,250])
    axis equal
    axis off
    title(['Noise = ',num2str(scan_noise)])
    
  end
  
  %% give overall title
  ttl = [V(1).fname(end-50:end)];
  ttl = regexprep(ttl,'_','\_');
  tax = axes('Position',[0,0.975,1,1]);
  h= text(0.5,0,ttl);
  set(tax,'xlim',[0,1],'ylim',[0,1])
  set(h,'FontSize',11,'HorizontalAlignment','center','FontWeight','bold')
  axis off

  disp(' ');
  disp('To print this set, type P')
  disp('To see another set, hit RETURN')
  disp('To quit the script, type Q')
  prn = input('  option? ','s')
  if(length(prn)>0)
    if(prn=='q' | prn=='Q')
      disp('Thank you and goodbye')
      return
    elseif(prn=='p' | prn=='P')
      if(~exist('pname'))
        pname = input('Enter filensme for printing with .ps extension ','s')
      end
      printstr = ['print -dpsc2 -painters -append',pname]
      eval(printstr)
      disp(['Figure printed to ',pname])
    end
  end
  
    
end
