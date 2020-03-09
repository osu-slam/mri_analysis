function slices_summary(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this function summarises the slice noise for each session
% function slice_noise_summary(a,th,wth,prn)
% inputs are - 
%    a - structure from slices_analyse
%    prn = name to print 
%% if function is called without any inputs, user is prompted for info
%
% written by Antonia Hamilton, Dartmouth College, 2007
  
  disp('Welcome the Slices Summary script')
  disp(' ')
  
  spm_defaults  %% assume spm2
  
  try
    load slices_def
  catch
    disp('Could not load slices_def.mat')
    disp('Please run the slices_defaults program')
    disp('and make sure you are in the right directory')
    return
  end
  
  if(nargin==2)
    
    a = varargin{1};
    prn = varargin{2};
    
    ses = 1;
    
    P = a(ses).P;
    noise = a(ses).noise;
    mask = a(ses).mask;
    ra = a(ses).ra;
    neigh = a(ses).neigh;

    %% all is good
  elseif(nargin==1)
    
    a = varargin{1};
    prn = '';
    
    ses = 1;
    
    P = a(ses).P;
    noise = a(ses).noise;
    mask = a(ses).mask;
    ra = a(ses).ra;
    neigh = a(ses).neigh;

  elseif(nargin==0)
    
    disp('Pick the saved noise file to analyse')
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
    
     disp('Enter a filename for printing the figure with .ps at the end ')
    prn = input('  or press ENTER for no printing ','s')
    if(length(prn)>0)
      disp(['Figure will be appended to ',prn])
    end
    disp('Press space to start')
  else
    disp('wrong number of inputs')
    return
  end
  
  V = spm_vol(P);
  %%%%%%%%%%%%%%%%%%%%%%
  %% now do the plotting
  
  th = sdef.th;
  wth = sdef.wth;
  
  %% plot summary of everything for this session
figure(1), clf

subplot(3,1,1)
imagesc(noise,[0,80])
ttl = [V(1).fname(end-50:end)];
ttl = regexprep(ttl,'_','\_');
title(ttl);

subplot(3,1,2)
[n,b] = hist(noise(:),[0:80]);
bar(b,n);
set(gca,'xlim',[0,80])
title('distribution of slice noise')
hold on
plot([th,th],[0,max(n)],'r-')
h=text(th,max(n),[num2str(100*sum(noise(:)>th)./prod(size(noise)),3),...
                  '% of slices are over ',num2str(th)]);
set(h,'Color',[1 0 0])

plot([wth,wth],[0,max(n)/2],'g-')
h=text(wth,max(n)/2,[num2str(100*sum(noise(:)>wth)./prod(size(noise)),3),...
                  '% of slices are over ',num2str(wth)]);
set(h,'Color',[0 1 0])

subplot(6,1,5)
plot(ra(:,1:3))
title('translation (mm)')
set(gca,'XTick',[])

subplot(6,1,6)
plot(ra(:,4:6)*180/pi)
title('rotation (deg)')

%% set to full page on Letter paper
set(1,'PaperPosition',[0.25,0.25,8,10.5])

pause(0.1)
if(length(prn)>0)
  prnstr = ['print -dpsc2 -painters -append ',prn]
  %prnstr = ['print -dpdf -append ',prn] prints to pdf instead
  eval(prnstr)
end
