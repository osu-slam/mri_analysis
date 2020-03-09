%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% inspect bad slices and reject the ones you don't like
%% you must run slices_analyse first so you have some data

%% written by Antonia Hamilton, Dartmouth College, 2007


clear all, close all

disp('Welcome the the noise rejection script')
disp('This script lets you inspect your session for noise')
disp('and decide which slices to reject and replace')
disp('The actual replacement happens later')
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
slices_summary(a(ses),'')

neednew = 1;
%% see if old data was saved
if(exist('r','var'))
  if(length(r)>=ses)
    if(length(r(ses).current)>0)
      disp('old data loaded')
      reject = r(ses).reject;
      current = r(ses).current;
      rep1 = r(ses).rep1;
      rep2 = r(ses).rep2;
      nx = r(ses).nx;
      ny = r(ses).ny;
      current = input(['Enter starting slice (you got to ',num2str(current),...
                       ' last time) '])
      neednew = 0;
    end
  end
end

if(neednew)
  %% find new orders etc
  
  %% rank slices by noise
  cnoise = reshape(noise,nslice*nscan,1);
  cnoise(isnan(cnoise)) = -20;
  [y,ii] = sort(cnoise);
  %% flip order of ii
  ii = flipud(ii);
  nrank(ii) = 1:length(ii);
  nrank(isnan(noise)) = NaN;
  
  %% get slices & scans in order
  [nx,ny] = ind2sub(size(noise),ii);
  
  current = 1;
  
  reject = zeros(size(noise));
  rep1 = zeros(size(noise));
  rep2 = zeros(size(noise));
end

rej = {'keep','reject'};

%%%%%%%%%%%
%% draw permenant bits of plot
figure(2), clf
set(2,'Position',[50,50,800,600])

sax = subplot(3,2,1); 
zax = subplot(3,2,2);

%% set up plots for slices
pax(1) = subplot(3,4,5);
pax(2) = subplot(3,4,6);
pax(3) = subplot(3,4,7);
pax(4) = subplot(3,4,8);
pax(5) = subplot(3,4,11);
pax(6) = subplot(3,4,12);

bw = 100;
xpos = linspace(50,750-bw,6);

%% set up buttons
b(1) = uicontrol(2,'Style','PushButton',...
                 'String','Previous',...
                 'Tag','Previous',...
                 'Userdata',0,...
                 'CallBack',['set(gcbo,''UserData'',1)'],...
                 'Position',[xpos(1),20,90,30]);

b(2) = uicontrol(2,'Style','PushButton',...
                 'String','Next',...
                 'Tag','Next',...
                 'Userdata',0,...
                 'CallBack',['set(gcbo,''UserData'',1)'],...
                 'Position',[xpos(2),20,90,30]);

b(3) = uicontrol(2,'Style','PushButton',...
                 'String','Keep',...
                 'Tag','Keep',...
                 'Userdata',0,...
                 'CallBack',['set(gcbo,''UserData'',1)'],...
                 'Position',[xpos(3),20,90,30]);

b(4) = uicontrol(2,'Style','PushButton',...
                 'String','Reject',...
                 'Tag','Reject',...
                 'Userdata',0,...
                 'CallBack',['set(gcbo,''UserData'',1)'],...
                 'Position',[xpos(4),20,90,30]);

b(5) = uicontrol(2,'Style','PushButton',...
                 'String','Save',...
                 'Tag','save',...
                 'Userdata',0,...
                 'CallBack',['set(gcbo,''UserData'',1)'],...
                 'Position',[xpos(5),20,90,30]);

b(6) = uicontrol(2,'Style','PushButton',...
                 'String','Quit',...
                 'Tag','Quit',...
                 'Userdata',0,...
                 'CallBack',['set(gcbo,''UserData'',1)'],...
                 'Position',[xpos(6),20,90,30]);


done = 0;
while(~done)

  %% set all buttons to zero
  set(b,'Userdata',0);
  
  %% determine current scans etc
  scan = ny(current);
  slice = nx(current);
  
  %% load mask
  M = spm_vol(mask);
  msk = squeeze(M.dat(:,:,slice));
  msk(isnan(msk)) = 0;
  nvox = sum(msk(:)>0.5);
  
  %% load up the slice
  bim = V(scan);
  [dat1,loc] = spm_read_vols(bim,1);  %% read with zero masking
  slice1 = squeeze(dat1(:,:,slice));
  mx(1) = max(dat1(:));
  md(1) = nanmean(dat1(msk(:)>0.5));

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

  %% find and load up proposed replacements
  %% best possible replacements are adjacent slices
  ppr = [scan-1,scan+1];
  wrn = '';
  if(any(ppr<1) | any(ppr>nscan) )
    wrn = ['Warning, hard to find good replacments'];
  elseif(any(noise(slice,ppr)>th))
    wrn = ['Warning, hard to find good replacments'];
  end
  if(length(wrn)>1) 
    
    possn = [-5:-1,1:5]+scan;
    possn = possn(possn>0);
    possn = possn(possn<=nscan);
    score = ((possn-scan).^2)*5+noise(slice,possn);
    
    [tmp,i] = sort(score);
    best = i(1:2);
    
    
    ppr = possn(best);
  end
    
  for i=1:2
    dat2 = spm_read_vols(V(ppr(i)),1); 
    prslice(:,:,i) = squeeze(dat2(:,:,slice));
    mx(i+1) = max(dat2(:));
    md(i+1) = nanmean(dat2(msk(:)>0.5));
  end
  
  %----------- start the drawing -------------------------------------  
  
  figure(2),
  %% set up slice info axis
  axes(sax), cla
  set(sax,'xlim',[0,100],'ylim',[0,6],'ydir','reverse')
  h=text(50,1,['Slice no: ',num2str(slice),'  Scan no: ',num2str(scan)]);
  set(h,'HorizontalAlignment','center','FontSize',12,'FontWeight','bold')
  h=text(50,2,regexprep([V(scan).fname(end-50:end)],'_','\_'));
  set(h,'HorizontalAlignment','center')
  h=text(50,3,['noise = ',num2str(noise(slice,scan),3),...
               '   rank = ',num2str(current),...
               '   nvox = ',num2str(nvox)]);
  set(h,'HorizontalAlignment','center','FontSize',11)
  h=text(50,4,['Image has ',num2str(sum(isfinite(neigh(scan,:)))),' good ' ...
               'neighbours']);
  set(h,'HorizontalAlignment','center','FontSize',11)
  hstatus=text(50,5,['Status is ',rej{reject(slice,scan)+1}]);  
  set(hstatus,'HorizontalAlignment','center','FontSize',12,'FontWeight', ...
              'bold');
  hwrn=text(50,6,wrn,'HorizontalAlignment','center','FontSize',11,'FontWeight', ...
              'bold','color',[1 0 0]);
  axis off
  
  %% set up session info
  axes(zax), cla 
  set(zax,'xlim',[0,100],'ylim',[0,7],'ydir','reverse') 
  h(1)=text(50,1,'This session has:'); 
  h(2)=text(50,2,[num2str(nscan),' scans x ',num2str(nslice),' slices = ' ...
                  ,num2str(nscan*nslice),' total']);
  h(3)=text(50,3,[num2str(sum(isnan(noise(:)))),' slices have movement or no data']);
  h(4)=text(50,4,[num2str(sum(reject(:))),' slices rejected so far (',...
                  num2str(100*sum(reject(:))./(nslice*nscan)),'%)']); 
  h(5)=text(50,5,[num2str(current),' slices examined so far']); 
  
  h(6)=text(50,6,['mean is: ',num2str(md,3)]);
  h(7)=text(50,7,['max is: ',num2str(mx)]);
  set(h,'HorizontalAlignment','center','FontSize',11)
  axis off
  
  %% plot the bad slice 
  axes(pax(1))
  imagesc(slice1',[0,max(mx)]) 
  title('bad slice')
  axis equal 
  axis off 
  
  %% plot the difference from the neighbours
  axes(pax(2)) 
  nmsk = 1-(msk>0.5);
  nmsk(msk<0.5) = NaN; 
  imagesc((df+nmsk)',[-250,250])
  axis equal
  axis off
  title('noise')
  
  %% plot movement and noise for this slice and its neighbours
  xx = -7:7;
  tscan = xx+scan;
  tscan = tscan(tscan>0);
  tscan = tscan(tscan<=nscan);
  mvt = ra(tscan,1:3);
  nnn = noise(slice,tscan);
  
  axes(pax(3)), cla
  plot(tscan,mvt-repmat(nanmean(mvt),length(mvt),1),'-')
  hold on
  plot([scan,scan],[-1,1],'k-')
  title('Mvt params')
  set(gca,'XLim',[min(tscan),max(tscan)])
  
  axes(pax(5)), cla
  plot(tscan,nnn,'r-o')
  hold on
  plot([scan,scan],[10,nnn(8)],'k-')
  plot(tscan,ones(length(tscan),1)*th,'g--')
  plot(scan,noise(slice,scan),'k*')
  plot(ppr(1),noise(slice,ppr(1)),'b*')
  plot(ppr(2),noise(slice,ppr(2)),'b*')
  title('noise')
  set(gca,'XLim',[min(tscan),max(tscan)])
  
  %% plot noise for all scans
  subplot(3,2,5), cla
  set(gca,'Position',[0.05    0.12    0.45    0.22])
  imagesc(noise,[0,80])
  hold on
  h=plot(scan,slice,'ko');
  set(h,'MarkerSize',10,'LineWidth',2)
  for i=1:nslice
    tmp = find(reject(i,:));
    h=plot(tmp,i*ones(length(tmp),1),'kx');
    set(h,'MarkerSize',3)
  end
  title('overall noise')
  axis off
  axis tight
  
  %% plot the possible replacements
  axes(pax(4))
  imagesc(squeeze(prslice(:,:,1))',[0,max(mx)])
  title(['Replacement 1 = ',num2str(ppr(1))])
  axis equal
  axis off
  
  axes(pax(6))
  imagesc(squeeze(prslice(:,:,2))',[0,max(mx)])
  title(['Replacement 2 = ',num2str(ppr(2))])
  axis equal
  axis off
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% read in buttons and decide what to do
  tmp = get(b,'UserData');
  hits = cat(2,tmp{:});
  
  while(all(hits==0))  %% stay in this little loop until a button is hit
    drawnow
    tmp = get(b,'UserData');
    hits = cat(2,tmp{:});
  end
  
  if(all(hits==0))
    %% do nothing
    pause(0.1)
  else
    switch(find(hits))
      
     case 1  %% prev
      current = current-1;
      
     case 2  %% next
      current = current+1;
      
     case 3 %% button is keep
      reject(slice,scan) = 0;
      rep1(slice,scan) = 0;
      rep2(slice,scan) = 0;
      axes(sax);
      set(hstatus,'String',['Status is ',rej{reject(slice,scan)+1}],'Color',[0,1,0]);
      pause(0.5);
      current = current+1;
      
     case 4 %% button is reject
      reject(slice,scan) = 1;
      rep1(slice,scan) = ppr(1);
      rep2(slice,scan) = ppr(2);
      axes(sax);
      set(hstatus,'String',['Status is ',rej{reject(slice,scan)+1}],'Color',[1,0,0]);
      pause(0.5);
      current = current+1;
      
     case 5  %% button is save
      disp('you hit save')
      axes(sax);
      set(hstatus,'String',['Saving'],'Color',[0,0,1])
      drawnow
      
      r(ses).reject = reject;
      r(ses).current = current;
      r(ses).rep1 = rep1;
      r(ses).rep2 = rep2;
      r(ses).nx = nx;
      r(ses).ny = ny;
      
      save(inname,'a','r')
      disp(['Reject markings saved to file: ',inname])
      set(hstatus,'String',['Saved to ',inname],'Color',[1,0,0])
      drawnow
      pause(1);
      
     case 6  %% button is quit
      disp('you hit quit')
      set(hstatus,'String',['Go to the command line to save'],'Color',[0,0,1])
      drawnow
      pause(0.5);
      done = 1;
      
    end
  end
end

%figure(3)
%imagesc(reject)
%title('Red = rejected slices')

disp('Do you want to save the data just to be sure')
opt = input('option? (y/n) ','s')
if(opt=='y')
     r(ses).reject = reject;
      r(ses).current = current;
      r(ses).rep1 = rep1;
      r(ses).rep2 = rep2;
      r(ses).nx = nx;
      r(ses).ny = ny;
      
      save(inname,'a','r')
      disp(['Reject markings saved to file: ',inname])
else
  disp('data not saved')
end



return      
  