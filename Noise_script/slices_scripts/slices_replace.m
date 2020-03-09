%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% replace bad slices with their neighbours

%% written by Antonia Hamilton, Dartmouth College, 2007


clear all, close all

disp('Welcome the the noise replacement script')
disp('This script replaces bad slices in your data')
disp('with their neighbours')

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
%slices_summary(a(ses),'')

%% make sure rejections have been marked for this session
ok = 0;
if(exist('r','var'))
  if(length(r(ses).current)>0)
    ok = 1;
  end
end
if(~ok)
  disp(['Cannot find rejection data for session ',num2str(ses)])
  disp('Cannot replace slices for this session')
  return
end

disp('old data loaded')
reject = r(ses).reject;
current = r(ses).current;
rep1 = r(ses).rep1;
rep2 = r(ses).rep2;
nx = r(ses).nx;
ny = r(ses).ny;
disp(['You checked the worst ',num2str(current),' slices'])

%% plot summary of slices to replace

figure(1), clf
set(1,'Position',[50,50,800,600])

%% set up axes
for i=1:4
  ax(i) = subplot(2,2,i);
end
set(ax(1),'Position',[0.05    0.5811    0.3270    0.3439])
set(ax(2),'Position',[0.48    0.5811    0.5    0.3439])
set(ax(3),'Position',[0.05,0.11, 0.3270    0.3439])
set(ax(4),'Position',[0.48    0.11    0.5    0.3439])


axes(ax(1))
set(gca,'Ylim',[0,8],'Xlim',[0,100],'YDir','reverse')
axis off
h(1)=text(50,1,'Session details:','FontWeight','bold');
h(2)=text(50,2,regexprep([V(1).fname(end-50:end)],'_','\_'));
h(3)=text(50,3,[num2str(nscan),' scans x ',num2str(nslice),' slices = ' ...
                ,num2str(nscan*nslice),' total']);
h(4)=text(50,4,[num2str(sum(isnan(noise(:)))),' slices have movement or no data']);
h(5)=text(50,5,[num2str(sum(reject(:))),' slices rejected (',...
                num2str(100*sum(reject(:))./(nslice*nscan)),'%)']);
h(6)=text(50,6,[num2str(current),' slices examined']);

h(7)=text(50,7,'Go to command line to replace bad slices','FontWeight','bold');
set(h,'HorizontalAlignment','center','FontSize',11)
axis off


axes(ax(2))
imagesc(noise,[0,80])
hold on
for i=1:nslice
  tmp = find(reject(i,:));
  h=plot(tmp,i*ones(length(tmp),1),'kx');
  set(h,'MarkerSize',3)
end
title('noise distribution');

rim = ones(size(reject))*-1;
for i=1:nslice
  rim(i,:) = reject(i,:)*4 - isnan(noise(i,:));
  ind1 = rep1(i,(rep1(i,:)>0));
  ind2 = rep2(i,(rep2(i,:)>0));
  rim(i,unique([ind1,ind2])) = 1;
  rim(i,ind1) = rim(i,ind1)+1;
  rim(i,ind2) = rim(i,ind2)+1;
end

axes(ax(4));
imagesc(rim,[-1,5])
title('red slices will be replaced by green/yellow slices')

axes(ax(3));
[n,b] = hist(noise(:),[0:80]);
bar(b,n,'r');
set(gca,'xlim',[0,80])
title('distribution of slice noise')
[n,b] = hist(noise(~reject),[0:80]);
hold on
bar(b,n,'b')
legend('rejected','kept')

%--------------------------------------------------
disp('Do you want to replace all the bad slices marked?')
disp('Old images will be saved as obold*.img')
disp('New images will be saved as bold*.img')
opt = input('Enter Y to replace, N to to do nothing ','s')

w=1;
if(opt~='y')
  disp('doing nothing')
  return
else
  
  %% first remove mat files
  disp('removing mat files')
  for scan=1:nscan
    mname = V(scan).fname(1:end-3);
    eval(['!rm ',mname,'mat'])
  end

  %% load up mask for plotting
  M = spm_vol(mask);

  %% store filenames
  fname = strvcat(V.fname)
  
  for scan=1:nscan
    if(any(reject(:,scan)))
      disp(['Starting on scan ',num2str(scan)])
      figure(2), clf
      
      %% read in bad data
      bname = deblank(fname(scan,:))
      B = spm_vol(bname);
      bdat = spm_read_vols(B);
      
      %% save data in obold.img
      Vo = B;
      [pth,nam,ext] = fileparts(Vo.fname);
      Vo.fname = [pth,filesep,'o',nam,ext];
      spm_write_vol(Vo,bdat);
   
      %% copy bad image to new image
      ndat = bdat;
      
      %% for each bad slice
      for slice = find(reject(:,scan))'

        %% read in replacement 1
        rname1 = deblank(fname(rep1(slice,scan),:));
        R1 = spm_vol(rname1);
        rdat1 = spm_read_vols(R1);
        slice1 = squeeze(rdat1(:,:,slice));

        %% read in replacement 2
        rname2 = deblank(fname(rep1(slice,scan),:));
        R2 = spm_vol(rname2);
        rdat2 = spm_read_vols(R2);
        slice2 = squeeze(rdat2(:,:,slice));

        %% calculate new slice        
        newslice = (rdat1(:,:,slice)+rdat2(:,:,slice))./2;
        ndat(:,:,slice) = newslice;

        lim = 2000;
        
        subplot(4,4,9), 
        imagesc(squeeze(bdat(:,:,slice))',[0,lim])
        title(['original ',num2str(scan),' ',num2str(slice)])
        axis off
        
        subplot(4,4,10), 
        imagesc(squeeze(rdat1(:,:,slice))',[0,lim])
        title(['rep1 ',num2str(rep1(slice,scan))])
        axis off
        
        subplot(4,4,11), 
        imagesc(squeeze(rdat2(:,:,slice))',[0,lim]);
        title(['rep2 ',num2str(rep2(slice,scan))])
        axis off
        
        subplot(4,4,12), 
        imagesc(squeeze(newslice)',[0,lim]);
        title('new')
        axis off
        pause(0.5);
             
      end
      
      ilim = [max(bdat(:)),max(rdat1(:)),max(rdat2(:)),max(ndat(:))];
      
      subplot(4,4,13)
      hist(bdat(M.dat>0.5),linspace(1,lim,50))
      set(gca,'XLim',[0,lim])
      title(num2str(ilim(1)))
        
      subplot(4,4,14)
      hist(rdat1(M.dat>0.5),linspace(1,lim,50))
      set(gca,'XLim',[0,lim])
      title(num2str(ilim(2)))
              
      subplot(4,4,15)
      hist(rdat2(M.dat>0.5),linspace(1,lim,50))
      set(gca,'XLim',[0,lim])
      title(num2str(ilim(3)))
      
      subplot(4,4,16)
      hist(ndat(M.dat>0.5),linspace(1,lim,50))
      set(gca,'XLim',[0,lim])
      title(num2str(ilim(4)))


      subplot(2,2,1)
      imagesc(squeeze(ndat(:,30,:))');
      title(num2str(scan))
      axis off
      set(gca,'YDir','Normal')

      subplot(2,2,2)
      imagesc(squeeze(ndat(30,:,:))');
      title(num2str(scan))
      axis off
      set(gca,'YDir','Normal')
      pause(1)

      %-Prepare handle for contrast image
      Vout = struct(...
          'fname',  B.fname, ...
          'dim',    B.dim,...
          'mat',    B.mat,...
          'pinfo',  [1,0,0]',...
          'descrip','replaced bold');
      
      %-Write image
      Vout            = spm_create_vol(Vout);
      %Vout.pinfo(1,1) = spm_add(V(scan),Vout);
      Vout = spm_write_vol(Vout,ndat);
      Vout            = spm_close_vol(Vout);
      Vout            = spm_create_vol(Vout,'noopen');
      Vout            = spm_close_vol(Vout);

      w=w+1;
      
      %% write new volume
      disp(['Writing new scan: ',num2str(scan)])
    end
  end
end

disp([num2str(w-1),' new images were written'])
disp('All marked slices have been replaced')
disp('Thank you for using the slices program')

return
%%---------------------------------------------------------------------------
