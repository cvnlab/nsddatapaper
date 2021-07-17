% define
subjix = 1;     % main subject of interest
roimainix = 1;  % main ROI of interest
sess = 10;   % session of interest
roistoload = {{'rh.floc-faces.nii.gz' 2} {'rh.floc-places.nii.gz' 2}};  % RH FFA-1 and RH PPA
betadirs = {'betas_assumehrf' 'betas_fithrf' 'betas_fithrf_GLMdenoise_RR'};
ncfun = @(x)100*(x.^2./(x.^2+1/3));
nsessions = [40 40 32 30 40 32 40 30];

% load exp design stuff
a1 = load('~/nsd/nsddata/experiments/nsd/nsd_expdesign.mat');

% load from all subjects...
a2 = {}; data = {}; noiseceiling = {};  % each of these is a cell vector of 8 subjects x 2 ROIs
for sss=1:8
  for rrr=1:length(roistoload)

    % load desired ROI
    a2{sss,rrr} = load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj%02d/func1pt8mm/roi/%s',sss,roistoload{rrr}{1}));
    [d1,d2,d3,ii] = computebrickandindices(a2{sss,rrr}.img==roistoload{rrr}{2});

    % load entire dataset (for all 3 versions)
    data0 = [];  % voxels x 750 x 40 x betatype
    for bb=1:3, bb
      for ss=1:nsessions(sss), ss
        file0 = sprintf('~/nsd/nsddata_betas/ppdata/subj%02d/func1pt8mm/%s/betas_session%02d.hdf5',sss,betadirs{bb},ss);
        temp = squish(double(h5read(file0,'/betas',[d1(1) d2(1) d3(1) 1],[range(d1)+1 range(d2)+1 range(d3)+1 750]))/300,3);
        data0(:,:,ss,bb) = temp(ii,:);
        if bb==1 && ss==1
          data0(end,750,nsessions(sss),3) = NaN;
        end
      end
    end

    % pick a random order for the voxels (note that this is deterministic)
    reorderix = picksubset(1:size(data0,1),size(data0,1));
    data0 = data0(reorderix,:,:,:);  % reorder the voxels
    data{sss,rrr} = data0;  % record

    % load ncsnr and compute noise ceiling values
    ncsnr = [];  % voxels x betatypes
    for bb=1:3
      file0 = sprintf('~/nsd/nsddata_betas/ppdata/subj%02d/func1pt8mm/%s/ncsnr.nii.gz',sss,betadirs{bb});
      temp = load_untouch_nii(file0);
      temp = temp.img(a2{sss,rrr}.img==roistoload{rrr}{2});
      ncsnr(:,bb) = temp(reorderix);
    end
    noiseceiling{sss,rrr} = ncfun(ncsnr);
  
  end

end

% carpet plots
rngs = {1:750*40 (sess-1)*750+(1:750)};
for bb=1:3
  for rr=1:length(rngs)
    for wantz=0:1
      temp = data{subjix,roimainix}(:,:,:,bb);    % voxels x 750 x 40
      if wantz
        temp = calczscore(temp,2);
      end
      temp = temp(:,rngs{rr});  % voxels x trials
      if wantz==0
        datarng = [-10 10];
      else
        datarng = [-3 3];
      end
      imwrite(uint8(255*makeimagestack(temp,datarng,0)),cmapsign4(256),sprintf('~/Dropbox/KKTEMP/b%d_rng%d_wantz%d.png',bb,rr,wantz));
    end
  end
end

% calc triplets
magicix = [];  % 3 x 10000, elements are 1-indexes indicating which trial numbers were used for the 10,000 distinct images
for q=1:10000
  magicix(:,q) = find(a1.masterordering==q);
end

% calc vmetric values
temp = calczscore(data{subjix,roimainix},2);
vmetric = [];  % voxels x 3
for bb=1:3
  temp2 = temp(:,:,:,bb);    % voxels x 750 x 40
  temp2 = reshape(temp2(:,flatten(magicix)),[],3,10000);
  vmetric(:,bb) = sqrt(mean(std(temp2,[],2).^2,3));
end

% CHOOSE A VOXEL
vxix = 99;    %  77:151, vxix
fprintf('after reordering, looking at voxel %d out of %d.\n',vxix,size(vmetric,1));
vxix = 64;  % this is the old 99

% report some information
subscript(find(a2{subjix,roimainix}.img==roistoload{1}{1}),reorderix(vxix))
vmetric(vxix,:)
snr0 = translatevmetric(vmetric(vxix,:))
ncfun(translatevmetric(vmetric(vxix,:)))

% after reordering, looking at voxel 99 out of 160.
% 
% ans =
% 
%       271898
% 
% 
% ans =
% 
%          0.833974608757999          0.77663357443613         0.729783169818409
% 
% 
% snr0 =
% 
%          0.661654220658782         0.811132342078095         0.936824387269245
% 
% 
% ans =
% 
%           56.7728124678659          66.3730702693165          72.4739312819372

% visualize for each beta version some stuff
%cmap0 = hot(7);
%cmap0 = cmap0(1:3,:);
cmap0 = [0 0 0; 1 0 0; 0 0 1];
nimages = 100;
dataZ = calczscore(data{subjix,roimainix},2);
figureprep([100 100 350 350],1); hold on;
for bb=1:3
  subplot(3,1,bb); hold on;
  temp = squish(dataZ(vxix,:,:,bb),2);  % 750 x 40
  datatemp = reshape(temp(magicix(:,1:nimages)),3,[]);  % 3 x 100   [note that the image number is just 1 through 10,000]
  for im=1:nimages
    set(straightline(mean(datatemp(:,im)),'h','k-',[im-.4 im+.4]),'Color',[.5 .5 .5]);
    h = scatter(im*ones(1,3),datatemp(:,im)',49,cmap0,'o','filled');
  end
  effectivesignalsd = snr0(bb)/sqrt(1+snr0(bb)^2);
  effectivenoisesd  = 1/sqrt(1+snr0(bb)^2);
  yy = -4:.1:4;
  plot(evalgaussian1d([0 effectivesignalsd 1 11],yy),yy,'r-');
  plot(evalgaussian1d([0 effectivenoisesd 1 11],yy),yy,'m-');
  for p=1:3
    subplot(3,1,p);
    axis([10.5 25.5 -4 4]);
    set(gca,'XTick',1:1:100);
  end 
  xlabel('Image');
  ylabel('Response (z-score units)');
end
figurewrite('examplevoxel',[],-2,'~/Dropbox/KKTEMP');

%%% MORE VISUALIZATIONS

% calc some stuff
temp = a1.masterordering(1:750*30);
shared515 = [];  % 1 x 515 vector of 1-indices. these indices are between 1-1000.
for q=1:1000
  if sum(temp==q)==3
    shared515 = [shared515 q];
  end
end
magic515ix = [];  % 3 x 515, elements are 1-indexes indicating which trial numbers were used for the 515 shared images
for q=1:length(shared515)
  magic515ix(:,q) = find(a1.masterordering==shared515(q));  
end

% prepare some data matrices
whb = 2;  % ONLY b2 !
prepmatrix = [];  % 8*2 x 515 x 3 reps
for sss=1:8
  for rrr=1:2
    temp2 = data{sss,rrr}(:,:,:,whb);  % voxels x 750 x sessions
    temp2 = reshape(temp2(:,flatten(magic515ix)),[],3,515);  % voxels x 3 reps x 515 images
    temp2 = squish(mean(temp2,1),2);  % 3 reps x 515 images
    prepmatrix((rrr-1)*8+sss,:,:) = reshape(temp2',1,515,3);
  end
end

% within and across-subject variability, for two ROIs
collabels = {'Within-subject variability' 'Across-subject variability'};
figureprep([100 100 650 275]); hold on;
for rr=1:2
  for typ=1:2  % within, across
    subplot(2,2,(rr-1)*2+typ); hold on;
    if typ==1
      grandmn = [];
      cmap0 = [.3 .3 .3; 1 0 0; 0 0 1];
      for rep=1:3
        temp0 = prepmatrix((rr-1)*8+1,:,rep);  % only subject 1
        plot(1:515,temp0,'-','Color',cmap0(rep,:)); 
        grandmn(rep,:) = temp0;
      end
      plot(1:515,mean(grandmn,1),'k-','LineWidth',2);
    else
      grandmn = [];
      cmap0 = jet(8);
      cmap0(3:5,:) = (cmap0(3:5,:)*5 + [0 0 0])/6;
      for sss=1:8
        temp0 = mean(prepmatrix((rr-1)*8+sss,:,:),3);
        plot(1:515,temp0,'-','Color',cmap0(sss,:));
        grandmn(sss,:) = temp0;
      end
      plot(1:515,mean(grandmn,1),'k-','LineWidth',2);
    end
    axis([0 50 -1 6]);
    straightline(0,'h','k-');
    if rr==1
      title(collabels{typ});
    end
    if typ==1
      ylabel('BOLD (% change)');
    end
    if rr==2
      xlabel('Shared515 image number');
    end
  end
end
figurewrite('withinacross',[],-2,'~/Dropbox/KKTEMP');

% confusion matrix
figureprep([100 100 275 275]); hold on;
temp = squish(permute(reshape(prepmatrix,[8 2 515 3]),[4 1 2 3]),3);
imagesc(calcconfusionmatrix(temp',[],2));
caxis([-1 1]);
colormap(cmapsign4);
colorbar;
axis image;
set(gca,'YDir','reverse');
straightline(.5+(0:3:48),'h','c-');
straightline(.5+(0:3:48),'v','c-');
straightline(.5+(0:3*8:48),'h','r-');
straightline(.5+(0:3*8:48),'v','r-');
axis off;
figurewrite('confusionmatrix',[],-2,'~/Dropbox/KKTEMP');

%%%%%%%%%%%%%%%%%%%%%%% JUNK:

%  set(errorbar2(1:nimages,mean(datatemp,1),std(datatemp,[],1),'v','k-'),'LineWidth',2);
%  errorbar3(1:nimages,mean(datatemp,1),repmat(vmetric(vxix,bb),1,nimages),'v',[1 .5 .5]);
%  errorbar3(1:nimages,mean(datatemp,1),std(datatemp,[],1),'v',[1 .5 .5]);
%   plot(1:nimages,mean(datatemp,1),'k-','LineWidth',2);
%   bar(1:nimages,mean(datatemp,1));
%  plot([.5 .5],vmetric(vxix,bb)*[-1 1],'b-','LineWidth',2);
%     thedata = temp(magicix(:,im));
%     resid = abs(thedata-mean(thedata));
%     for rep=1:3
%       h = scatter(im,thedata(rep),'ko','filled');
% %      set(h,'CData',cmaplookup(resid(rep),0,.5,0,copper(256)));
%     end
%   pause;
% end
