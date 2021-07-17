%% %%%%% DERIVE CANONICAL PCS
%% %%%%% [avg across stimuli, select good R2, aggregate across subjects, 
%% %%%%%  unitlength normalize, cross-validate to decide 3, derive 3 from all data]

% setup
nsdsetup;

% load the GLM SIMPLE results for all 8 subjects
clear a1 a2;
for p=1:size(sessmatrix,1)
  ppdir = ppdirs{sessmatrix(p,1)};
  a1(p) = load(sprintf('%s/GLMdenoise_nsdSIMPLE.mat',ppdir));
  a2(p) = load(sprintf('%s/GLMdenoise_nsdSIMPLEFIR.mat',ppdir));
end

% use bet to get a quick-and-dirty brain mask
masks = [];  % X x Y x Z x 8 subjects
for p=1:size(sessmatrix,1)
  ppdir = ppdirs{sessmatrix(p,1)};
  file0 = tempname;
  unix_wrapper(sprintf('bet %s %s -v -m -f 0.1', ...
    [ppdir '/preprocessVER1/mean.nii'], ...
    [file0 '.nii']));
  a0 = load_untouch_nii(gunziptemp([file0 '_mask.nii.gz']));
  masks(:,:,:,p) = a0.img;
end

% collect R2
R2 = cat(5,a2.R2);            % X x Y x Z x [FULL/ODD/EVEN] x 8 subjects
meanvol = cat(4,a1.meanvol);  % X x Y x Z x 8
firs = squish(cat(6,a2.firs),3);  % XYZ x [FULL/ODD/EVEN] x 31 TRs x 8 subjects

% define
outdir0 = '~/Dropbox/KKTEMP/hrffigs';
mkdirquiet(outdir0);

% inspect maps
thresh0 = 10;
for p=1:8
  imwrite(uint8(255*makeimagestack(meanvol(:,:,:,p),1)),         gray(256),sprintf('%s/meanvol_%d.png',outdir0,p));
  imwrite(uint8(255*makeimagestack(R2(:,:,:,1,p),[0 40])),       hot(256), sprintf('%s/fullr2_%d.png',outdir0,p));
  imwrite(uint8(255*makeimagestack(masks(:,:,:,p),[0 1])),       gray(256),sprintf('%s/betmask_%d.png',outdir0,p));
  imwrite(uint8(255*makeimagestack(R2(:,:,:,1,p)>thresh0 & masks(:,:,:,p), ...
                                                  [0 1])),       gray(256),sprintf('%s/threshr2_%d.png',outdir0,p));
end

% check numbers of voxels
for p=1:8
  fprintf('%d ',length(find(R2(:,:,:,1,p)>thresh0 & masks(:,:,:,p))));
end
% 17387 15369 9874 10094 11258 5841 7201 2834

% choose a random set of 20000 from each subject [subjects are contributing equally]
num = 20000;
vxselect = [];
for p=1:8
  ix = find(R2(:,:,:,1,p)>thresh0 & masks(:,:,:,p));
  ix = ix(randintrange(1,length(ix),[1 num]));
  vxselect(p,:) = ix;
end

% extract this subset from the firs and normalize based on the ODD data
firsSELECT = [];  % voxels x [FULL/ODD/EVEN] x 31 TRs x 8 subjects
for p=1:8
  vlen = vectorlength(firs(vxselect(p,:),2,:,p),3);
  firsSELECT(:,:,:,p) = bsxfun(@rdivide,firs(vxselect(p,:),:,:,p),vlen);
end

% using ODD timecourses, perform PCA
[u,s,v] = svd(squish(permute(firsSELECT(:,2,:,:),[1 4 2 3]),3),0);

% for increasing numbers of PCs, calculate cross-validated R2
% between the reconstructed timecourse and the left-out timecourse
metricR2 = [];  % 8 subjects x number-of-PCs with R2 values
for p=1:size(v,2)
  for q=1:8
    recon0 = squish(firsSELECT(:,2,:,q),2) * (v(:,1:p)*v(:,1:p)');  % voxels x time
    metricR2(q,p) = calccod(flatten(recon0),flatten(firsSELECT(:,3,:,q)));
  end
end

% inspect performance
figureprep([100 100 500 500]); hold on;
plot(metricR2')
plot(mean(metricR2,1),'k-','LineWidth',3);
xlabel('Number of PCs');
ylabel('Cross-validated variance explained');
figurewrite('crossvalidation',[],-1,outdir0);

% decide that we want 3!

% now, using all of the timecourses, perform PCA
[u,s,v] = svd(squish(permute(firsSELECT(:,1,:,:),[1 4 2 3]),3),0);

% manually flip?
v(:,1) = -v(:,1);
%v(:,2) = -v(:,2);
v(:,3) = -v(:,3);

% inspect canonical set of PCs
temp = flatten(diag(s(1:3,1:3))).^2;
temp = temp / sum(temp);

figureprep([100 100 1000 500],1);

subplot(1,2,1); hold on;
plot(0:30,temp(1)*v(:,1),'ro-');
plot(0:30,temp(2)*v(:,2),'go-');
plot(0:30,temp(3)*v(:,3),'bo-');
straightline(0,'h','k-');
xlabel('Time from stimulus onset (s)');
ylabel('Signal');

subplot(1,2,2); hold on;
plot(0:30,v(:,1),'ro-');
plot(0:30,v(:,2),'go-');
plot(0:30,v(:,3),'bo-');
straightline(0,'h','k-');
xlabel('Time from stimulus onset (s)');
ylabel('Signal');

figurewrite('pcs',[],-1,outdir0);

% save results
save('~/nsd/ppdata/hrfbasis.mat','vxselect','s','v','metricR2');

%%%%%% make a 4/3 version

% interpolate
vnew = [];
for p=1:size(v,2)
  vnew(:,p) = interp1(0:30,v(:,p)',0:(4/3):30,'pchip');
end

% orthogonalize
for p=2:size(vnew,2)
  vnew(:,p) = projectionmatrix(vnew(:,1:p-1))*vnew(:,p);
end

% truncate
vnew = vnew(:,1:23);

% make unit length
vnew = unitlength(vnew,1);

% set
v = vnew;

% visualize
figureprep([100 100 500 500]); hold on;
plot(0:4/3:30,v(:,1),'ro-');
plot(0:4/3:30,v(:,2),'go-');
plot(0:4/3:30,v(:,3),'bo-');
straightline(0,'h','k-');
xlabel('Time from stimulus onset (s)');
ylabel('Signal');
figurewrite('pcs_3pertrial',[],-1,outdir0);

% save results
save('~/nsd/ppdata/hrfbasis_3pertrial.mat','v');
