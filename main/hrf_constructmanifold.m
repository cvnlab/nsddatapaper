%% %%%%% From the hrfbasis fit (to all runs), visualize and record coordinates on unit sphere
%% [do a GLM fit in the special hrfbasis(3PCs) and then look at distribution of the loadings for good voxels]

nsdsetup;

ns = [];
allpoints = {};
for p=1:size(sessmatrix,1)
  for q=1:size(sessmatrix,2)
    if sessmatrix(p,q)~=0

      % use bet to get a quick-and-dirty brain mask
      ppdir = ppdirs{sessmatrix(p,q)};
      file0 = tempname;
      unix_wrapper(sprintf('bet %s %s -v -m -f 0.1', ...
        [ppdir '/preprocessVER1SECONDSTAGE/mean.nii'], ...
        [file0 '.nii']));
      a0 = load_untouch_nii(gunziptemp([file0 '_mask.nii.gz']));
      mask = a0.img;

      % load    
      load(sprintf('%s/GLMdenoise_nsdBASIConoffhrfbasis.mat',glmdirs{sessmatrix(p,q)}),'R2','meanvol','modelmd');

      % we only want the full fit
      ii = 1;

      % get data
      r0 = squish(R2(:,:,:,ii),3);        % XYZ x 1
      m0 = squish(modelmd(:,:,:,:,ii),3); % XYZ x 3
  
      % flip each voxel such that loading on first component is positive
      m0 = bsxfun(@times,m0,sign(m0(:,1)));
      
      % map to the unit sphere
      m0 = unitlength(m0,2);
  
      % voxel selection (reasonable R2 and in the brain)
      thresh0 = 10;
      good = r0>thresh0 & mask(:)==1;
      
      % get points
      XXX0 = m0(good,2);
      YYY0 = m0(good,3);
      ZZZ0 = m0(good,1);

      % calc
      bins = -1.5:.02:1.5;
      [n,xx,yy] = hist2d(XXX0,YYY0,bins,bins);

      % record
      if q==1
        allpoints{p} = [XXX0 YYY0 ZZZ0];
        ns(:,:,p) = n/sum(n(:));
      end

      % visualize
      figureprep; hold on;
      imagesc(xx(1,:),yy(:,1),n);
      colormap(hot);
      axis equal tight; 
%      straightline(0,'h','w-');
%      straightline(0,'v','w-');
      xlabel('loading on PC2');
      ylabel('loading on PC3');
      drawellipse(0,0,0,1,1,[],[],'w-');
      figurewrite(sprintf('hist2d_p%02d_q%02d',p,q),[],-1,'~/Dropbox/KKTEMP/hist2d');
    
    end
  end
end

% visualization of aggregate across all nsd01 sessions
fig = figureprep([],1); hold on;
imagesc(xx(1,:),yy(:,1),mean(ns,3));
colormap(hot);
axis equal tight; 
%straightline(0,'h','w-');
%straightline(0,'v','w-');
xlabel('loading on PC2');
ylabel('loading on PC3');
drawellipse(0,0,0,1,1,[],[],'w-');
figurewrite(sprintf('hist2dall'),[],[],'~/Dropbox/KKTEMP/hist2d',1);

%% %%%%% Parametrize manifold via clicks and great-circle interpolation

% define line
[xypoint,xyline,pixelmask] = roiline(fig);
xypoint = [
         0.686363636363637        -0.255932203389831
         0.560724191063175        -0.153559322033898
         0.407164869029276       -0.0604930662557783
         0.221032357473036        0.0279198767334361
        0.0814329738058555         0.065146379044684
        -0.072126348228043        0.0837596302003081
        -0.221032357473035        0.0791063174114019
        -0.369938366718028        0.0418798151001538
         -0.51884437596302       -0.0279198767334363
        -0.644483821263482       -0.0930662557781203
        -0.765469953775038        -0.195439137134052
        -0.863189522342065        -0.311771956856703
];

% make nice unit-sphere points (1st column is inward, 2nd and 3rd columns are x- and y-)
pts = [sqrt(1-sum(xypoint.^2,2)) xypoint];  % clicked-points x 3

% uniformly distribute 200 particles across the surface of a unit sphere.
% let's take V(:,3) to be the inward dimension
[V,Tri,~,Ue] = ParticleSampleSphere('N',100*2);
fv=struct('faces',Tri,'vertices',V); 
fv_new=SubdivideSphericalMesh(fv,2);   % subdivide to become denser
V = fv_new.vertices;
V = V(V(:,3)>0,:);  % T x 3

% how dense should we make it??
angularsep = 6;

% construct path
thepath = constructpathonsphere(pts,angularsep);  % new-points x 3

% do a more realistic rendering
  % do the histogram
counts = [];  % subjects x locations with relative proportions
for p=1:length(allpoints)
  [mm,ii] = min(acos(allpoints{p}*V'),[],2);
  counts0 = [];
  for q=1:size(V,1)
    counts0(q) = sum(ii==q);
  end
  counts(p,:) = counts0 / sum(counts0);
end
figureprep([],1); hold on;
h = scatter3(V(:,1),V(:,2),V(:,3),'ro','filled');
set(h,'CData',mean(counts,1));
set(h,'SizeData',13^2);
colormap(hot);
cax = caxis;
caxis([0 cax(2)]);
sc = 1.05;  % expand a little for visibility
h3 = scatter3(sc*thepath(:,2),sc*thepath(:,3),sc*thepath(:,1),'co','filled');
h2 = scatter3(sc*pts(:,2),    sc*pts(:,3),    sc*pts(:,1),'bo');
axis equal;
xlabel('PC2');
ylabel('PC3');
zlabel('PC1');
title(sprintf('there are %d points',size(thepath,1)));
figurewrite('manifold',[],-1,'~/Dropbox/KKTEMP');

% visualize these HRFs. what are they like?
a1 = load('~/nsd/ppdata/hrfbasis.mat');
tc = [];  % time (0-30) x HRFs
for p=1:size(thepath,1)
  tc(:,p) = a1.v(:,1:3)*thepath(p,:)';  % unit-length basis functions scaled by the coordinates on the sphere
end
figureprep([],1); hold on;
cmap0 = cmapsinglecolor(size(tc,2),[0 0 0],cmaptab10([],10));
hrf1s = [];    % 31 time points x 24 hrfs
hrf43s = [];   % 23 time points x 24 hrfs
for p=1:size(tc,2)
  tc0 = interp1(0:30,tc(:,p),0:.1:30,'pchip');
  plot(0:.1:30,tc0/max(tc0),'r-','Color',cmap0(p,:));  % upsample and ensure max is 1
  hrf1s(:,p) =  interp1(0:.1:30,tc0/max(tc0),0:30,'pchip');
  hrf43s(:,p) = interp1(0:.1:30,tc0/max(tc0),0:4/3:30,'pchip');
end
hrf = getcanonicalhrf(3,1);
plot(0:30,hrf(1:31),'g-');
  % how similar are nearby hrfs?
cmatrix = calcconfusionmatrix(hrf1s,[],2);
% figure;imagesc(cmatrix)
% colormap(gray)
% caxis([0 1])
% colorbar;
% axis image tight;
temp = reshape(1:numel(cmatrix),size(cmatrix));
avgr2 = mean(cmatrix(diag(temp,1))'.^2 * 100);
title(sprintf('avg adjacent r2 = %.3f',avgr2));
figurewrite('hrfs',[],-1,'~/Dropbox/KKTEMP');

% save master hrfs into official files
hrfs = hrf1s;
save('~/nsd/ppdata/hrfmanifold.mat','xypoint','pts','angularsep','thepath','hrfs');  % hrfs_func1mm.mat
hrfs = hrf43s;
save('~/nsd/ppdata/hrfmanifold_3pertrial.mat','hrfs');  % hrfs_func1pt8mm.mat
