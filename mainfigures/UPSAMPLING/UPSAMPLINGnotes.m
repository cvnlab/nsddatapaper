%% %%%%%%%%%%%%% Create views in ITK-SNAP

% subj01.

% anterior, right, inferior is the corner that matches 1mm and 1.8mm
%
% from center of first voxel (1,1,1) in 1mm, we achieve correspondence every +9 1mm slices or +5 1.8mm slices.

(.,24,.)
186 - 162 = 24
162 is 9*18.
104 - 162/1.8 = 14

(.,.,73) for 1mm is good
1+9*8 = 73
1+72/1.8 = 41

(82,.,.) for 1mm is good
145-9*7=82
81-63/1.8=46

1mm: (82,24,73)
1.8mm: (46,14,41)

settings:
0-1300 mean
0-40 hot R2
0-1500 T1

maximize window.

3.00 zoom.
center on cursor

take snapshots.

%%%%% MAKE UPSAMPLED VERSION:

vols = { '~/nsd/nsddata/ppdata/subj01/func1pt8mm/T1_to_func1pt8mm.nii.gz' ...
         '~/nsd/nsddata/ppdata/subj01/func1pt8mm/mean.nii.gz' ...
         '~/nsd/nsddata/ppdata/subj01/func1pt8mm/R2.nii.gz'};
for pp=1:length(vols)
  a1 = load_untouch_nii(vols{pp});
  
  %corner is matched at anterior right inferior
  %niftis are LPI
  % matrixvals = {
  % [145 186 148] [81 104 83] 227021 226601
  % [146 190 150] [82 106 84] 239633 239309
  % [145 190 146] [81 106 82] 240830 243023
  % [152 177 143] [85 99 80] 228495 227262
  % [141 173 139] [79 97 78] 197594 198908
  % [152 202 148] [85 113 83] 253634 259406
  % [139 170 145] [78 95 81] 198770 200392
  % [143 184 139] [80 103 78] 224364 224398
  % };
  
  ii1 = linspacefixeddiff(size(a1.img,1),-1/1.8,145);
  ii2 = linspacefixeddiff(size(a1.img,2),-1/1.8,186);
  ii3 = linspacefixeddiff(1,1/1.8,148);
  [aa1,aa2,aa3] = ndgrid(ii1,ii2,ii3);
  
  f = ba_interp3_wrapper(double(a1.img),[aa1(:) aa2(:) aa3(:)]','cubic');
  
  nsd_savenifti(flipdim(flipdim(reshape(f,size(aa1)),1),2),[1 1 1],sprintf('~/Desktop/vol%d.nii.gz',pp));
end

%% %%%%%%%%%%%%% Make quantitative plots

%% LOAD

dirname = {'func1mm' 'func1pt8mm'};
resval = [1 1.8];
subjix = 1;

meanvols = {};
for rr=1:2
  files0 = matchfiles(sprintf('~/nsd/nsddata/ppdata/subj%02d/%s/mean_session*.nii.gz',subjix,dirname{rr}));
  for zz=length(files0):-1:1, zz
    meanvols{rr}(:,:,:,zz) = single(getfield(load_untouch_nii(files0{zz}),'img'));
  end
end

r2vols = {};
for rr=1:2
  files0 = matchfiles(sprintf('~/nsd/nsddata/ppdata/subj%02d/%s/R2_session*.nii.gz',subjix,dirname{rr}));
  for zz=length(files0):-1:1, zz
    r2vols{rr}(:,:,:,zz) = single(getfield(load_untouch_nii(files0{zz}),'img'));
  end
end

meanvolM = {};
meanvolM{1} = single(getfield(load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj01/func1mm/mean.nii.gz')),'img'));
meanvolM{2} = single(getfield(load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj01/func1pt8mm/mean.nii.gz')),'img'));
meanvolM{3} = single(getfield(load_untouch_nii(sprintf('~/Desktop/vol2.nii.gz')),'img'));

r2volM = {};
r2volM{1} = single(getfield(load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj01/func1mm/R2.nii.gz')),'img'));
r2volM{2} = single(getfield(load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj01/func1pt8mm/R2.nii.gz')),'img'));
r2volM{3} = single(getfield(load_untouch_nii(sprintf('~/Desktop/vol3.nii.gz')),'img'));

%% PROCEED

figureprep([100 100 370*.7 480*1.2],1);

offset = 13;
xx = 82+offset:-1:82-offset;
xxalt = linspacefixeddiff(145,-1.8,size(r2volM{2},1));
%%irange = 82-offset-1:.1:82+offset+1;

subplot(3,1,1); hold on;
% 82,24,73; ITKSNAP has 1-indexed LPI
plot(xx,r2volM{1}(xx,24,73),'ro-','LineWidth',2);                             % 1mm is ro-
%%plot(irange,interp1(xx,r2volM{1}(xx,24,73),irange,'cubic'),'r-');
plot(xxalt,flipud(r2volM{2}(:,14,41)),'cx','LineWidth',2);                    % 1.8mm is co-
plot(xx,r2volM{3}(xx,24,73),'bo-','LineWidth',2);                             % 1.8mm upsampled is bo-
%%plot(irange,interp1(xx,r2volM{3}(xx,24,73),irange,'cubic'),'b-');
set(gca,'XDir','reverse');
ax = axis;
axis([min(xx) max(xx) -5 55]);
ax = axis;
xlabel('Left-to-right position (mm)');
ylabel('Variance explained (R^2)');

subplot(3,1,2); hold on;
for p=1:40
  plot(xx,r2vols{1}(xx,24,73,p),'r-','Color',[1 .7 .7]);
  mn = mean(r2vols{1}(xx,24,73,:),4);
  se = std(r2vols{1}(xx,24,73,:),[],4)/sqrt(40);
  set(errorbar2(xx,mn,se,'v','r-'),'LineWidth',2);
  plot(xx,mn,'r-','LineWidth',2);
end
set(gca,'XDir','reverse');
axis(ax);
xlabel('Left-to-right position (mm)');
ylabel('Variance explained (R^2)');

subplot(3,1,3); hold on;
for p=1:40
  plot(xxalt,flipud(r2vols{2}(:,14,41,p)),'c-','Color',[.7 1 1]);
  mn = mean(r2vols{2}(:,14,41,:),4);
  se = std(r2vols{2}(:,14,41,:),[],4)/sqrt(40);
  set(errorbar2(xxalt,flipud(mn),flipud(se),'v','c-'),'LineWidth',2);
  plot(xxalt,flipud(mn),'c-','LineWidth',2);
end
set(gca,'XDir','reverse');
axis(ax);
xlabel('Left-to-right position (mm)');
ylabel('Variance explained (R^2)');

figurewrite('profiles',[],-2,'~/Desktop/');
