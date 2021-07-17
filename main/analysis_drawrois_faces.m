%%%%% define

% define
subjid = 'subj04';
hemi = 'rh';

%%%%% for each hemisphere

% prep
close all;

% calc
dir0 = sprintf('~/Dropbox/nsddata/freesurfer/%s/',subjid);

% prep
valsLH = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/%s/label/lh.flocfacestval.mgz',subjid));
valsRH = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/%s/label/rh.flocfacestval.mgz',subjid));
vals = cat(1,valsLH,valsRH);
labels = zeros(size(vals));

% perhaps load the existing labels?
if 0
  labels = cvnloadmgz(sprintf('%s/label/*.floc-faces-raw.mgz',dir0));
end

% more
view0 = { {{[10 -70 180] [-10 -70 180]} {'lh' 'rh'}} 'sphere' 0 1000 0 [1 1]};
cvnlookup(subjid,11,   vals,[-10 10],cmapsign4(256),3);  % inflated thresholded
cvnlookup(subjid,view0,vals,[-10 10],cmapsign4(256),3);  % sphere thresholded

%% restart from here for another ROI

% draw ROI (find islands, approximate the t>0 mark but don't have to be super tight)
cvnlookup(subjid,view0,vals,[-10 10],cmapsign4(256));  % sphere unthresholded
pause(2);
Rmask = drawroipoly(himg,Lookup);

% shrink-wrap to 0
Rmask = Rmask & vals > 0;

% make mutual exclusive (in the mask AND not labeled yet)
Rmask = Rmask & labels==0;

% record it
valtouse = 5;
labels = labels + Rmask*valtouse;

% SPECIAL CASE: deletion
if 0
  labels(labels==2 & Rmask) = 0;
end

% SPECIAL CASE: addition
if 0
  labels(Rmask) = 1;
end

%% when done with all the ROIs in this hemisphere

% save file
file0 = sprintf('%s/label/%s.floc-faces-raw.mgz',dir0,hemi);
if isequal(hemi,'lh')
  nsd_savemgz(labels(1:length(valsLH)),    file0,dir0);
else
  nsd_savemgz(labels(length(valsLH)+1:end),file0,dir0);
end

%% when done with both hemispheres

% these images are for reference (if we want to draw again)

% write full image
cvnlookup(subjid,view0,vals,[-10 10],cmapsign4(256),[],[],0);
imwrite(rgbimg,sprintf('~/Dropbox/KKTEMP/sphere/%s_full.png',subjid));

% write ROIs
alllabels = cvnloadmgz(sprintf('%s/label/*.floc-faces-raw.mgz',dir0));
cvnlookup(subjid,view0,alllabels,[0 max(alllabels)],jet,0.5,[],0);
imwrite(rgbimg,sprintf('~/Dropbox/KKTEMP/sphere/%s_rois.png',subjid));

%% when everything looks good, save a nicer version from the raw version

% define
labels = {'IOG' 'pFus' 'mFus' 'mTL' 'aTL'};
groups = {};
groups{1,1} = {1   2   []     [] 3};
groups{1,2} = {1   2   [3 4]  [] 5};
groups{2,1} = {1   2   3      [] 4};
groups{2,2} = {1   2   3      [] 4};
groups{3,1} = {1   2   3      [] 4};
groups{3,2} = {1   2   3      [] 4};
groups{4,1} = {1   2   3      4  []};
groups{4,2} = {1   2   3      4  []};
groups{5,1} = {1   2   3      4  5};
groups{5,2} = {1   2   3      4  5};
groups{6,1} = {1   2   3      [] []};
groups{6,2} = {1   2   3      [] 4};
groups{7,1} = {1   2   3      4  5};
groups{7,2} = {1   2   3      [] 4};
groups{8,1} = {1   2   3      4  5};
groups{8,2} = {1   2   3      4  5};

% perform the union and save the new version
hemis = {'lh' 'rh'};
for subjix=1:8
  dir0 = sprintf('~/Dropbox/nsddata/freesurfer/subj%02d/',subjix);
  for hemi=1:2
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.floc-faces-raw.mgz',subjix,hemis{hemi});
    out0 =  sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.floc-faces.mgz',    subjix,hemis{hemi});
    vals = cvnloadmgz(file0);
    newvals = zeros(size(vals));
    for rr=1:length(labels)
      newvals(ismember(vals,groups{subjix,hemi}{rr})) = rr;
    end
    nsd_savemgz(newvals,out0,dir0);
  end
end

%% create ctab file

for subjix=1:8
  copyfile(sprintf('~/Dropbox/KKTEMP/%s.mgz.ctab','floc-faces'), ...
           sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/',subjix));
end

%%%%% proceed to TRACE OVER AND REGEN

% go through and trace.

zz = 3;

subjid = sprintf('subj%02d',zz);   % which subject
cmap   = jet(256);   % colormap for ROIs
rng    = [0 5];      % should be [0 N] where N is the max ROI index
roilabels = {'OFA' 'FFA-1' 'FFA-2' 'mTL-faces' 'aTL-faces'};  % 1 x N cell vector of strings
mgznames = {'flocfacestval'};%{'curvature'};           % quantities of interest (1 x Q)
crngs = {[-10 10]};%{[-1 1]};          % ranges for the quantities (1 x Q)
cmaps = {cmapsign4(256)};%{copper};   % colormaps for the quantities (1 x Q)
  %roivals = [];                                                                    % start with a blank slate?
  %roivals = cvnloadmgz(sprintf('%s/%s/label/?h.visualsulc.mgz',cvnpath('freesurfer'),subjid));  % load in an existing file?
threshs = {[]};
roivals = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/?h.floc-faces.mgz',zz));

% do it
cvndefinerois;

% save to folder called try2faces/

%%%%%% NEXT

% NOTE: Apr 19 2020, regen from try2faces/

% whittle based on t>0
hemis = {'lh' 'rh'};
for zz=1:8
  for hemi=1:2
    data0 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.flocfacestval.mgz',zz,hemis{hemi}));
    file0b = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.floc-faces.mgz',zz,hemis{hemi});
    file0 = sprintf('~/Desktop/try2faces/%d/%s.floc-faces.mgz',zz,hemis{hemi});
    roi0 = cvnloadmgz(file0);
    for roi=1:5
      roi0(roi0==roi) = (data0(roi0==roi)>0) * roi;
    end
    nsd_savemgz(roi0(:),file0b,sprintf('~/nsd/nsddata/freesurfer/subj%02d/',zz));
  end
end

%%%%% proceed to analysis_drawrois_facesandwords.m
