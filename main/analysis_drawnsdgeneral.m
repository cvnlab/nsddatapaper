%% This document shows the manual drawing of a "nsdgeneral" ROI on the fsaverage surface.
%% This atlas ROI is saved as lh and rh .mgz files to the fsaverage/label directory.


% NOTE: this relies on the creation of *.b3R2.mgz which is done by junk/examples_funcviz.m


%% Plot R2 values

% Load surface-based R2 data (created by examples_funcviz.m). These values 
% reflect variance explained by the 'betas_fithrf_GLMdenoise_RR' GLM model
% and are relative to the fsaverage surface.
inputfile = sprintf('~/nsd/ppdata/??.b3R2.mgz');
data = cvnloadmgz(inputfile);  % 327684 x 1 x 1 x [8 subjects + 1 group-average = 9]

% Look at group-average R2 values on the fsaverage sphere occip view.
data0 = (posrect(data(:,1,1,end))/100).^.5;
[rawimg,Lookup,rgbimg,himg] = cvnlookup('fsaverage',1,data0,[0 1],hot(256),(20/100).^.5);


%% Draw and save ROIs

% Note: repeat this process for the left and right hemispheres!

% Draw the ROI
Rmask = drawroipoly(himg,Lookup);  % make sure we are drawing on the spherical surface!

% Save out .mgz file into the fsaverage/label directory
hemi = 'lh';
valstruct = valstruct_create('fsaverage');
valstruct = setfield(valstruct,'data',Rmask);
fsdir = [nsd_datalocation '/freesurfer/fsaverage'];
outputfile = sprintf([nsd_datalocation '/freesurfer/fsaverage/label/%s.nsdgeneral.mgz'],hemi);
nsd_savemgz(valstruct_getdata(valstruct,hemi),outputfile,fsdir);

% Also create a .ctab file to document what these regions are:
%   freesurfer/fsaverage/label/nsdgeneral.mgz.ctab


% obsolete. will be handled by analysis_surfacevisualizations.m
% %% Make a visualization
% 
% % Using the files from above, save out a .png visualization
% inputfile = sprintf([nsd_datalocation '/freesurfer/fsaverage/label/lh.nsdgeneral.mgz']);
% outputfile = sprintf([nsd_datalocation '/figures/fsaverageflat_nsdgeneral_fsaverage.png'])
% [rawimg,Lookup,rgbimg,himg] = cvnlookup('fsaverage',10,inputfile,[0 1],gray,0.5);
% imwrite(rgbimg,outputfile);
% 
