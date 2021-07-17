% inherit from analysis_drawrois_facesandwords.m

% these images are our final inspections and quality assessments

%%%%% first, go back to _part02.m and do some manual ones:

do the following:
angle (4), q
ecc (5), o, 
kastner (3), q
eccroi (9), q

filenames:
subj01_prf-visualrois_on_angle.png
subj01_prf-eccrois_on_eccentricity.png
subj01_prf-visualrois_on_Kastner2015.png
subj01_prf-visualrois_on_eccrois.png

%%%%% next, do some tailored ones

outputdir = '~/Dropbox/KKTEMP/prf-visualandecc';
mkdirquiet(outputdir);

todo = {{'prf-visualrois'} {'prf-eccrois'}};
for zz=1:length(todo)
  strs = todo{zz};

  for subjix=1:8

    % load visualsulc
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.visualsulc.mgz',subjix);
    file0 = cvnloadmgz(file0);

    % write visualsulc
    cvnlookup(sprintf('subj%02d',subjix),1,file0,[0 14],jet(256),0.5,[],0);
    imwrite(rgbimg,sprintf('%s/subj%02d_visualsulc.png',outputdir,subjix));

    % load ROIs
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.%s.mgz',subjix,strs{1});
    file0 = cvnloadmgz(file0);

    % write ROIs
    extraopts = {'roiname',strs{1},'roiwidth',0,'drawroinames',true};
    cvnlookup(sprintf('subj%02d',subjix),1,file0,[0 7],jet,0.5,[],0,extraopts);
    imwrite(rgbimg,sprintf('%s/subj%02d_%s.png',outputdir,subjix,strs{1}));
      % more (fsaverage inflated)
    cvnlookup(sprintf('subj%02d',subjix), ...
      {'occip' 'sphere' 0 1000 1 [1 1]}, ...
      file0,[0 7],jet,0.5,[],0);
    imwrite(rgbimg,sprintf('%s/fsaverage_subj%02d_%s.png',outputdir,subjix,strs{1}));

  end

end

% let's look at both at the same time!
for subjix=1:8
  extraopts = {'roiname',{'prf-visualrois' 'prf-eccrois'},'roicolor',{'r' 'b'},'drawroinames',true};
  cvnlookup(sprintf('subj%02d',subjix),1,[],[],[],100,[],0,extraopts);
  imwrite(rgbimg,sprintf('%s/subj%02d_visualandecc.png',outputdir,subjix));
end

%% then proceed to convert surface to volume (examples_surfacetovolume.m)
