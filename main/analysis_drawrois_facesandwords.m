%% when done with all subjects, write out inflated views

% these images are our final inspections and quality assessments

  % OBSOLETE
  %mkdirquiet('~/Dropbox/KKTEMP/floc-faces-raw');

outputdir = '~/Dropbox/KKTEMP/floc-facesandwords';
mkdirquiet(outputdir);

todo = {{'floc-faces' 'flocfacestval'} {'floc-words' 'flocwordtval'}};
for zz=1:length(todo)
  strs = todo{zz};

  for subjix=1:8

    % load visualsulc
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.visualsulc.mgz',subjix);
    file0 = cvnloadmgz(file0);

    % write visualsulc
    cvnlookup(sprintf('subj%02d',subjix),11,file0,[0 14],jet(256),0.5,[],0);
    imwrite(rgbimg,sprintf('%s/subj%02d_visualsulc.png',outputdir,subjix));

    % load tvalues
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.%s.mgz',subjix,strs{2});
    file0 = cvnloadmgz(file0);

    % write tvalue
    cvnlookup(sprintf('subj%02d',subjix),11,file0,[-10 10],cmapsign4(256),[],[],0);
    imwrite(rgbimg,sprintf('%s/subj%02d_%s.png',outputdir,subjix,strs{2}));

    % OBSOLETE
    %   % load ROIs
    %   file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.floc-faces-raw.mgz',subjix);
    %   file0 = cvnloadmgz(file0);
    % 
    %   % write ROIs
    %   cvnlookup(sprintf('subj%02d',subjix),11,file0,[0 max(file0)],jet,0.5,[],0);
    %   imwrite(rgbimg,sprintf('~/Dropbox/KKTEMP/floc-faces-raw/subj%02d_roi.png',subjix));

    % load ROIs
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.%s.mgz',subjix,strs{1});
    file0 = cvnloadmgz(file0);

    % write ROIs
    extraopts = {'roiname',strs{1},'roiwidth',0,'drawroinames',true};
    cvnlookup(sprintf('subj%02d',subjix),11,file0,[0 5],jet,0.5,[],0,extraopts);
    imwrite(rgbimg,sprintf('%s/subj%02d_%s.png',outputdir,subjix,strs{1}));
      % more (fsaverage inflated)
    cvnlookup(sprintf('subj%02d',subjix), ...
      {'ventral-lateral' 'inflated' 1 1000 1 [1 1]}, ...
      file0,[0 5],jet,0.5,[],0);
    imwrite(rgbimg,sprintf('%s/fsaverage_subj%02d_%s.png',outputdir,subjix,strs{1}));

  end

end

% let's look at both at the same time!
for subjix=1:8
  extraopts = {'roiname',{'floc-faces' 'floc-words'},'roicolor',{'r' 'b'},'drawroinames',true};
  cvnlookup(sprintf('subj%02d',subjix),11,[],[],[],100,[],0,extraopts);
  imwrite(rgbimg,sprintf('%s/subj%02d_facesandwords.png',outputdir,subjix));
end

%% then proceed to convert surface to volume (analysis_surfaceroistovolume.m)
