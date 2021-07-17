% setup
cd /home/stone/generic/Dropbox/nsdfaceprf/
hemis = {'lh' 'rh'};
setfsnsd;

% load stuff
clear rois rps angs eccs preds r2s szs prfangs prfeccs prfr2s prfszs;  % 8 subjects x 2 hemis
for p=1:8
  for hh=1:2
    rois{p,hh}  = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.prf-visualrois.mgz',p,hemis{hh}));
    rps{p,hh}   = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.randomprojection.mgz',p,hemis{hh}));

if 0
    angs{p,hh}  = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrast_angle.mgz',p,hemis{hh}));
    eccs{p,hh}  = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrast_eccentricity.mgz',p,hemis{hh}));
    preds{p,hh} = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrast_predr.mgz',p,hemis{hh}));
    r2s{p,hh}   = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrast_R2.mgz',p,hemis{hh}));
    szs{p,hh}   = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrast_size.mgz',p,hemis{hh}));

    angs{p,hh}  = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastBASELINE_angle.mgz',p,hemis{hh}));
    eccs{p,hh}  = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastBASELINE_eccentricity.mgz',p,hemis{hh}));
    preds{p,hh} = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastBASELINE_predr.mgz',p,hemis{hh}));
    r2s{p,hh}   = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastBASELINE_R2.mgz',p,hemis{hh}));
    szs{p,hh}   = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastBASELINE_size.mgz',p,hemis{hh}));
end

    angs{p,hh}  = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastNEW_angle.mgz',p,hemis{hh}));
    eccs{p,hh}  = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastNEW_eccentricity.mgz',p,hemis{hh}));
    preds{p,hh} = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastNEW_predr.mgz',p,hemis{hh}));
    r2s{p,hh}   = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastNEW_R2.mgz',p,hemis{hh}));
    szs{p,hh}   = cvnloadmgz(sprintf('freesurfer/subj%02d/label/%s.contrastNEW_size.mgz',p,hemis{hh}));

    prfangs{p,hh}  = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.prfangle.mgz',p,hemis{hh}));
    prfeccs{p,hh}  = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.prfeccentricity.mgz',p,hemis{hh}));
    prfr2s{p,hh}   = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.prfR2.mgz',p,hemis{hh}));
    prfszs{p,hh}   = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.prfsize.mgz',p,hemis{hh}));
  end
end

%%%%%%% plot maps

% setup
mkdirquiet('~/Dropbox/KKTEMP/maps');
cd ~/Dropbox/KKTEMP/maps
viewpoint = 13; % native flat
viewpoint = 1;  % native sphere

% make maps
for p=1:8
  Lookup = [];
  eo = {'hemibordercolor' 'w' 'rgbnan' 1 'drawscalebar' false};

  cvnlookup(sprintf('subj%02d',p),viewpoint,posrect(cat(1,prfr2s{p,:})),[0 100],hot,[],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_prfR2.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,mod(180-prfangs{p,1},360),prfangs{p,2}),[0 360],cmfanglecmapRH,[],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_prfang.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,prfeccs{p,:}),[0 12], cmfecccmap(4), [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_prfecc.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,prfszs{p,:}), [0 6],  cmfecccmap(7), [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_prfsize.png',p));

  cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,preds{p,:}),  [0 1],  hot,           [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_nsdpredr.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,posrect(cat(1,r2s{p,:})),   [0 50],hot,[],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_nsdR2.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,mod(180-angs{p,1},360),angs{p,2}),[0 360],cmfanglecmapRH,[],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_nsdang.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,eccs{p,:}),   [0 12], cmfecccmap(4), [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_nsdecc.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,szs{p,:}),    [0 6],  cmfecccmap(7), [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_nsdsize.png',p));

  cvnlookup(sprintf('subj%02d',p),viewpoint,[],[],[],                               10000,Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_curvature.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,rois{p,:}),  [0 7],  jet(256),       [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_rois.png',p));
  cvnlookup(sprintf('subj%02d',p),viewpoint,subscript(cat(1,rps{p,:}),{':' 1 1 1}),   [-.04 .04],cmapsign4,   [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_randomproj.png',p));

      %   cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,prfangs{p,:}),[0 360],cmfanglecmapRH,[],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_prfangLVF.png',p));
      %   cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,prfangs{p,:}),[0 360],cmfanglecmap,  [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_prfangRVF.png',p));
      %   cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,angs{p,:}),   [0 360],cmfanglecmapRH,[],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_nsdangLVF.png',p));
      %   cvnlookup(sprintf('subj%02d',p),viewpoint,cat(1,angs{p,:}),   [0 360],cmfanglecmap,  [],Lookup,0,eo); imwrite(rgbimg,sprintf('subj%02d_nsdangRVF.png',p));

end

%%%%%%% other visualizations

% 0 Unknown
% 1 V1v	    
% 2 V1d	   
% 3 V2v	   
% 4 V2d	   
% 5 V3v	    
% 6 V3d	    
% 7 hV4

% 0 Unknown
% 1 ecc0pt5
% 2 ecc1
% 3 ecc2
% 4 ecc4
% 5 ecc4+

% plot visual field center estimates (NSD)
thresh = 0;
for hh=1:2
  figureprep([0 0 4000 3000]); hold on;
  for p=1:8     % 8 subjects (columns)
    for r=1:6   % 6 ROIs (rows)
      subplot(6,8,(r-1)*8+p); hold on;
      ix1 = find(rois{p,hh} == r & r2s{p,hh} > thresh);
      ix2 = find(rois{p,hh} == r & r2s{p,hh} <= thresh);
      yval = sin(angs{p,hh}/180*pi) .* eccs{p,hh};
      xval = cos(angs{p,hh}/180*pi) .* eccs{p,hh};
      axis([-8 8 -8 8]);
      axis square;
      straightline(0,'h','k-');
      straightline(0,'v','k-');
      scatter(xval(ix2),yval(ix2),16,'b.');
      scatter(xval(ix1),yval(ix1),16,'r.');
      drawrectangle(0,0,8.4,8.4,'k-');
      drawellipse(0,0,0,4.2,4.2,[],[],'k-');
    end
  end
  figurewrite(sprintf('nsdcenters%d',hh),[],-2,'~/Dropbox/KKTEMP');
end

% plot visual field center estimates (PRF)
thresh = 0;
for hh=1:2
  figureprep([0 0 4000 3000]); hold on;
  for p=1:8     % 8 subjects (columns)
    for r=1:6   % 6 ROIs (rows)
      subplot(6,8,(r-1)*8+p); hold on;
      ix1 = find(rois{p,hh} == r & r2s{p,hh} > thresh);
      ix2 = find(rois{p,hh} == r & r2s{p,hh} <= thresh);
      yval = sin(prfangs{p,hh}/180*pi) .* prfeccs{p,hh};
      xval = cos(prfangs{p,hh}/180*pi) .* prfeccs{p,hh};
      axis([-8 8 -8 8]);
      axis square;
      straightline(0,'h','k-');
      straightline(0,'v','k-');
      scatter(xval(ix2),yval(ix2),16,'b.');
      scatter(xval(ix1),yval(ix1),16,'r.');
      drawrectangle(0,0,8.4,8.4,'k-');
      drawellipse(0,0,0,4.2,4.2,[],[],'k-');
    end
  end
  figurewrite(sprintf('prfcenters%d',hh),[],-2,'~/Dropbox/KKTEMP');
end

% scatter one vs other (ang)
thresh = 0;
for hh=1:2
  figureprep([0 0 4000 3000]); hold on;
  for p=1:8     % 8 subjects (columns)
    for r=1:6   % 6 ROIs (rows)
      subplot(6,8,(r-1)*8+p); hold on;
      ix1 = find(rois{p,hh} == r & r2s{p,hh} > thresh);
      ix2 = find(rois{p,hh} == r & r2s{p,hh} <= thresh);
      scatter(angs{p,hh}(ix2),prfangs{p,hh}(ix2),16,'b.');
      scatter(angs{p,hh}(ix1),prfangs{p,hh}(ix1),16,'r.');
      axis([0 360 0 360]);
      axis square;
      axissquarify;
      axis([0 360 0 360]);
    end
  end
  figurewrite(sprintf('compareang%d',hh),[],-2,'~/Dropbox/KKTEMP');
end

% scatter one vs other (ecc)
thresh = 0;
for hh=1:2
  figureprep([0 0 4000 3000]); hold on;
  for p=1:8     % 8 subjects (columns)
    for r=1:6   % 6 ROIs (rows)
      subplot(6,8,(r-1)*8+p); hold on;
      ix1 = find(rois{p,hh} == r & r2s{p,hh} > thresh);
      ix2 = find(rois{p,hh} == r & r2s{p,hh} <= thresh);
      scatter(eccs{p,hh}(ix2),prfeccs{p,hh}(ix2),16,'b.');
      scatter(eccs{p,hh}(ix1),prfeccs{p,hh}(ix1),16,'r.');
      axis([0 8 0 8]);
      axis square;
      axissquarify;
      axis([0 8 0 8]);
    end
  end
  figurewrite(sprintf('compareecc%d',hh),[],-2,'~/Dropbox/KKTEMP');
end

========================

% duplicate maps to mapscrop!!!

% cd to mapscrop directory.
cd('/research/papers/2023 nsdfaceprf/analysis/PARAMETERS/mapscrop');

% crop the two hemispheres and concatenate
thecrop = cell(1,8);
thecrop{1} = imcropfiles('subj01*.png',thecrop{1});
thecrop{2} = imcropfiles('subj02*.png',thecrop{2});
thecrop{3} = imcropfiles('subj03*.png',thecrop{3});
thecrop{4} = imcropfiles('subj04*.png',thecrop{4});
thecrop{5} = imcropfiles('subj05*.png',thecrop{5});
thecrop{6} = imcropfiles('subj06*.png',thecrop{6});
thecrop{7} = imcropfiles('subj07*.png',thecrop{7});
thecrop{8} = imcropfiles('subj08*.png',thecrop{8});
save('thecrop.mat','thecrop');

  %1     {'subj01_curvature.png' }
  %2     {'subj01_nsdR2.png'     }
  %3     {'subj01_nsdang.png'    }
  %4     {'subj01_nsdecc.png'    }
  %5     {'subj01_nsdpredr.png'  }
  %6     {'subj01_nsdsize.png'   }
  %7     {'subj01_prfR2.png'     }
  %8     {'subj01_prfang.png'    }
  %9     {'subj01_prfecc.png'    }
  %10    {'subj01_prfsize.png'   }
  %11     {'subj01_randomproj.png'}
  %12     {'subj01_rois.png'      }
% consolidate images into a large column (for each subject)
for p=1:8
  imwrite(concatimages(sprintf('subj%02d*.png',p),{[8 3 9 4 7 2 1]},1),sprintf('../allsubj%02d.png',p));
end

% make an alpha with the borders of the ROIs
for p=1:8
  aa = double(imread(sprintf('subj%02d_rois.png',p)));
  bb = double(detectedges(  aa(:,:,1)/255*100 + aa(:,:,2)/255*10 + aa(:,:,3)  ,.1) > 0.1);  % binary, 0/1
  bb = repmat(bb,[7 1]);
  imwrite(uint8(255*bb),sprintf('../allsubj%02d_rois.png',p),'Alpha',double(bb));  % white borders with transparency
end
