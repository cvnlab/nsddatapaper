% inherited from analysis_prf_maps.m
%
% this script just makes a bunch of maps.

%% %%%%%%%% Prep special colormaps

cmfecccmapA = cmfecccmap(4);  % to green
cmfecccmapB = cmfecccmap(4);  % to green
cmfecccmapC = cmfecccmap(7);  % to cyan

vals = linspace(0,12,257);
vals = vals(1:256) + diff(vals(1:2))/2;
whval =  [0.5 1 2];
whband = [  1 2 2];
keepix = [];
for p=1:length(whval)
  [~,ix] = min(abs(vals-whval(p)));
  keepix = [keepix ix+(-whband(p):whband(p))];
%  keepix = [keepix; find(abs(vals-whval(p)) < whval(p)*0.2)];
end
cmfecccmapB(setdiff(1:256,keepix),:) = 0;

figure; colormap(cmfecccmapA); colorbar;
figure; colormap(cmfecccmapB); colorbar;
figure; colormap(cmfecccmapC); colorbar;

%% %%%%%%%% Generate maps of prf results

% define
files = ...
...%  png name          mgz name        color range    color map       transform function           threshold
  { ...
   {'prfR2'            'prfR2'           [0 1]          hot(256)        @(x)(posrect(x)/100).^.5     []} ...
   {'prfangleLVF'      'prfangle'        [0 360]        cmfanglecmapRH  @(x)x                        []} ...
   {'prfangleRVF'      'prfangle'        [0 360]        cmfanglecmap    @(x)x                        []} ...
   {'prfeccentricity'  'prfeccentricity' [0 12]         cmfecccmapA     @(x)x                        []} ...
   {'prfeccentricityB' 'prfeccentricity' [0 12]         cmfecccmapB     @(x)x                        []} ...
   {'prfexponent'      'prfexponent'     [0 1]          copper(256)     @(x)x                        []} ...
   {'prfgain'          'prfgain'         [-8 8]         cmapsign4(256)  @(x)x                        []} ...
   {'prfmeanvol'       'prfmeanvol'      [0 1]          gray(256)       @(x)(posrect(x)/4095).^.5    []} ...
   {'prfsize'          'prfsize'         [0 6]          cmfecccmapC     @(x)x                        []} ...
   {'Kastner2015'      'Kastner2015'     [0 25]         jet(256)        @(x)x                        0.5} ...
   {'visualsulc'       'visualsulc'      [0 14]         jet(256)        @(x)x                        0.5} ...
  };
hemis = {'lh' 'rh'};
fsdir = [nsd_datalocation '/freesurfer/subj%02d'];
nsubj = 8;
viewswanted = [1 3 4 13];                                                         % refers to order in cvnlookup.m
viewnames = {'sphere-occip' 'inflated-ventral' 'inflated-parietal' 'full-flat'};  % filenames for these views [how about fsaverage flat??]
outputdir = '~/Dropbox/KKTEMP/prffigures'; mkdirquiet(outputdir);

% loop over subject
for ss=1:nsubj

  % loop over file
  for vv=1:length(files)    %[1 9]%4 5]% 10 11]  %

    % calc
    pngname = files{vv}{1};
    mgzname = files{vv}{2};
    crng =    files{vv}{3};
    cmap =    files{vv}{4};
    tfun =    files{vv}{5};
    thresh =  files{vv}{6};
  
    % load data
    alldata = [];
    for hh=1:length(hemis)
      file0 = sprintf([nsd_datalocation '/freesurfer/subj%02d/label/%s.%s.mgz'],ss,hemis{hh},mgzname);
      alldata = [alldata; cvnloadmgz(file0)];
    end
    
    % special handling of R2 (compute an R2 mask)
    if vv == 1     
      r2 = alldata;
      r2thresh = findtailthreshold(r2(:),0);   %% plot MORE threshold levels?
      r2mask = r2 > r2thresh;
      %continue;
    end
    
    % special stuff for atlases
    if ismember(mgzname,{'Kastner2015' 'visualsulc'})
      roimask = {};
      for hh=1:length(hemis)
        [roimask{hh}] = cvnroimask(sprintf('subj%02d',ss),hemis{hh},[mgzname '*']);
      end
      roimask = cellfun(@(x,y) [x;y],roimask{1},roimask{2},'UniformOutput',0);
      roicolors = {};
      for q=crng(1):crng(2)
        roicolors{end+1} = cmaplookup(q,crng(1),crng(2),[],cmap);
      end
      opts = {'roimask',roimask,'roicolor',roicolors,'overlayalpha',0};  %%,'drawroinames',true ? roidescription0??
    else
      opts = {};
    end
    
    % loop over view
    for viewix=1:length(viewswanted)

      % unthresholded map (OR atlas border-only)
      outputfile = sprintf([outputdir '/subj%02d-%s_%s.png'],ss,viewnames{viewix},pngname);
      cvnlookup(sprintf('subj%02d',ss),viewswanted(viewix),feval(tfun,alldata),crng,cmap,[],[],0,opts);
      imwrite(rgbimg,outputfile);

      % thresholded map
      outputfile = sprintf([outputdir '/subj%02d-%s_%s_thresh.png'],ss,viewnames{viewix},pngname);
      if isempty(thresh)  % if empty, use special r2mask
        cvnlookup(sprintf('subj%02d',ss),viewswanted(viewix),feval(tfun,alldata),crng,cmap,[],[],0,{'overlayalpha',r2mask});
      else                % if non-empty, use that value as the threshold
        cvnlookup(sprintf('subj%02d',ss),viewswanted(viewix),feval(tfun,alldata),crng,cmap,thresh,[],0);
      end
      imwrite(rgbimg,outputfile);
      
    end  % view
    
  end  % file
  
end  % subject
