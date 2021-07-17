%% %%%%%%%% Generate maps of prf results

% define
resultsdir = [nsd_datalocation '/ppdata/subj%02d/func1mm/'];
files = ...
...% name            fname            mgz name          color range    color map       transform function           circular
  { ...
   {'R2'            'prfR2'           'prfR2'           [0 1]          hot(256)        @(x)(posrect(x)/100).^.5     0} ...
   {'angle'         'prfangleLVF'     'prfangle'        [0 360]        cmapangLVF      @(x)x                        1} ...
   {'angle'         'prfangleRVF'     'prfangle'        [0 360]        cmapangRVF      @(x)x                        1} ...
   {'eccentricity'  'prfeccentricity' 'prfeccentricity' [0 6].^.5      jet(256)        @(x)posrect(x).^.5           0} ...
   {'exponent'      'prfexponent'     'prfexponent'     [0 1]          copper(256)     @(x)x                        0} ...
   {'gain'          'prfgain'         'prfgain'         [-8 8]         cmapsign4(256)  @(x)x                        0} ...
   {'meanvol'       'prfmeanvol'      'prfmeanvol'      [0 1]          gray(256)       @(x)(posrect(x)/4095).^.5    0} ...
   {'size'          'prfsize'         'prfsize'         [0 4].^.5      jet(256)        @(x)posrect(x).^.5           0} ...
  };
hemis = {'lh' 'rh'};
fsdir = [nsd_datalocation '/freesurfer/subj%02d'];
nsubj = 8;
ndepth = 3;
viewswanted = [1 3];
viewnames = {'sphere-occip' 'inflated-ventral'};

%%% internal constants
%%rawniifun = @(x) flipdim(rotatematrix(x,1,2,1),2);

% do it
for ss=1:nsubj

  for vv=1:length(files)
  
    % load data
    sourcedata = sprintf('%s/prf_%s.nii.gz',sprintf(resultsdir,ss),files{vv}{1});
    a1 = getfield(load_untouch_nii(sourcedata),'img');

    % map to surface
    alldata = [];
    for hh=1:length(hemis)

      % map data to subject-native surfaces using linear and average across depth
      data1 = [];
      for ii=1:ndepth
        if files{vv}{7}
          data1(:,:,ii,1) = nsd_mapdata(ss,'func1pt0',sprintf('%s.layerB%d',hemis{hh},ii),cos(a1/180*pi),'linear',[],[],'single');
          data1(:,:,ii,2) = nsd_mapdata(ss,'func1pt0',sprintf('%s.layerB%d',hemis{hh},ii),sin(a1/180*pi),'linear',[],[],'single');
        else
          data1(:,:,ii) = nsd_mapdata(ss,'func1pt0',sprintf('%s.layerB%d',hemis{hh},ii),a1,'linear',[],[],'single');
        end
      end
      if files{vv}{7}
        data1 = mod(atan2(mean(data1(:,:,:,2),3),mean(data1(:,:,:,1),3))/pi*180,360);
      else
        data1 = mean(data1,3);
      end

      % save official .mgz file
      outputfile = sprintf([nsd_datalocation '/freesurfer/subj%02d/label/%s.%s.mgz'],ss,hemis{hh},files{vv}{3});
      nsd_savemgz(data1,outputfile,sprintf(fsdir,ss));

      % record
      alldata = cat(1,alldata,data1);

    end
    if vv == 1
      r2 = alldata;
      r2thresh = findtailthreshold(r2(:),0);
    end

    for viewix=1:length(viewswanted)

      % unthresholded map
      outputfile =    sprintf([nsd_datalocation '/figures/subj%02d-%s_%s.png'],      ss,viewnames{viewix},files{vv}{2});
%%      outputfileNII = sprintf(['~/nsd/ppdata/figuresnii/subj%02d-%s_%s.nii.gz'],ss,viewnames{viewix},files{vv}{2});
      [rawimg,Lookup,rgbimg,himg] = cvnlookup(sprintf('subj%02d',ss),viewswanted(viewix),feval(files{vv}{6},alldata),files{vv}{4},files{vv}{5},[],[],0);
      imwrite(rgbimg,outputfile);
%%      nsd_savenifti(rawniifun(rawimg),[1 1 1],outputfileNII);

      % thresholded map
      outputfile = sprintf([nsd_datalocation '/figures/subj%02d-%s_%s_thresh.png'],ss,viewnames{viewix},files{vv}{2});
      [rawimg,Lookup,rgbimg,himg] = cvnlookup(sprintf('subj%02d',ss),viewswanted(viewix),feval(files{vv}{6},alldata),files{vv}{4},files{vv}{5},[],[],0,{'overlayalpha',r2 > r2thresh});
      imwrite(rgbimg,outputfile);
      
%       % curvature
%       if vv==1
%         temp = cvnreadsurfacemetric(sprintf('subj%02d',ss),[],'curv','white','orig');
% %%        outputfileNII = sprintf(['~/nsd/ppdata/figuresnii/subj%02d-%s_%s.nii.gz'],ss,viewnames{viewix},'curvature');
%         [rawimg,Lookup,rgbimg,himg] = cvnlookup(sprintf('subj%02d',ss),viewswanted(viewix),double(temp.data < 0),[-1 2],gray(256),[],[],0);
% %%        nsd_savenifti(rawniifun(rawimg),[1 1 1],outputfileNII);
%       end

    end
    
  end
  
end

% touch
unix_wrapper(sprintf('ssh kendrick@stone.cmrr.umn.edu touch %s',[nsd_datalocation '/figures']));
