% The goal here is to generate basic surface visualizations
% relevant to the NSD dataset.

%% DEFINE

  %          name           color range   color map           threshold outname
atlasesMASTER = { ...
            {'curvature'    []            []                  1000      []} ...
            {'Kastner2015'  [0 25]        jet(256)            0.5       []} ...
            {'HCP_MMP1'     [-.5 180.5]   colormap_hcp_mmp()  0.5       []} ...
            {'corticalsulc' [0 28]        jet(256)            0.5       []} ...
            {'streams'      [0 7]         jet(256)            0.5       []} ...
            {'nsdgeneral'   [0 1]         copper(256)           0.5       []} ...
            {'gEVC'         [0 1]         copper(256)           0.5       []} ...
            {'gVTC'         [0 1]         copper(256)           0.5       []} ...
         };
hemis = {'lh' 'rh'};
outputdir = '~/nsd/nsddata/inspections/surfacevisualizations/'; mkdirquiet(outputdir);
extraopts = {'rgbnan',1,'hemibordercolor',[1 1 1],'scalebarcolor','k'};  %% NOTE
viewswanted = [  13     3                 5                6                 14];
viewnames =   { 'flat' };   %'inflatedventral' 'inflatedmedial' 'inflatedlateral' 'inflatedsuperior'};

%% PLOT CURVATURE+ATLASES ON THE FSAVERAGE SURFACE

subjectid = 'fsaverage';
dataid = subjectid;
atlases = atlasesMASTER;
analysis_surfacevisualizations_helper;

%% PLOT CURVATURE+CORTICALSULC ON SUBJECT-NATIVE SURFACES

for p=1:8
  subjectid = sprintf('subj%02d',p);
  dataid = subjectid;
  atlases = atlasesMASTER([1 5]);
  analysis_surfacevisualizations_helper;
end

%% PLOT SUBJECT-NATIVE RESULTS ON THE FSAVERAGE SURFACE

% For individual NSD subjects, values are automatically transferred (by cvnlookup.m) using 
% nearest-neighbor to the fsaverage surface.

todo = { ...
%  name        fileloc                                                                                     transform function                 color range   color map    threshold   outname
  {'curvature' ['~/nsd/nsddata/freesurfer/subj%02d/*/%s.%s.mgz']                                           @(x)x<0                            [-1 2]        gray(256)    []          []} ...
  {'mean'      ['~/nsd/nsddata/freesurfer/subj%02d/*/%s.%s.mgz']                                           @(x)x                              [0 2000]      gray(256)    []          []} ...
  {'valid'     ['~/nsd/nsddata/freesurfer/subj%02d/*/%s.%s.mgz']                                           @(x)x                              [0 1]         jet(256)     []          []} ...
  {'R2'        ['~/nsd/nsddata/freesurfer/subj%02d/*/%s.%s.mgz']                                           @(x)signedarraypower(x/100,0.5)    [0 1]         hot(256)     []          []} ...
  {'b3R2'      ['~/nsd/nsddata/freesurfer/subj%02d/*/%s.%s.mgz']                                           @(x)signedarraypower(x/100,0.5)    [0 1]         hot(256)     []          []} ...
  {'surfaceimperfections' ['~/nsd/nsddata/freesurfer/subj%02d/*/%s.%s.mgz']                                @(x)x>0.01                         [0 1]         winter(256)  1/8/2       []} ...
  {'signaldropout'        ['~/nsd/nsddata/freesurfer/subj%02d/*/%s.%s.mgz']                                @(x)x>2                            [0 1]         winter(256)  1/8/2       []} ...
  {'ncsnr'     ['~/nsd/nsddata_betas/ppdata/subj%02d/nativesurface/betas_assumehrf/%s.%s.mgh']             @(x)100*(x.^2./(x.^2+1/3))         [0 75]        jet(256)     []          'b1nc'} ...
  {'ncsnr'     ['~/nsd/nsddata_betas/ppdata/subj%02d/nativesurface/betas_fithrf/%s.%s.mgh']                @(x)100*(x.^2./(x.^2+1/3))         [0 75]        jet(256)     []          'b2nc'} ...
  {'ncsnr'     ['~/nsd/nsddata_betas/ppdata/subj%02d/nativesurface/betas_fithrf_GLMdenoise_RR/%s.%s.mgh']  @(x)100*(x.^2./(x.^2+1/3))         [0 75]        jet(256)     []          'b3nc'} ...
};

for tt=1:length(todo)
  name =    todo{tt}{1};
  fileloc = todo{tt}{2};
  tfun =    todo{tt}{3};
  rng =     todo{tt}{4};
  cmap =    todo{tt}{5};
  thresh =  todo{tt}{6};
  outname = todo{tt}{7};

  % load data and map to fsaverage
  vals = {};
  for p=1:8
    for hh=1:length(hemis)
      inputfile = matchfiles(sprintf(fileloc,p,hemis{hh},name));
      assert(length(inputfile)==1);
      inputfile = inputfile{1};
      data = cvnloadmgz(inputfile);
      vals{p,hh} = nsd_mapdata(p,[hemis{hh} '.white'],'fsaverage',data);
    end
  end

  % plot subject data on fsaverage
  for p=1:8
    subjectid = 'fsaverage';
    dataid = sprintf('subj%02d',p);
    atlases = {{{name tfun(catcell(1,vals(p,:)))} rng cmap thresh outname}};
    analysis_surfacevisualizations_helper;
  end
  
  % plot group-average data on fsaverage
  subjectid = 'fsaverage';
  dataid = 'groupavg';
  atlases = {{{name vflatten(mean(tfun(reshape(cell2mat(vals),[],8,2)),2))} rng cmap thresh outname}};
  analysis_surfacevisualizations_helper;

end
