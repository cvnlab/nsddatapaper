% this is a one-off to create a really nice big figure.
% goal: fsaverage. all sessions, all subjects. mean and R2 (from onoff)

% NOTE: this was all achieved BEFORE we reset the cvnlookup cache. thus, the annotations have to be updated too.

%% %%%%% Map various 1mm volume results to subject surfaces and then to fsaverage.

% define
dir0 = sprintf('%s/ppdata/subj%%02d/func1mm/',nsd_datalocation);
files =     {'mean_session%02d' 'R2_session%02d'};
hemis = {'lh' 'rh'};
nsubj = 8;
ndepth = 3;

% init alldata (cell matrix of dimensions files x subjects, each with 163k vertices x SESSION x 2 hemis)
alldata = {};

% do it
for ff=1:length(files)       % for each file
  for ss=1:nsubj           % for each subject

    % load the volume
    cnt = 1;
    a1 = [];  % X Y Z SESSION
    while 1
      sourcedata = sprintf('%s/%s.nii.gz',sprintf(dir0,ss),sprintf(files{ff},cnt));
      if ~exist(sourcedata,'file')
        break;
      end
      a1 = cat(4,a1,getfield(load_untouch_nii(sourcedata),'img'));
      cnt = cnt + 1;
    end

    for hh=1:length(hemis)     % for each hemisphere

      % map data to surfaces (cubic, conversion to single) and average across depth.
      data1 = [];  % V x SESSION
      for ii=1:ndepth
        data1(:,:,ii) = nsd_mapdata(ss,'func1pt0',sprintf('%s.layerB%d',hemis{hh},ii),a1,'cubic',[],[],'single');
      end
      data1 = mean(data1,3);  % there could be NaNs!

      % map to fsaverage
      alldata{ff,ss}(:,:,hh) = nsd_mapdata(ss,sprintf('%s.white',hemis{hh}),'fsaverage',data1,[],[],[],'single');

    end
    
  end
end

%% %%%%% Generate visualizations

% note: mimic autoqc_glm_nsdBASIC.m

% define
filenames = {'meansession' 'R2session'};
clims = {[0 1] [0 1]};
cmaps = {gray(256) hot(256)};
datafun = {@(x)normalizerange(x,0,1,prctile(x(:),1),prctile(x(:),99),1) @(x)signedarraypower(x/100,0.5)};
outputdir = '~/Dropbox/KKTEMP/grandviz';
mkdirquiet(outputdir);

% do it
for ff=1:length(files)
  for ss=1:8
    for p=1:size(alldata{ff,ss},2)

      % get and transform data
      data = vflatten(alldata{ff,ss}(:,p,:));
      if ~isempty(datafun{ff})
        data = feval(datafun{ff},data);
      end

      % visualize it
      outputfile = sprintf('%s/%s_subj%02d_sess%02d.png',outputdir,filenames{ff},ss,p);
      [rawimg,Lookup,rgbimg,himg] = cvnlookup('fsaverage',10,data,clims{ff},cmaps{ff},[],[],0,{'rgbnan',0.2,'hemibordercolor',[.2 .2 .2]});  % dark background
      imwrite(rgbimg,outputfile);

    end
  end
end

%% %%%%% Prettify with annotations


% MAKE SURE TO USE fsaverageflatannotationsBEFORERESET.png



cd ~/Dropbox/KKTEMP/grandviz

% define
todo = {'R2session' 'meansession'};
colors = {[255 255 255] [0 255 255]};
nsess = [40 40 32 30 40 32 40 30];

% read annotations
%[~,~,alpha] = imread('/research/figures/surfaceannotations/fsaverageflatannotations.png');  % uint8 0-255 where 255 is the labeling
[~,~,alpha] = imread('~/Dropbox/KKTEMP/fsaverageflatannotations.png');  % uint8 0-255 where 255 is the labeling
alpha = double(alpha)/255;  % 0-1 where 1 is the labeling

% do it
for subjix=1:8
  for p=1:nsess(subjix)

    % R x C, 0-1 (1 means opaque)
    textalpha = drawtexts(1000,0,0,'FreeMono',0.05,[1 1 1],[0 0 0],sprintf('subj%02d, session %2d',subjix,p),{'FontWeight','bold'});
    % R x C (cropped)    [write a function to crop tight???]
    if subjix==1 && p==1
      colix = find(sum(textalpha,1)~=0);  % which columns
      rowix = find(sum(textalpha,2)~=0);  % which rows
    end
    textalpha = textalpha(:,colix(1):colix(end));
    textalpha = textalpha(rowix(1):rowix(end),:);
    % create solid color, R x C x 3
    textalphasolid = repmat(reshape(double([255 255 255]),1,1,3),[size(textalpha,1) size(textalpha,2)]);

    for tt=1:length(todo)

      % define
      file0 = sprintf('%s_subj%02d_sess%02d.png',todo{tt},subjix,p);
      outfile0 = sprintf('%s_subj%02d_sess%02d_processed.png',todo{tt},subjix,p);

      % load image
      im = imread(file0);

      % create solid color, R x C x 3
      alphasolid = repmat(reshape(colors{tt},1,1,3),[size(alpha,1) size(alpha,2)]);

      % mix image with alphasolid
      im = uint8(double(im).*(1-alpha) + alphasolid.*alpha);

      % mix image with textalphasolid     [write function to deal with this placement???]
      r1 = round(170-size(textalpha,1)/2); c1 = round(size(im,2)/2 - size(textalpha,2)/2);
      im(r1-1+(1:size(textalpha,1)),c1-1+(1:size(textalpha,2)),:) = ...
        double(im(r1-1+(1:size(textalpha,1)),c1-1+(1:size(textalpha,2)),:)).*(1-textalpha) + textalphasolid.*textalpha;

      % write image
      imwrite(im,outfile0);

    end
    
  end
end

%% %%%%% Generate movie

ffmpeg -framerate 30 -pattern_type glob -i 'meansession_*_processed.png' -crf 18 -c:v libx264 -pix_fmt yuv420p grandmeansurface.mp4
ffmpeg -framerate 30 -pattern_type glob -i 'R2session_*_processed.png'   -crf 18 -c:v libx264 -pix_fmt yuv420p grandR2surface.mp4

superhq 1080 in handbrake
