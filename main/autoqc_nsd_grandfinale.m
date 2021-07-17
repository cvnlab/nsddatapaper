function autoqc_nsd_grandfinale(glmdirname,curix)

% setup
nsdsetup;

alltempfiles = {};

switch glmdirname

case {'glmdata' 'glmdataLOW'}

  if isequal(glmdirname,'glmdata')
    min0 = -620;
  else
    min0 = -340;
  end

  % what is the 1-index of the max nsd session we have?
  if ~exist('curix','var') || isempty(curix)
    curix = lastel(find(any(sessmatrix~=0,1)));  % actually, this needs to be +2 when we get the aux experiments
  end

  % calc
  outputdir = sprintf('%s/ppdata/NSD-%s-glmBASIC',nsddir,glmdirname);
  dirs0 = matchfiles(sprintf('%s/ppdata/NSD???-%s-glmBASIC',nsddir,glmdirname));

  % make output dir
  mkdirquiet(outputdir);

  % orthMID
  masterfile = {};
  for p=0:curix  % START FROM SCREENING
    order = {[1 5] [2 6] [3 7] [4 8]};
    for q=1:8
      a0 = matchfiles(sprintf('%s/orthMID%02d.png',dirs0{q},p));
      if ~isempty(a0)
        masterfile{q} = a0{1};

        % SPECIAL TITLES!!!
        imtemp = imread(masterfile{q});
        imtemp = cat(1,255*ones(1,size(imtemp,2)),imtemp);
        imtemp = cat(2,255*ones(size(imtemp,1),1),imtemp);
        im = drawtexts(500,0,0,'FreeMono',0.04,[1 1 1],[0 0 0],sprintf('subj%02d, session %2d',q,p),{'FontWeight','bold'});
        im = uint8(255*placematrix(zeros(30,size(imtemp,2)),im));
        imtemp = cat(1,im,imtemp);
        fileX = [tempname '.png'];
        imwrite(imtemp,fileX);
        masterfile{q} = fileX;
        alltempfiles = [alltempfiles {fileX}];
      end
    end
    im0 = concatimages(masterfile,order,0,min0);
    imwrite(uint8(im0),gray(256),sprintf('%s/orthMID%03d.png',outputdir,p));
  end

  % varexporthMID
  masterfile = {};
  for p=1:curix  % START FROM NSD01
    order = {[1 5] [2 6] [3 7] [4 8]};
    for q=1:8
      a0 = matchfiles(sprintf('%s/varexporthMID%02d.png',dirs0{q},p));
      if ~isempty(a0)
        masterfile{q} = a0{1};

        % SPECIAL TITLES!!!
        imtemp = imread(masterfile{q});
        imtemp = cat(1,255*ones(1,size(imtemp,2)),imtemp);
        imtemp = cat(2,255*ones(size(imtemp,1),1),imtemp);
        imtemp = cmaplookup(imtemp,-.5,255.5,0,hot(256));  % RGB
        im = drawtexts(500,0,0,'FreeMono',0.04,[1 1 1],[0 0 0],sprintf('subj%02d, session %2d',q,p),{'FontWeight','bold'});
        im = uint8(255*placematrix(zeros(30,size(imtemp,2)),im));
        im = cmaplookup(im,-.5,255.5,0,gray(256));  % RGB
        imtemp = cat(1,im,imtemp);
        fileX = [tempname '.png'];
        imwrite(imtemp,fileX);
        masterfile{q} = fileX;
        alltempfiles = [alltempfiles {fileX}];
      end
    end
    im0 = concatimages(masterfile,order,0,min0);
    imwrite(uint8(im0),sprintf('%s/varexporthMID%03d.png',outputdir,p));  % not indexed, RGB!
  end

case {'structurals'}

  % calc
  outputdir = sprintf('%s/ppdata/NSD-structurals',nsddir);
  dirs0 = matchfiles(sprintf('%s/ppdata/NSD???-structurals',nsddir));

  % make output dir
  mkdirquiet(outputdir);

  % define
  typs = {'T1' 'T2' 'EPIsyn' 'SWI' 'TOF'};
  namelabels = {'T1 ' 'T2 ' 'EPI' 'SWI' 'TOF'};
  
  % loop
  masterfile = {};
  for p=1:length(typs)
    order = {[1 5] [2 6] [3 7] [4 8]};
    for q=1:8
      a0 = matchfiles(sprintf('%s/inspection/%s.png',dirs0{q},typs{p}));
      out0 = sprintf('%s/inspection/%s_processed.png',dirs0{q},typs{p});
      if ~isempty(a0)

        %%% [hmm.. can some part of this be made into a function??]
        imtemp = imread(a0{1});
        imtemp = cmaplookup(imtemp,-.5,255.5,0,gray(256));  % RGB, values up to [1 1 1]
        im = drawtexts(500,0,0,'FreeMono',0.04,[1 1 1],[0 0 0],sprintf('subj%02d, %s',q,namelabels{p}),{'FontWeight','bold'});
        alphaim = placematrix(zeros(24,150),mean(im,3));  % 1 means show text
        r1 = 1; c1 = 1;
        imtemp(r1-1+(1:size(alphaim,1)),c1-1+(1:size(alphaim,2)),:) = ...
          double(imtemp(r1-1+(1:size(alphaim,1)),c1-1+(1:size(alphaim,2)),:)).*repmat(1-alphaim,[1 1 3]) + ...
          repmat(reshape([1 1 1],1,1,3),[size(alphaim,1) size(alphaim,2)]).*repmat(alphaim,[1 1 3]);
        imwrite(uint8(255*imtemp),gray(256),out0);
        
        masterfile{q} = out0;
      end
    end
    im0 = concatimages(masterfile,order);
    imwrite(uint8(im0),gray(256),sprintf('%s/%s.png',outputdir,typs{p}));
  end

end

for zz=1:length(alltempfiles)
  delete(alltempfiles{zz});
end
