function autoqc_fmri(datadir,expecteddim)

ppdir = regexprep(datadir,'rawdata','ppdata');
mkdirquiet(ppdir);

% init
pptc = {};    % cell type x run, elements are pp x volumes  (percentile of raw data)
histall = {}; % cell type x run, elements are 1 x bins with counts
mmtc = {};    % cell 1 x type, elements are pp x runs  (percentile of averaged mean-abs-diff)
mnvol = {};   % cell 1 x type, element is the mean of the mean volume
tsnr = [];    % [tsnrtype1/tsnrtype2] x runs with median tSNR
tsnrvol = []; % X x Y x Z x (tsnrtype1/tsnrtype2)

% define
types = {'*bold*' '*field_map*' '*field_map*'};  % how to match EPI and the FIELDMAP; EPI, phase, mag

% load
for tt=1:length(types)

  % calc
  switch tt
  case 1
    argix = 1;
  case 2
    argix = 2;
  case 3
    argix = 2;
  end

  % get the directories
  files0 = matchfiles(sprintf('%s/dicom/%s',datadir,types{tt}));
  
  % deal with field map case
  switch tt
  case 1
  case 2  % phase
    files0 = files0(2:2:end);
  case 3  % mag (and remember to take only the first half of the slices)
    files0 = files0(1:2:end);
  end
  
  % check number of directories
  assert(isequal(length(files0),length(expecteddim{argix}{4})));

  % do it
  for pp=1:length(files0)
  
    % load
    [vols,volsizes,inplanematrixsizes,trs] = dicomloaddir(files0{pp});
    
    % drop some in the case of fieldmap mag
    if tt==3
      vols{1} = vols{1}(:,:,1:end/2);
    end
    
    % checks
    assert(isequal(round(1000*volsizes{1}),1000*expecteddim{argix}{1}));  % voxel size
    assert(isequal(trs{1},expecteddim{argix}{2}));                        % TR
    assert(isequal(sizefull(vols{1},3),expecteddim{argix}{3}));           % matrix dimensionality
    assert(isequal(size(vols{1},4),expecteddim{argix}{4}(pp)));           % number of volumes
    
    % percentiles of raw data over time
    pptc{tt,pp} = single(prctile(squish(vols{1},3),[5 25 50 75 95 100],1));
    
    % histograms of all of the raw data
    bb = 0:15:4095;
    histall{tt,pp} = hist(vols{1}(:),bb);

    % mean volume
    mnvol{tt}(:,:,:,pp) = mean(vols{1},4);
    valid = mnvol{tt}(:,:,:,pp) > 0.1 * prctile(flatten(mnvol{tt}(:,:,:,pp)),99);

    % percentiles of averaged mean-abs-diff volume
    temp1 = mean(abs(diff(vols{1},[],4)),4);
    mmtc{tt}(:,pp) = single(prctile(squish(temp1,3),[5 25 50 75 95],1));
    
    % tsnr
    if tt==1

      % type 1 [quadratic detrended, standard tSNR metric]
      temp2 = mnvol{tt}(:,:,:,pp) ./ ...
                reshape(std(projectionmatrix(constructpolynomialmatrix(size(vols{1},4),0:2))*squish(vols{1},3)',[],1),sizefull(mnvol{tt},3));
      tsnr(1,pp) = median(temp2(valid));
      
      % type 2 [mad-based tSNR]
      temp3 = mnvol{tt}(:,:,:,pp) ./ temp1;
      tsnr(2,pp) = median(temp3(valid));
      
      % save tsnr volumes for the first run
      if pp==1
        tsnrvol = single(cat(4,temp2,temp3));
      end

    end

  end
  mnvol{tt} = single(mean(mnvol{tt},4));

end

% save
mkdirquiet(ppdir);
save([ppdir '/autoqc_fmri.mat'],'pptc','histall','mmtc','mnvol','tsnr','tsnrvol');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK:

%%ssvol = {};   % cell 1 x type, element is chunk the mean volume and then average across runs
%  % internal constant
%tchunk = 5;


%     % homogeneity
%     chunks = {};
%     for dim=1:3
%       chunksize = round(size(vols{1},dim)/tchunk);
%       chunks{dim} = chunking(1:size(vols{1},dim),chunksize);
%     end
%     for xx=1:tchunk
%       for yy=1:tchunk
%         for zz=1:tchunk
%           ssvol{tt}(xx,yy,zz,pp) = mean(flatten(mnvol{tt}(chunks{1}{xx},chunks{2}{yy},chunks{3}{zz},pp)));
%         end
%       end
%     end
%   ssvol{tt} = mean(ssvol{tt},4);
  % ,'ssvol');
