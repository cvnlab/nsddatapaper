%% This is just run one-time to generate some files.

*** note:
intereseting stimulus statistics:
model wise cross correlation...
  contrast,
  foreground, etc.

%%%

% define
types = {'Faces' 'Words' [] [] 'Bodies' 'Foreground'};
labeltypes = {'Automated' 'Semi-automated'};
imageres = 50;  % resolution for the prf model (8.4/50 = 0.168 deg)
modelnames = {'category' 'prf' 'fullprf'};
manualchoices = {[1 2 5 6] [1 2 4 5 6]};

%%%%%% STIMULUS PRE-PROCESSING (ONE-TIME) [BOUNDING BOXES]

% This was run to create:
% - versionA. manualchoices = {[1 2 5 6] [1 2 4 5 6]};

% define
modeltype = 2;

% loop over category type
for p=1:length(types)

  % loop over annotation type
  for q=1:length(labeltypes)

    % load annotation information
    switch modeltype
    case 1
      % no need to pre-compute

    case {2 3}
    
      % load
      switch q
      case 1
        if ismember(p,[1 2])
          a1 = loadmulti(sprintf('~/nsd/nsdextensions/NSD_Annotation_Efforts_1.0/%s/%s/binaryimage_mode2.mat',labeltypes{q},types{p}),'binaryimagemat');
            % binaryimagemat: [425x425x73000 int8] with number of bounding boxes at each location
        elseif p==5
          a1 =  load(sprintf('~/nsd/nsdextensions/BodyAnnotationsCOCO/animals_annotation_73k.mat'));
            % a1 = 
            % 
            %   struct with fields:
            % 
            %             bmap: [73000x425x425 int8]
            %          outcell: {73000x4 cell}
            %     output_count: [73000x1 int8]
          a1b = load(sprintf('~/nsd/nsdextensions/BodyAnnotationsCOCO/persons_annotation_73k.mat'));
            % a1b = 
            % 
            %   struct with fields:
            % 
            %             bmap: [73000x425x425 int8]
            %          outcell: {73000x4 cell}
            %     output_count: [73000x1 int8]
          a1 = permute(a1.bmap,[2 3 1]) + permute(a1b.bmap,[2 3 1]);
          clear a1b;
        elseif p==6
          a1 = load('~/nsd/nsdextensions/ForegroundMapCOCO/foreground_map_73k.mat');
            % a1 = 
            % 
            %   struct with fields:
            % 
            %     bmap: [73000x425x425 int8]
          a1 = permute(a1.bmap,[2 3 1]);
        end
      case 2
        assert(ismember(p,[1 2]));
        a1 = loadmulti(sprintf('~/nsd/nsdextensions/NSD_Annotation_Efforts_1.0/%s/%s/Combined/Outputbigmat_mode2v7.mat',labeltypes{q},types{p}),'outputbmat');
          % outputbmat: [425x425x1000x9 uint8]. For example, outputbmat(i,j,k,l) is the
          %   number of l faces in the kth image at the (i,j) position, averaged across raters, and then multiplied by 10.
          % unique(a1.outputbmat(:)) is [0 5 10 15 20].
      end
      
      % downsample
      stim = zeros(imageres,imageres,size(a1,3),'single');  % NOTE: single
      for ii=1:size(a1,3)
        if q==1
          temp = double(double(a1(:,:,ii))>0);
        else
          temp = double(any(double(a1(:,:,ii,manualchoices{p}))>0,4));
        end
        stim(:,:,ii) = imresize(temp,[imageres imageres],'cubic');
      end
      stim = normalizerange(stim,0,1,0,1);  % force to range 0-1
      
      % clean up
      clear a1;
      
      % save
      save(sprintf('~/ext/figurefiles/nsd/stim%s_p%d_q%d_versionA.mat',modelnames{modeltype},p,q),'stim','-v7.3','-nocompression');

    end
  end
end

%%%%%% STIMULUS PRE-PROCESSING (ONE-TIME) [SALIENCE]

%%% OOPS. THIS IS TRANSPOSED. CORRECTION IS mod(270-ang,360)

% load
ims = h5read('/home/surly-raid4/kendrick-data/nsd/nsdextensions/salience/salience_maps.hdf5','/img');  % uint8

% process
stim = zeros(imageres,imageres,size(ims,3),'single');  % NOTE: single
for p=1:size(ims,3)
  stim(:,:,p) = imresize(double(ims(:,:,p)),[imageres imageres],'cubic');
end
stim = normalizerange(stim,0,1,0,255);  % force to range 0-1

% save
save('~/ext/figurefiles/nsd/stimsalience.mat','stim','-v7.3','-nocompression');

%%%%%% STIMULUS PRE-PROCESSING (ONE-TIME) [CONTRASTGRID]

% This was run to create:
% - versionA.

% load and prepare images
stim = zeros(imageres+1,imageres+1,73000,'single');  % 50 x 50 x 73000 (ALT: 51 x 51 x 73000)
npx = 800/imageres;
for p=1:73
  statusdots(p,73);
  ims = h5read(nsdstimfile,'/imgBrick',[1 1 1 (p-1)*1000+1],[3 425 425 1000]);
  for q=1:1000
    im0 = permute(ims(:,:,:,q),[3 2 1]);
    im0 = imresize(single(rgb2gray(im0)),[800 800],'cubic');
    im0(im0<0) = 0;
    im0 = (im0/255).^2;                            % convert to [0,1] and square to match display gamma
    im0 = placematrix((127/255)^2*ones(800+npx,800+npx),im0,[]);  % embed within gray
    for rowix=1:imageres+1
      for colix=1:imageres+1
        rowii = (rowix-1)*npx + (1:npx);
        colii = (colix-1)*npx + (1:npx);
        stim(rowix,colix,(p-1)*1000+q) = std(flatten(im0(rowii,colii)));  % standard deviation of pixels
      end
    end
  end
end
clear ims;

% max is .502247

% deal with units
stim = 2*stim;  % in the limit, [0 1] produces std of 0.5. So we multiply by 2 so that maximum stimulation is 1
stim = normalizerange(stim,0,1,0,1);

% save
save(sprintf('~/ext/figurefiles/nsd/stim%s_versionA.mat','contrastgridNEW'),'stim','-v7.3','-nocompression');
  % contrastgrid: old version ignoring the background
  % contrastgridNEW: new version offset by 1
  