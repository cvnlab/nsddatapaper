% define
sessix = 20;
runix = 1;
sliceix = 42;
epislice = [1 14 27 12 25 10 23 8 21 6 19 4 17 2 15 28 13 26 11 24 9 22 7 20 5 18 3 16 1 14 27 12 25 10 23 8 21 6 19 4 17 2 15 28 13 26 11 24 9 22 7 20 5 18 3 16 1 14 27 12 25 10 23 8 21 6 19 4 17 2 15 28 13 26 11 24 9 22 7 20 5 18 3 16];
offsets = -10:10:0;

% load data
dataBALL = {};
onsetsALL = {};
for subjix=1:8

  % load raw fMRI data
  file0 = sprintf('/home/surly-raid3/kendrick-data/nsd/nsddata_rawdata/data/sub-%02d/ses-nsd%02d/func/*nsdcore*run-%02d_bold.nii.gz',subjix,sessix,runix);
  file0 = matchfiles(file0);
  data = load_untouch_nii(file0{1});

  % extract axial slice
  dataB = permute(double(data.img(:,:,sliceix,:)),[1 2 4 3]);  % R x C x TR
  dataBALL{subjix} = flipdim(rotatematrix(dataB,1,2,1),2);

  % load experiment timing
  file0 = sprintf('/home/surly-raid3/kendrick-data/nsd/nsddata_rawdata/data/sub-%02d/ses-nsd%02d/func/*nsdcore*run-%02d_events.tsv',subjix,sessix,runix);
  file0 = matchfiles(file0);
  event0 = importdata(file0{1});
  onsetsALL{subjix} = cellfun(@str2double,event0.textdata(2:end,1));

end

% loop over subjects
figureprep([0 0 750 1000]);
for subjix=1:8

  % load
  dataB = dataBALL{subjix};
  onsets = onsetsALL{subjix};
  
  % calculate the data's offset from 0 s
  offset0 = (epislice(sliceix)-1)*(1.6/max(epislice));

  % interpolate data to 0:1:... in order to fit the simple model
  dataC = tseriesinterp(dataB,1.6,1,3,[],offset0);
  dataC = dataC(:,:,1:300);

  % compute tsnr
  tsnr = mean(dataB,3) ./ ...
            reshape(std(projectionmatrix(constructpolynomialmatrix(size(dataB,3),0:2))*squish(dataB,2)',[],1),sizefull(dataB,2));

  % construct simple ON-OFF design matrix
  design = zeros(300,1);
  for q=1:length(onsets)
    design(round(onsets(q)+1),1) = 1;
  end

  % fit simple ONOFF GLM model
  results = GLMestimatemodel(design,squish(dataC,2),3,1,'assume',[],0);
  R2 = reshapesquare(results.R2);
  [mx,ix] = max(R2(:));

  % location of max ONOFF R2
  [i1,i2] = ind2sub([120 120],ix);

  % visualize raw slice
  subplot(8,4,(mod2(subjix,2)-1)*2 + (ceil(subjix/2)-1)*8 + [1]); hold on;
  im0 = cmaplookup(dataB(:,:,1),0,2000,[],gray(256));
  imagesc(im0); axis image; set(gca,'YDir','reverse'); axis off;
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  for pp=1:length(offsets)
    scatter(i2,i1+offsets(pp),'c.');
  end

  % visualize R2
  subplot(8,4,(mod2(subjix,2)-1)*2 + (ceil(subjix/2)-1)*8 + [5]); hold on;
  im0 = cmaplookup(R2,0,40,[],hot(256));
  imagesc(im0); axis image; set(gca,'YDir','reverse'); axis off;
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  for pp=1:length(offsets)
    scatter(i2,i1+offsets(pp),'c.');
  end

  % visualize time-series data
  for pp=1:length(offsets)
  
    alllocs = (mod2(subjix,2)-1)*2 + (ceil(subjix/2)-1)*8 + [2 6];
    subplot(8,4,alllocs(pp)); hold on;
  
    dataD = flatten(dataB(i1+offsets(pp),i2,:));  % 1 x TR
    dataE = flatten(dataC(i1+offsets(pp),i2,:));  % 1 x TR
  
    h = plot(offset0 + linspacefixeddiff(0,1.6,length(dataD)),dataD,'r-','LineWidth',1);
    ax = axis;
    set(straightline(onsets,'v','k-'),'Color',[.75 .75 .75],'LineWidth',0.5);
    uistack(h,'top');

    hrf = getcanonicalhrf(3,1);
    pred = conv(design,hrf);
    pred = pred(1:300);
    X = cat(2,pred(:),constructpolynomialmatrix(300,0:3));
    h = X\dataE(:);
    pred2 = X*h;
    plot(0:length(pred2)-1,pred2,'k-','LineWidth',1);
    %xlabel('Time (s)');
    %ylabel('Signal intensity (raw scanner units)');
    if pp==length(offsets)
    else
      set(gca,'XTickLabel',[]);
    end
    axis([0 300 ax(3:4)]);
  
    title(sprintf('ON-OFF R^2 = %.0f%%, tSNR = %.0f',R2(i1+offsets(pp),i2),tsnr(i1+offsets(pp),i2)),'FontSize',9,'FontWeight','normal');

  end

end
figurewrite('rawinspection',[],-2,'~/Dropbox/KKTEMP/');
