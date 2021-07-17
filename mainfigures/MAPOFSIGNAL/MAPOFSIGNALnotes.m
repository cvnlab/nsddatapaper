%%%%%%%%%%%%%%% postprocess the bad stuff into outlines

todo = {'subj01' 'subj02' 'subj03' 'subj04' 'subj05' 'subj06' 'subj07' 'subj08' 'groupavg'};
things = {'signaldropout' 'surfaceimperfections'};
threshs = [4/8-1/16 1/8-1/16];
cmaps = {winter(256) winter(256)};
colors = {hsv2rgb([30/360 100/100 90/100]) hsv2rgb([120/360 100/100 90/100])};
for yy=1:length(things)
  for zz=1:length(todo)
    infile = sprintf('~/Dropbox/nsddata/inspections/surfacevisualizations/fsaverageflat_%s_%s.png',things{yy},todo{zz});
    outfile = sprintf('/research/papers/2022 nsd/figures/MAPOFSIGNAL/figures/fsaverageflat_%s_%s_outline.png',things{yy},todo{zz});
    
    im = double(imread(infile));
    
    xx = linspace(0,1,257);
    xx = (xx(1:end-1)+xx(2:end))/2;
    cmap0 = double(uint8(255*cmaps{yy}));
    cmap0 = cmap0(:,1)*256^2 + cmap0(:,2)*256 + cmap0(:,3);

    imtemp = im(:,:,1)*256^2 + im(:,:,2)*256 + im(:,:,3);
    vals = NaN*zeros(size(imtemp));
    for p=1:256
      vals(imtemp==cmap0(p)) = xx(p);
    end

    im = detectedges(double(vals>=threshs(yy)),1.5);
    alpha = im / max(im(:));
    baseim = repmat(reshape(colors{yy},1,1,3),[size(alpha,1) size(alpha,2)]);
    imwrite(uint8(255*baseim),outfile,'Alpha',alpha);
  end
end

%%%%%%%%%%%%%%% quantitative plot

%% LOAD

todo = { ...
  {'ncsnr'     ['~/nsd/nsddata_betas/ppdata/subj%02d/nativesurface/betas_assumehrf/%s.%s.mgh']             @(x)100*(x.^2./(x.^2+1/3))         [0 75]        jet(256)     []          'b1nc'} ...
  {'ncsnr'     ['~/nsd/nsddata_betas/ppdata/subj%02d/nativesurface/betas_fithrf/%s.%s.mgh']                @(x)100*(x.^2./(x.^2+1/3))         [0 75]        jet(256)     []          'b2nc'} ...
  {'ncsnr'     ['~/nsd/nsddata_betas/ppdata/subj%02d/nativesurface/betas_fithrf_GLMdenoise_RR/%s.%s.mgh']  @(x)100*(x.^2./(x.^2+1/3))         [0 75]        jet(256)     []          'b3nc'} ...
};
hemis = {'lh' 'rh'};
extensions = {'' '_split1' '_split2'};

vals = [];  % vertices x subj x hemis x betatype x full/test/retest
for tt=1:length(todo)
  name =    todo{tt}{1};
  fileloc = todo{tt}{2};
  tfun =    todo{tt}{3};
  rng =     todo{tt}{4};
  cmap =    todo{tt}{5};
  thresh =  todo{tt}{6};
  outname = todo{tt}{7};

  % load data and map to fsaverage
  for p=1:8
    for hh=1:length(hemis)
      for ee=1:length(extensions)
        inputfile = matchfiles(sprintf(fileloc,p,hemis{hh},[name extensions{ee}]));
        assert(length(inputfile)==1);
        inputfile = inputfile{1};
        data = cvnloadmgz(inputfile);
        vals(:,p,hh,tt,ee) = tfun(nsd_mapdata(p,[hemis{hh} '.white'],'fsaverage',data));
      end
    end
  end
end

% load ROI
roimask1 = cvnroimask('fsaverage','lh','nsdgeneral');
roimask2 = cvnroimask('fsaverage','rh','nsdgeneral');
roimask = cat(1,roimask1{1},roimask2{1});

%% VISUALIZE

% big histograms
bins = .25:.5:100;
%bins = .5:1:99.5;
%bins = 1:2:99;
comp = {[1 2] [2 3] [1 3]};
figureprep([100 100 500 900]); hold on;
for cc=1:3
  for ss=1:8
    subplot(8,3,(ss-1)*3+cc); hold on;
    [n,x,y] = hist2d(flatten(vals(:,ss,:,comp{cc}(1),1)),flatten(vals(:,ss,:,comp{cc}(2),1)),bins);
    imagesc(x(1,:),y(:,1),log(1+n));
    colormap(flipud(bone(256)));
    set(gca,'YDir','normal');
    axissquarify([],'r-');
    axis([0 100 0 100]);
    straightline([25 50 75],'h','k:');
    straightline([25 50 75],'v','k:');
    if ss<=7
      set(gca,'XTick',[]);
    else
      set(gca,'XTick',0:25:100);
    end
    if cc>=2
      set(gca,'YTick',[]);
    else
      set(gca,'YTick',0:25:100);
    end
    if ss==8
      xlabel(sprintf('Noise ceiling for b%d',comp{cc}(1)));
      ylabel(sprintf('Noise ceiling for b%d',comp{cc}(2)));
    end
  end
end
figurewrite('histograms',[],-2,'~/Dropbox/KKTEMP/');

% reliability
test1 = subscript(squish(permute(vals(:,:,:,:,2),[1 3 2 4]),2),{logical(roimask) ':' ':'});  % vertices x subj x betatype
test2 = subscript(squish(permute(vals(:,:,:,:,3),[1 3 2 4]),2),{logical(roimask) ':' ':'});  % vertices x subj x betatype
reliabilities = calccorrelation(test1(:,:,3),test2(:,:,3),1)
min(reliabilities)
max(reliabilities)
%   Columns 1 through 4
% 
%          0.987465388249121         0.989233269883145         0.970685416300537         0.966845667856931
% 
%   Columns 5 through 8
% 
%          0.986597788023607         0.969685776047881         0.970748514556372         0.946420144787961
%
% 
%          0.946420144787961
% 
%          0.989233269883145
bins = .25:.5:100;
figureprep([100 100 1200 500]);
for p=1:8
  subplot(1,8,p); hold on;
%  title(sprintf('subject %d (r = %.2f)',p,reliabilities(p)));
  title(sprintf('r = %.2f',reliabilities(p)));
  [n,x,y] = hist2d(test1(:,p,3),test2(:,p,3),bins,bins);
  imagesc(x(1,:),y(:,1),log(1+n));
  colormap(flipud(bone(256)));
  set(gca,'YDir','normal');
  axissquarify([],'r-');
  axis([0 100 0 100]);
  straightline([25 50 75],'h','k:');
  straightline([25 50 75],'v','k:');
  if p==1
    set(gca,'XTick',0:25:100);
    set(gca,'YTick',0:25:100);
    xlabel('b3 noise ceiling (%), split 1');
    ylabel('b3 noise ceiling (%), split 2');
  else
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
  end
  %scatter(test1(:,p,3),test2(:,p,3),'r.');
end
figurewrite('reliability',[],-2,'~/Dropbox/KKTEMP/');

% summary
summary = squish(median(subscript(squish(permute(vals(:,:,:,:,1),[1 3 2 4]),2),{logical(roimask) ':' ':'}),1),2);
figureprep([100 100 400 120]); hold on;
bar(summary,1);
ylabel('Noise ceiling (%)');
ax = axis;
axis([0 9 0 50]);
set(gca,'XTick',1:8);
figurewrite('summary',[],-2,'~/Dropbox/KKTEMP/');
save('~/Dropbox/KKTEMP/results.mat','summary');

% p-values
[h,pval1] = ttest(summary(:,1),summary(:,2))
[h,pval2] = ttest(summary(:,2),summary(:,3))
[h,pval3] = ttest(summary(:,1),summary(:,3))
% pval1 =
%       0.000405092581134296
% pval2 =
%       7.57913779013044e-07
% pval3 =
%        6.2886688028074e-06

% Compute equivalent SNR1 trials!
ncsnrtemp = sqrt((summary/100)*(1/3) ./ (1-summary/100))           % median noise ceiling, converted back to ncsnr
totalnsdsessions = [40 40 32 30 40 32 40 30];
etrials = repmat((totalnsdsessions*750)',[1 3]) .* ncsnrtemp.^2    % number of trials * ncsnr^2
sum(etrials,1)
% ncsnrtemp =
% 
%          0.358246803282803         0.436626255377687         0.490176752072449
%           0.33599907159801         0.427050679599555         0.510127365586935
%          0.275612577796804         0.360503435134597         0.430892258881179
%           0.28948661685029         0.328749939799874         0.392465204000473
%          0.338273942469937           0.4550386518031         0.542167603966064
%            0.3004481941528         0.353647351264954         0.438968554129977
%          0.262474209070296         0.313352897763734         0.370504736878823
%          0.214662194253004          0.24070001393563         0.322861254215315
% 
% 
% etrials =
% 
%           3850.22316187043          5719.27460655423          7208.19744816886
%           3386.86128344174          5471.16848839326          7806.89787361999
%           1823.09503295519          3119.10544185226           4456.0353303294
%            1885.5562800471          2431.72176566447          3465.65106790049
%           3432.87780462462          6211.80523904348           8818.3713237091
%           2166.45881687228          3001.59477736122          4624.64139635911
%           2066.78131281232          2945.70115610787          4118.21280148938
%           1036.79679693407          1303.57117594378          2345.38626315345
% 
% 
% ans =
% 
%           19648.6504895577          30203.9426509206          42843.3935047298

% one trial at SNR==1 (SNR1)
% if SNR==2, that's like 4 SNR1 trials
% you can use your SNR1 trials to either do reps or fresh images

% common currency is, what is the (median) ncsnr? for NSD and BOLD5000

% what is noise ceiling median, mean across subjects?
mean(summary(:,3))
sqrt(mean(summary(:,3))/100)
% 36.1991070241837
% 0.601656937333758

% what is per-subject equivalent number of trials?
mean(etrials,1)
%          2456.08131119472          3775.49283136507          5355.42418809122

% 8 subjects, 5355 trials each, total of 42843 trials
% 3.5 subjects, ??BOLD5000?

%%% UPDATES (with Jacob's additions)

% We have 213000 stimulus trials in total across 8 NSD subjects
% In BOLD5000, there are 5254*3 + 3108 = 18870 trials.

%%%%%%%%%%%%%%%

figureprep([100 100 470 200]);
drawcolorbar([0 75],0:25:75,jet(256),'Noise ceiling (%)',1);
figurewrite('nccolorbar',[],-2,'~/Dropbox/KKTEMP/');

%%%%%%%%%%%%%%%

notes:
- massive montages of all individual subjects? started an illustrator document but... hold off...

%%%%%%%%%%%%%%%%%%%%%%%% JUNK BELOW:

temp = im-rgb2gray(im);
bad = std(temp,[],3)==0;

std(im-rgb2gray(im),[],3);
>> figure;imagesc(temp)

  im = detectedges(std(im-rgb2gray(im),[],3),0.5);
  alpha = im / max(im(:));
  baseim = repmat(reshape([1 0 0],1,1,3),[size(alpha,1) size(alpha,2)]);
  imwrite(uint8(255*baseim),'test.png','Alpha',alpha);
